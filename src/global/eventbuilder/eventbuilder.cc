// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <JANA/JApplication.h>
#include <JANA/JEventUnfolder.h>
#include <JANA/JException.h>
#include <JANA/JLogger.h>

#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/CaloHitContribution.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociation.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/CalorimeterHitCollection.h>
#include <edm4eic/MCRecoCalorimeterHitAssociationCollection.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/eventbuilder/TimeAlignment_factory.h"
#include "factories/eventbuilder/TimeCoincidence_factory.h"
#include "factories/eventbuilder/TrkCoincidence_factory.h"
#include "factories/eventbuilder/CalCoincidence_factory.h"
#include "factories/eventbuilder/EventBuilder_factory.h"
#include "factories/eventbuilder/EventPrefilter_factory.h"
#include "services/eventbuilder/TimesliceBuffer_service.h"
#include "services/onnx/ONNXRuntime_service.h"


// ---------------------------------------------------------------------------
// EventUnfolder — minimal JANA unfolder (no analysis).
//
// Reads the candidate-event list from EventBuilder_factory ("EventCandidates",
// an edm4hep::EventHeaderCollection carrying t0/t0sigma) and materialises one
// PhysicsEvent per candidate. All event-finding / coincidence / topology logic
// lives in EventBuilder_factory; only a JEventUnfolder can create finer-level
// child events in JANA, so this thin glue is unavoidable.
// ---------------------------------------------------------------------------
struct EventUnfolder : public JEventUnfolder {

  // Acceptance-window controls. ALL detectors: half-width is computed per-hit
  // from edm4eic::TrackerHit::getTimeError(), which each detector's own
  // TrackerHitReconstruction_factory already populates from its digitization
  // config (10 ns Si, 8 ns B0, 10 ns MPGD, ~30 ps TOF) — the same sigma the
  // digitizer actually smears the time with, so the gate is exact by
  // construction. (A constant-override escape hatch for Si/B0 used to live
  // here, to emulate a real MAPS shutter/integration window the digitization
  // doesn't model -- removed as unused: nothing in the pipeline ever set it,
  // and the per-hit getTimeError() this defaulted to is the exact gate
  // anyway. Re-add if a systematic study of that scenario is ever needed.)
  // The aligned hit-time distribution of a real event is ASYMMETRIC: the
  // r/c correction assumes beta=1 straight lines, so massive/curling
  // particles arrive strictly LATE by (r/c)(1/beta-1) + arc-vs-chord — up to
  // a few ns at the outer radii (measured with `make gatecheck`, 2026-07).
  // A symmetric +-nsigma*sigma gate clips that tail and costs tracking
  // efficiency, worst for the 25 ps TOF. trk_late_nsigma_window_asym widens
  // ONLY the + side.
  //
  // late is a SIGMA MULTIPLE (of that hit's own getTimeError()), not a flat
  // ns constant: a flat ns add-on is not comparable across detectors of very
  // different resolution — 2 ns is ~0.2 sigma of a 10 ns Si/B0 hit but ~65
  // sigma of a 30 ps TOF hit, i.e. effectively unbounded for TOF while barely
  // widening Si/B0 at all. Scaling by the hit's own sigma keeps the extra
  // LATE allowance proportional for every detector the gate runs over.
  // Default below is an UNMEASURED PLACEHOLDER: the old flat-ns defaults
  // (trk=2.0 ns, cal=10.0 ns) don't translate to one nsigma across detectors
  // of such different resolution, so this needs its own efficiency/purity
  // retune (same `make gatecheck` / status --trigger --tracking --calo scan
  // the old defaults came from) once rebuilt.
  Parameter<float> trk_late_nsigma_window_asym{this, "eventbuilder:unfolder:trk_late_nsigma_window_asym", 0.5,
      "extra LATE acceptance, in units of that hit's/measurement's own sigma "
      "(getTimeError()), added to the + side of every per-hit and "
      "per-measurement gate (beta<1 arrival delay; 0 = legacy symmetric). "
      "UNMEASURED PLACEHOLDER — retune via an efficiency/purity scan."};

  Parameter<float> trk_nsigma_window{this, "eventbuilder:unfolder:trk_nsigma_window", 3.0,
                                 "N-sigma half-width for per-detector hit acceptance"};

  // Calo HITS (not clusters -- see below) get their OWN cal_nsigma_window/
  // cal_late_nsigma_window_asym rather than sharing the tracker/TOF ones
  // above: hit.getTimeError() only started being physically real in 2026-07
  // (previously hardcoded to 0 in CalorimeterHitReco.cc/
  // CalorimeterClusterRecoCoG.cc — see those files), and calo's timing
  // physics (shower development + light-collection spread) is a different
  // regime from tracker crossing-time + r/c propagation. Widening the SHARED
  // parameters to fit calo would have loosened the tracker/TOF gates too,
  // which already perform well at the defaults above.
  //
  // Calo CLUSTERS no longer use these at all: they're gated upstream by
  // CalCoincidence_factory (own nsigma/late_nsigma_window_asym under
  // eventbuilder:calcoincidence:*, plus a YR-based per-cluster energy floor)
  // before EventUnfolder ever sees them — see m_clu_in_names above. These two
  // parameters now apply only to the raw CalorimeterHit gate further down
  // (no coincidence factory exists for that type).
  //
  // The flat-ns defaults this replaced (nsigma=5/cal_late_ns=10 giving 29.1%
  // eff at 96.2% purity, see prior revision) were measured in absolute ns,
  // not against a sigma-proportional term, so cal_late_nsigma_window_asym's
  // default below is likewise an UNMEASURED PLACEHOLDER pending its own
  // retune — see trk_late_nsigma_window_asym's comment above for why a flat ns
  // number can't just be divided by "the" cluster sigma (there isn't a single
  // one across families).
  Parameter<float> cal_late_nsigma_window_asym{this, "eventbuilder:unfolder:cal_late_nsigma_window_asym", 0.5,
      "extra LATE acceptance, in units of that hit's own sigma "
      "(getTimeError()), added to the + side of the calo HIT gate only "
      "(shower-development/light-collection delay; separate from "
      "trk_late_nsigma_window_asym above, and from CalCoincidence_factory's "
      "cluster-level version). UNMEASURED PLACEHOLDER — retune via an "
      "efficiency/purity scan."};
  Parameter<float> cal_nsigma_window{this, "eventbuilder:unfolder:cal_nsigma_window", 5.0,
      "N-sigma half-width for the calo HIT gate (separate from trk_nsigma_window "
      "above, and from CalCoincidence_factory's cluster-level nsigma)"};

  // MCParticles carried into each child. Cloning the frame's FULL MCParticle
  // list per candidate multiplied file size by the candidate count (a 10-frame
  // dis_nc test wrote ~1 GB, almost all of it duplicated background MC).
  // Default: keep only particles produced within +-mc_time_window of the
  // candidate's t0 (the matched collision's products and coincident background;
  // fakes get almost nothing). Hit->MC truth associations are unaffected: they
  // reference the original Timeslice-level collection, not these clones.
  // Set to -1 to restore the old clone-everything behavior.
  Parameter<float> mc_time_window{this, "eventbuilder:unfolder:mc_time_window", 50.0,
      "half-width (ns) around t0 for MCParticles carried into the child; -1 = all"};

  // Cross-frame Si/B0 recovery window. Slow-detector hits of the adjacent
  // frames are fetched from the TimesliceBuffer_service (deposited there by
  // TrkTimeAlignment_factory) and gated per candidate like current-frame hits.
  // forward_frames adds no output latency: it only waits for the deposit of
  // frames that are already in flight (hence the max_inflight requirement).
  Parameter<int32_t> backward_frames{this, "eventbuilder:unfolder:backward_frames", 1,
      "Past frames whose Si/B0 hits are included when gating each candidate. 0 = off."};
  // Default 1: the measured hit-time residuals around t0 have a pronounced
  // LATE (Landau-like) tail — slow particles and delayed response arrive
  // after the collision — so spillover lands preferentially in the NEXT
  // frame. Forward recovery is therefore at least as important as backward.
  // Needs a concurrent pipeline (nthreads >= 2); when that is not available
  // checkConfig() degrades to 0 with a loud warning instead of aborting.
  Parameter<int32_t> forward_frames{this, "eventbuilder:unfolder:forward_frames", 1,
      "Future frames whose Si/B0 hits are included when gating each candidate. 0 = off. "
      "Requires nthreads >= 2 and jana:max_inflight_timeslices > forward_frames; "
      "degrades to 0 with a warning when unsupported."};

  // Diagnostic for the raw-hit truth carry-through: per candidate, prints
  // each detector's input/kept/output sizes for the trk + raw/link/asso
  // chains. The failure mode this exists for is SILENT (empty unfolder
  // outputs get no branch, and same-named child-level detector factories
  // then satisfy the writer with an ungated full-frame recompute instead --
  // see copyGatedByRawHit's comment), so nothing in the output file or the
  // default logs distinguishes "gated correctly" from "gate ate everything".
  Parameter<bool> debug_gating{this, "eventbuilder:unfolder:debug_gating", false,
      "Print per-candidate input/kept/output sizes for the raw-hit truth gating."};

  std::shared_ptr<TimesliceBuffer_service> m_ts_buffer;
  bool m_warned_candidate_overflow = false;
  bool m_checked_config = false;

  // Adjacent-frame slow-det hits, fetched once per parent (Unfold is sequential).
  uint64_t m_adjacent_parent = std::numeric_limits<uint64_t>::max();
  std::vector<TimesliceBuffer_service::FrameHits> m_adjacent;

  // -- tracker RecHits (10 detectors) ------------------------------------------
  // Fast detectors (TOF, MPGD): read from TrkCoincidence_factory output
  // (*TrkCoincRecHits, current frame) -- already candidate-t0 gated
  // (nsigma x per-hit getTimeError() + late asym, eventbuilder:trkcoincidence:*).
  // Slow detectors (Si, B0): read from TimeCoincidence_factory output
  // (*TimeCoincRecHits, current frame); adjacent-frame hits are merged in
  // Unfold() from the TimesliceBuffer_service (backward_frames/forward_frames).
  // Both groups used to be gated a second time, inline, right here with
  // trk_nsigma_window/trk_late_nsigma_window_asym -- now only Measurement2D
  // below still needs that (no coincidence factory exists for that type).
  std::vector<std::string> m_trk_in_names = {
      "TOFBarrelTrkCoincRecHits",          "TOFEndcapTrkCoincRecHits",
      "MPGDBarrelTrkCoincRecHits",         "OuterMPGDBarrelTrkCoincRecHits",
      "BackwardMPGDEndcapTrkCoincRecHits", "ForwardMPGDEndcapTrkCoincRecHits",
      "SiBarrelVertexTimeCoincRecHits",     "SiBarrelTrackerTimeCoincRecHits",
      "SiEndcapTrackerTimeCoincRecHits",    "B0TrackerTimeCoincRecHits"};
  std::vector<std::string> m_trk_out_names = {
      "TOFBarrelRecHits",       "TOFEndcapRecHits",          "MPGDBarrelRecHits",
      "OuterMPGDBarrelRecHits", "BackwardMPGDEndcapRecHits", "ForwardMPGDEndcapRecHits",
      "SiBarrelVertexRecHits",  "SiBarrelTrackerRecHits",    "SiEndcapTrackerRecHits",
      "B0TrackerRecHits"};

  // -- raw hits / links / hit associations (10 detectors each) -----------------
  // All detectors: current frame only. Cross-frame RecHit recovery is handled
  // upstream in TrkTimeAlignment_factory; raw hits carry no r/c correction
  // so cross-frame buffering at the raw level is not meaningful.
  //
  // TOF stays on the plain "RawHitFrame" chain ON PURPOSE, even though
  // BTOF.cc/ECTOF.cc also define a "Shared" (charge-sharing + pulse-shape)
  // frame chain: the gate below keys each raw hit against the getRawHit()
  // relation of THIS candidate's kept RecHits, and TOFBarrel/EndcapTimeAlign
  // RecHits derive from the plain chain -- so the raw/link/asso inputs must
  // be the SAME collections or no key ever matches (verified 2026-07-24: a
  // trial redirect of only these inputs to the Shared chain made every kept
  // TOF hit's key miss, silently zeroing the TOF truth output).
  std::vector<std::string> m_raw_in_names = {
      "TOFBarrelRawHitFrame",          "TOFEndcapRawHitFrame",
      "MPGDBarrelRawHitFrame",         "OuterMPGDBarrelRawHitFrame",
      "BackwardMPGDEndcapRawHitFrame", "ForwardMPGDEndcapRawHitFrame",
      "SiBarrelVertexRawHitFrame",     "SiBarrelRawHitFrame",
      "SiEndcapTrackerRawHitFrame",    "B0TrackerRawHitFrame"};
  std::vector<std::string> m_raw_out_names = {
      "TOFBarrelRawHits",       "TOFEndcapRawHits",          "MPGDBarrelRawHits",
      "OuterMPGDBarrelRawHits", "BackwardMPGDEndcapRawHits", "ForwardMPGDEndcapRawHits",
      "SiBarrelVertexRawHits",  "SiBarrelRawHits",           "SiEndcapTrackerRawHits",
      "B0TrackerRawHits"};

  std::vector<std::string> m_link_in_names = {
      "TOFBarrelRawHitLinkFrame",          "TOFEndcapRawHitLinkFrame",
      "MPGDBarrelRawHitLinkFrame",         "OuterMPGDBarrelRawHitLinkFrame",
      "BackwardMPGDEndcapRawHitLinkFrame", "ForwardMPGDEndcapRawHitLinkFrame",
      "SiBarrelVertexRawHitLinkFrame",     "SiBarrelRawHitLinkFrame",
      "SiEndcapTrackerRawHitLinkFrame",    "B0TrackerRawHitLinkFrame"};
  std::vector<std::string> m_link_out_names = {
      "TOFBarrelRawHitLinks",       "TOFEndcapRawHitLinks",          "MPGDBarrelRawHitLinks",
      "OuterMPGDBarrelRawHitLinks", "BackwardMPGDEndcapRawHitLinks", "ForwardMPGDEndcapRawHitLinks",
      "SiBarrelVertexRawHitLinks",  "SiBarrelRawHitLinks",           "SiEndcapTrackerRawHitLinks",
      "B0TrackerRawHitLinks"};

  std::vector<std::string> m_asso_in_names = {
      "TOFBarrelRawHitAssociationFrame",          "TOFEndcapRawHitAssociationFrame",
      "MPGDBarrelRawHitAssociationFrame",         "OuterMPGDBarrelRawHitAssociationFrame",
      "BackwardMPGDEndcapRawHitAssociationFrame", "ForwardMPGDEndcapRawHitAssociationFrame",
      "SiBarrelVertexRawHitAssociationFrame",     "SiBarrelRawHitAssociationFrame",
      "SiEndcapTrackerRawHitAssociationFrame",    "B0TrackerRawHitAssociationFrame"};
  std::vector<std::string> m_asso_out_names = {
      "TOFBarrelRawHitAssociations",       "TOFEndcapRawHitAssociations",
      "MPGDBarrelRawHitAssociations",      "OuterMPGDBarrelRawHitAssociations",
      "BackwardMPGDEndcapRawHitAssociations", "ForwardMPGDEndcapRawHitAssociations",
      "SiBarrelVertexRawHitAssociations",  "SiBarrelRawHitAssociations",
      "SiEndcapTrackerRawHitAssociations", "B0TrackerRawHitAssociations"};

  // -- calo hits + per-hit truth (11 detectors) --------------------------------
  // Same pattern as the tracker block above, one level below the clusters the
  // eventbuilder already carries: per-candidate time-gated CalorimeterHits
  // (the GNN's calo input features) and their MCRecoCalorimeterHitAssociations
  // (per-hit signal/background truth for segmentation). Every detector plugin
  // provides the frame-level mirror chain (…RecHitFrame /
  // …RawHitAssociationFrame, Timeslice level) -- see e.g. BEMC.cc.
  // EcalLumiSpec is the one calo family in the GNN input list with no frame
  // mirror; it is omitted here.
  std::vector<std::string> m_calhit_in_names = {
      "EcalBarrelScFiRecHitFrame",   "EcalBarrelImagingRecHitFrame",
      "EcalEndcapNRecHitFrame",      "EcalEndcapPRecHitFrame",
      "B0ECalRecHitFrame",           "HcalBarrelRecHitFrame",
      "HcalEndcapNRecHitFrame",      "HcalEndcapPInsertRecHitFrame",
      "LFHCALRecHitFrame",           "EcalFarForwardZDCRecHitFrame",
      "HcalFarForwardZDCRecHitFrame"};
  std::vector<std::string> m_calhit_out_names = {
      "EcalBarrelScFiRecHits",   "EcalBarrelImagingRecHits",
      "EcalEndcapNRecHits",      "EcalEndcapPRecHits",
      "B0ECalRecHits",           "HcalBarrelRecHits",
      "HcalEndcapNRecHits",      "HcalEndcapPInsertRecHits",
      "LFHCALRecHits",           "EcalFarForwardZDCRecHits",
      "HcalFarForwardZDCRecHits"};
  std::vector<std::string> m_calasso_in_names = {
      "EcalBarrelScFiRawHitAssociationFrame",   "EcalBarrelImagingRawHitAssociationFrame",
      "EcalEndcapNRawHitAssociationFrame",      "EcalEndcapPRawHitAssociationFrame",
      "B0ECalRawHitAssociationFrame",           "HcalBarrelRawHitAssociationFrame",
      "HcalEndcapNRawHitAssociationFrame",      "HcalEndcapPInsertRawHitAssociationFrame",
      "LFHCALRawHitAssociationFrame",           "EcalFarForwardZDCRawHitAssociationFrame",
      "HcalFarForwardZDCRawHitAssociationFrame"};
  std::vector<std::string> m_calasso_out_names = {
      "EcalBarrelScFiRawHitAssociations",   "EcalBarrelImagingRawHitAssociations",
      "EcalEndcapNRawHitAssociations",      "EcalEndcapPRawHitAssociations",
      "B0ECalRawHitAssociations",           "HcalBarrelRawHitAssociations",
      "HcalEndcapNRawHitAssociations",      "HcalEndcapPInsertRawHitAssociations",
      "LFHCALRawHitAssociations",           "EcalFarForwardZDCRawHitAssociations",
      "HcalFarForwardZDCRawHitAssociations"};

  // -- TOF Measurement2D (LGAD cluster hits) — the ACTS measurement inputs -----
  // Central tracking consumes TOF as Measurement2D ("TOFBarrelClusterHits" /
  // "TOFEndcapClusterHits", from LGADHitClustering), NOT as TrackerHits.
  // Project the frame-level mirrors into each child, gated around t0, so the
  // CKF chain on unfolded candidates finds them instead of re-running the
  // digitization chain (which aborts on missing sim inputs — see README).
  std::vector<std::string> m_meas_in_names = {
      "TOFBarrelClusterHitFrame", "TOFEndcapClusterHitFrame"};
  std::vector<std::string> m_meas_out_names = {
      "TOFBarrelClusterHits", "TOFEndcapClusterHits"};

  // -- calo clusters + cluster->particle associations (4 detectors) ------------
  // Read from CalCoincidence_factory output (*CoincClusters) -- already
  // candidate-t0 gated (nsigma x per-cluster getTimeError() + late asym) AND
  // energy-floored per the YR minimum detectable-photon energy
  // (eventbuilder:calcoincidence:*). EcalLumiSpec stays excluded, same as before.
  std::vector<std::string> m_clu_in_names = {
      "B0ECalCoincClusters", "EcalBarrelCoincClusters",
      "EcalEndcapNCoincClusters", "EcalEndcapPCoincClusters",
      "HcalBarrelCoincClusters", "HcalEndcapPInsertCoincClusters",
      "LFHCALCoincClusters", "EcalFarForwardZDCCoincClusters",
      "HcalFarForwardZDCCoincClusters"};
  std::vector<std::string> m_clu_out_names = {
      "B0ECalClusters", "EcalBarrelClusters", "EcalEndcapNClusters", "EcalEndcapPClusters",
      "HcalBarrelClusters", "HcalEndcapPInsertClusters", "LFHCALClusters",
      "EcalFarForwardZDCClusters", "HcalFarForwardZDCClusters"};
  std::vector<std::string> m_cluasso_in_names = {
      "B0ECalClusterAssociationFrame", "EcalBarrelScFiClusterAssociationFrame",
      "EcalEndcapNClusterAssociationFrame", "EcalEndcapPClusterAssociationFrame",
      "HcalBarrelClusterAssociationFrame", "HcalEndcapPInsertClusterAssociationFrame",
      "LFHCALClusterAssociationFrame", "EcalFarForwardZDCClusterAssociationFrame",
      "HcalFarForwardZDCClusterAssociationFrame"};
  std::vector<std::string> m_cluasso_out_names = {
      "B0ECalClusterAssociations", "EcalBarrelClusterAssociations",
      "EcalEndcapNClusterAssociations", "EcalEndcapPClusterAssociations",
      "HcalBarrelClusterAssociations", "HcalEndcapPInsertClusterAssociations",
      "LFHCALClusterAssociations", "EcalFarForwardZDCClusterAssociations",
      "HcalFarForwardZDCClusterAssociations"};

  // -- candidate list from EventBuilder_factory --------------------------------
  // All inputs live in the parent time-frame (Timeslice level); outputs are
  // written into the child PhysicsEvent. Hence .level = Timeslice on every
  // input so JANA fetches them from the parent rather than the (empty) child.
  PodioInput<edm4hep::EventHeader> m_candidates_in{
      this, {.name = "EventCandidatesFiltered", .level = JEventLevel::Timeslice}};

  // Debug pass-through of both candidate-list stages (see TimeAlign/TimeCoinc
  // note below — same level-boundary problem, same fix). EventCandidates is the
  // pre-prefilter list from EventBuilder_factory; EventCandidatesFiltered is
  // what m_candidates_in already reads. Re-exported as PhysicsEvent-level subset
  // copies purely so they show up in the output ROOT file for inspection.
  PodioInput<edm4hep::EventHeader> m_candidates_raw_in{
      this, {.name = "EventCandidates", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_candidates_raw_out{this, "EventCandidates"};
  PodioOutput<edm4hep::EventHeader> m_candidates_filtered_out{this, "EventCandidatesFiltered"};

  // Child-level header naming:
  //   EventHeader      — the child PhysicsEvent's OWN header, standard
  //                      edm4hep semantics only (eventNumber, runNumber,
  //                      timeStamp = t0 in ps, weight = 1). Safe for any
  //                      standard tool: no physics payload hidden in weights.
  //   TimesliceHeader  — passthrough of the parent time-frame's EventHeader,
  //                      kept for provenance.
  //   EventBuilderInfo — the eventbuilder-specific payload (t0, t0sigma,
  //                      significance, truth labels...), EventHeader-typed but
  //                      deliberately NOT named EventHeader so nothing
  //                      standard mistakes its weights vector for MC
  //                      generator weights. Layout below at the fill site.
  PodioInput<edm4hep::EventHeader> m_event_header_in{
      this, {.name = "EventHeader", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_timeslice_header_out{this, "TimesliceHeader"};
  PodioOutput<edm4hep::EventHeader> m_child_header_out{this, "EventHeader"};
  PodioOutput<edm4hep::EventHeader> m_info_out{this, "EventBuilderInfo"};
  // Frame-level GNN prefilter scores (weights[0..6] = 7-class softmax, class
  // 0 = BKG), produced by EventPrefilter_factory when a model runs; copied
  // verbatim onto each child PhysicsEvent so the prefilter cut/class-tag is a
  // reversible offline operation needing no ONNX reload. Optional/empty when
  // the prefilter is off or in veto mode with no model.
  PodioInput<edm4hep::EventHeader> m_scores_in{
      this, {.name = "PrefilterScores", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_scores_out{this, "PrefilterScores"};

  PodioInput<edm4hep::MCParticle> m_mcparticles_in{
      this, {.name = "MCParticles", .level = JEventLevel::Timeslice}};
  PodioOutput<edm4hep::MCParticle> m_mcparticles_out{this, "MCParticles"};

  VariadicPodioInput<edm4eic::TrackerHit> m_trk_in{
      this, {.names = m_trk_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::TrackerHit> m_trk_out{this, m_trk_out_names};

  // Declared BEFORE m_raw_in/m_link_in/m_asso_in below on purpose: JANA
  // populates a component's declared inputs in declaration order, each via a
  // plain "does a databundle with this name already exist" lookup that (for
  // an is_optional input) silently gives up rather than forcing computation
  // if not -- see PodioInput<T>::Populate() in JHasInputs.h. TOF's raw-hit
  // truth is fed from the "Shared" charge-sharing chain (TOFBarrel/EndcapShared
  // RawHit*Frame in BTOF.cc/ECTOF.cc), which is the SAME upstream factory this
  // Measurement2D input below pulls through TOFBarrel/EndcapClusterHitFrame.
  // Populating m_meas_in first forces that factory to actually run, so by the
  // time m_raw_in/m_link_in/m_asso_in look for its outputs a few lines down,
  // they're already there (bug found + fixed 2026-07-24: every TOFBarrel/
  // Endcap RawHit/RawHitAssociation/RawHitLink came out silently empty with
  // the declaration order reversed, while TOFBarrelClusterHits -- populated
  // via this same m_meas_in -- worked fine every time).
  //
  // B0ECAL, BEMC, EEMC are all registered at Timeslice level; FEMC (EcalEndcapP)
  // has no plugin yet so its optional inputs will be nullptr (handled below).
  VariadicPodioInput<edm4eic::Measurement2D> m_meas_in{
      this, {.names = m_meas_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::Measurement2D> m_meas_out{this, m_meas_out_names};

  VariadicPodioInput<edm4eic::RawTrackerHit> m_raw_in{
      this, {.names = m_raw_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::RawTrackerHit> m_raw_out{this, m_raw_out_names};

  VariadicPodioInput<podio::Link<edm4eic::RawTrackerHit, edm4hep::SimTrackerHit>> m_link_in{
      this, {.names = m_link_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<podio::Link<edm4eic::RawTrackerHit, edm4hep::SimTrackerHit>> m_link_out{
      this, m_link_out_names};

  VariadicPodioInput<edm4eic::MCRecoTrackerHitAssociation> m_asso_in{
      this, {.names = m_asso_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::MCRecoTrackerHitAssociation> m_asso_out{this, m_asso_out_names};

  VariadicPodioInput<edm4eic::CalorimeterHit> m_calhit_in{
      this, {.names = m_calhit_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::CalorimeterHit> m_calhit_out{this, m_calhit_out_names};

  VariadicPodioInput<edm4eic::MCRecoCalorimeterHitAssociation> m_calasso_in{
      this, {.names = m_calasso_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::MCRecoCalorimeterHitAssociation> m_calasso_out{
      this, m_calasso_out_names};

  VariadicPodioInput<edm4eic::Cluster> m_clu_in{
      this, {.names = m_clu_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::Cluster> m_clu_out{this, m_clu_out_names};

  VariadicPodioInput<edm4eic::MCRecoClusterParticleAssociation> m_cluasso_in{
      this, {.names = m_cluasso_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::MCRecoClusterParticleAssociation> m_cluasso_out{
      this, m_cluasso_out_names};

  // -- TimeAlign/TimeCoinc debug pass-through -----------------------------------
  // m_trk_in/m_clu_in above only exist at Timeslice level, so the (PhysicsEvent-
  // level) PODIO writer can never see them there: JEventProcessorPODIO enumerates
  // collections via event->GetAllCollectionNames(), which only inspects the
  // current JEvent's own factory set and never walks up to the parent. Re-export
  // them here, under their original names, as full (ungated) subset copies so
  // they land in the output ROOT file for inspection. Since these are Timeslice-
  // wide collections, every PhysicsEvent child of the same parent carries an
  // identical copy — acceptable for debugging, not meant for production output.
  VariadicPodioOutput<edm4eic::TrackerHit> m_trk_align_out{this, m_trk_in_names};
  VariadicPodioOutput<edm4eic::Cluster> m_clu_align_out{this, m_clu_in_names};

  // Ungated slow-detector *TimeAlignRecHits (pre-coincidence-gate). Produced by
  // TrkTimeAlignment_factory but only ever consumed internally by
  // TimeCoincidence_factory; m_trk_in above only reads the post-gate
  // *TimeCoincRecHits version for Si/B0 (indices 6-9). Re-exported here purely
  // to debug the coincidence gate itself (compare against *TimeCoincRecHits).
  std::vector<std::string> m_slow_align_names = {
      "SiBarrelVertexTimeAlignRecHits", "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits", "B0TrackerTimeAlignRecHits"};
  VariadicPodioInput<edm4eic::TrackerHit> m_slow_align_in{
      this, {.names = m_slow_align_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::TrackerHit> m_slow_align_out{this, m_slow_align_names};

  EventUnfolder() {
    SetTypeName(NAME_OF_THIS);
    SetParentLevel(JEventLevel::Timeslice);
    SetChildLevel(JEventLevel::PhysicsEvent);
  }

  // Slow-detector (Si/B0) gate half-width for one hit -- per-hit getTimeError().
  // (not const: Parameter::operator() is non-const)
  double slowHitWindow(double hit_sigma) {
    return double(trk_nsigma_window()) * hit_sigma;
  }

  // Fast (trigger) detectors are indices 0-5; slow (Si/B0) are 6-9.
  static constexpr size_t kNumFastDet = 6;

  // Validate the forward_frames configuration once. Waiting for a future
  // frame's deposit only terminates if that frame can be processed while this
  // one sits in the (sequential) unfold stage. When the pipeline cannot
  // support it (single-threaded), degrade to 0 with a loud warning rather
  // than aborting -- forward recovery is a default-on physics improvement,
  // not a hard requirement.
  int32_t m_fwd_effective = 0;
  void checkConfig() {
    if (m_checked_config)
      return;
    m_checked_config = true;
    const int32_t n_fwd = forward_frames();
    m_fwd_effective = n_fwd;
    if (n_fwd <= 0)
      return;
    auto* app = GetApplication();
    const int nthreads = app->GetNThreads();
    // jana:max_inflight_timeslices defaults to nthreads but is NOT registered
    // unless the user sets it -- GetParameterValue on an unset parameter
    // throws, so probe first and fall back to JANA's effective default.
    int inflight = nthreads;
    if (app->GetJParameterManager()->Exists("jana:max_inflight_timeslices"))
      inflight = app->GetParameterValue<int>("jana:max_inflight_timeslices");
    if (nthreads < 2 || inflight <= n_fwd) {
      m_fwd_effective = 0;
      jerr << "EventUnfolder: forward_frames=" << n_fwd
           << " needs nthreads >= 2 and jana:max_inflight_timeslices > forward_frames"
           << " (got nthreads=" << nthreads << ", max_inflight=" << inflight << ")"
           << " -- running with forward_frames=0; late-tail hits spilling into"
           << " the next frame will NOT be recovered" << jendl;
    }
  }

  Result Unfold(const JEvent& parent, JEvent& child, int child_idx) override {
    const auto& cands = *m_candidates_in();
    const int n_cand  = static_cast<int>(cands.size());
    if (n_cand == 0 || child_idx >= n_cand)
      return Result::KeepChildNextParent;

    auto cand         = cands.at(child_idx);
    const double t0   = cand.getWeights(0); // ns
    const double t0sigma  = cand.getWeights(1); // ns
    const double phys = cand.getWeights(2);

    // -- adjacent-frame Si/B0 hits (cross-frame boundary recovery) -------------
    // Fetched once per parent from the TimesliceBuffer_service. Keyed by
    // absolute frame number, so this is correct regardless of the order in
    // which frames were processed upstream.
    if (!m_ts_buffer)
      m_ts_buffer = GetApplication()->GetService<TimesliceBuffer_service>();
    const uint64_t parent_nr = parent.GetEventNumber();
    if (parent_nr != m_adjacent_parent) {
      m_adjacent_parent = parent_nr;
      m_adjacent.clear();
      checkConfig();
      const int32_t n_back = std::max(0, backward_frames());
      const int32_t n_fwd  = std::max(0, m_fwd_effective);
      for (int32_t d = -n_back; d <= n_fwd; ++d) {
        if (d == 0 || (d < 0 && parent_nr < static_cast<uint64_t>(-d)))
          continue;
        auto frame = m_ts_buffer->fetch(parent_nr + d);
        if (frame)
          m_adjacent.push_back(std::move(*frame));
      }
    }

    // -- tracker RecHits: per-detector resolution-matched window ---------------
    // Fast dets: subset of the parent's current-frame collections, windowed
    // per-hit by that hit's own getTimeError() (see timeResolution_Silicon
    // comment above -- fast dets have no group-level constant to fall back on,
    // nor do they need one). Slow dets: owned clones, merging the current frame
    // with the adjacent frames above (hits from other frames cannot be subset
    // references), windowed per-hit like the fast dets (or by the
    // timeResolution_Silicon constant override when set).
    // Per-detector set of the RawTrackerHits backing THIS candidate's own
    // gated RecHits below -- used to gate the raw hits/links/associations
    // carry-through further down to just this candidate, instead of the
    // whole frame (see copyGatedByRawHit's comment).
    std::vector<std::unordered_set<uint64_t>> kept_raw(m_trk_in().size());

    for (size_t det = 0; det < m_trk_in().size(); ++det) {
      auto& out = m_trk_out().at(det);
      const auto* in = m_trk_in().at(det);

      if (det < kNumFastDet) {
        // Owned clones, not a subset: fast-det inputs are Timeslice-level
        // collections never written to the child tree, so a subset output
        // degrades to an unreadable <name>_objIdx stub offline (this is why
        // TOFBarrel/MPGD RecHits looked "missing" in every output while the
        // slow dets below -- always cloned -- came through; found 2026-07-24).
        if (in == nullptr)
          continue;
        // in is already candidate-t0 gated by TrkCoincidence_factory (see
        // m_trk_in_names above) -- clone through, no re-gating here.
        for (const auto& hit : *in) {
          out->push_back(hit.clone(true));
          kept_raw[det].insert(rawHitKey(hit.getRawHit()));
        }
        continue;
      }

      std::vector<edm4eic::MutableTrackerHit> merged;
      const double late_nsigma_window_asym = double(trk_late_nsigma_window_asym());
      if (in != nullptr)
        for (const auto& hit : *in) {
          const double sigma = double(hit.getTimeError());
          const double w = slowHitWindow(sigma);
          const double dt = hit.getTime() - t0;
          if (dt >= -w && dt <= w + late_nsigma_window_asym * sigma) {
            merged.push_back(hit.clone());
            kept_raw[det].insert(rawHitKey(hit.getRawHit()));
          }
        }
      // BUG-002 fix: adjacent-frame hits are cloned WITHOUT their relations
      // (cloneRelations=false), unlike the current-frame clone above. A
      // default hit.clone() keeps the relation pointing at the ORIGINAL
      // RawTrackerHit, which here lives in one of the m_adjacent frames --
      // and m_adjacent is cleared/replaced (line ~415, above) the moment
      // the unfolder moves to the next parent timeslice. If this child's
      // TrackerHit is still unwritten when that happens (podio's writer can
      // lag the unfolder), prepareForWrite() later dereferences a
      // now-destroyed RawTrackerHit and segfaults
      // (edm4eic::RawTrackerHit::getObjectID() on a dangling object) --
      // observed specifically on pileup data, where denser cross-frame
      // content makes this race far more likely to land inside the window.
      // The current-frame clone above is safe as-is: its source collection
      // (m_trk_in(), owned by the parent event) outlives this child by
      // construction, matching the fast-det subset path just above, which
      // shows zero faults across every test run.
      const size_t si = det - kNumFastDet;
      for (const auto& frame : m_adjacent)
        if (si < frame.size())
          for (const auto& hit : frame[si]) {
            const double sigma = double(hit.getTimeError());
            const double w = slowHitWindow(sigma);
            const double dt = hit.getTime() - t0;
            if (dt >= -w && dt <= w + late_nsigma_window_asym * sigma)
              merged.push_back(hit.clone(false));
          }
      std::sort(merged.begin(), merged.end(),
                [](const auto& a, const auto& b) { return a.getTime() < b.getTime(); });
      for (auto& hit : merged)
        out->push_back(hit);
    }

    // -- TOF Measurement2D (LGAD cluster hits): gated like the fast hits -------
    // Per-measurement time error = sqrt(covariance.zz) (loc0/loc1/time
    // ordering); fall back to the 25 ps LGAD resolution when unset.
    for (size_t d = 0; d < m_meas_in().size(); ++d) {
      auto& out = m_meas_out().at(d);
      const auto* in = m_meas_in().at(d);
      if (in == nullptr)
        continue;
      for (const auto& m : *in) {
        const double terr  = std::sqrt(std::max(0.0f, m.getCovariance().zz));
        const double sigma = terr > 0.0 ? terr : 0.025;
        const double win   = double(trk_nsigma_window()) * sigma;
        const double dt    = m.getTime() - t0;
        // clone, not subset -- same objIdx-stub trap as the fast dets above
        if (dt >= -win && dt <= win + double(trk_late_nsigma_window_asym()) * sigma)
          out->push_back(m.clone(true));
      }
    }

    // -- raw hits / links / hit associations: gated to THIS candidate ----------
    // Was an unconditional copyAll (every raw hit/link/association in the
    // WHOLE frame, handed identically to every candidate unfolded from it)
    // until 2026-07-24 -- see copyGatedByRawHit's comment for why that's
    // wrong for any per-candidate consumer (e.g. GNN segmentation truth).
    copyGatedByRawHit(m_raw_in(), m_raw_out(), kept_raw,
                       [](const auto& h) -> const edm4eic::RawTrackerHit& { return h; });
    copyGatedByRawHit(m_link_in(), m_link_out(), kept_raw,
                       [](const auto& l) { return l.getFrom(); });
    // Associations get an explicit loop instead of copyGatedByRawHit: the
    // clone's weight is overwritten with the sim particle's generatorStatus.
    // Rationale: the association's simHit relation points into the parent
    // frame's SimTrackerHit collection, which is never written to the child
    // tree -- so the assoc->simHit->MCParticle chain that turns a hit into a
    // signal/background label is UNRESOLVABLE offline (and re-exporting the
    // frame collections per candidate would both bloat the file and break
    // ObjectID index alignment against the child's gated MCParticles clone).
    // This is the one place the whole chain is still in memory, so bake the
    // label in here. The SBM production encodes provenance in
    // generatorStatus: <2000 = signal event (unshifted), 2000-6999 = machine
    // background (per-source offset), >=7000 = cross-class pileup (per-class
    // offset) -- see sro/scripts/signal.sh. The original weight is a
    // constant 1.0 from SiliconTrackerDigi, so nothing meaningful is lost;
    // consumers who need an unmodified truth channel still have the
    // *RawHitLinks written alongside. Offline join: match an association to
    // a RecHit by rawHit ObjectID (collectionID, index) equality -- both
    // clones keep that relation, and key equality needs no resolution.
    for (size_t i = 0; i < m_asso_out().size(); ++i) {
      if (i >= m_asso_in().size() || m_asso_in().at(i) == nullptr || i >= kept_raw.size())
        continue;
      const auto& keep = kept_raw[i];
      for (const auto& a : *m_asso_in().at(i)) {
        if (keep.count(rawHitKey(a.getRawHit())) == 0)
          continue;
        auto clone = a.clone(true);
        const auto sim = a.getSimHit();
        if (sim.isAvailable() && sim.getParticle().isAvailable())
          clone.setWeight(static_cast<float>(sim.getParticle().getGeneratorStatus()));
        m_asso_out().at(i)->push_back(clone);
      }
    }

    // -- calo hits + per-hit truth: gated to THIS candidate --------------------
    // Mirrors the tracker block above at the hit level. Hits are time-gated
    // inline with cal_nsigma_window/cal_late_nsigma_window_asym (no
    // coincidence factory exists for raw CalorimeterHit); the clusters
    // further down are a DIFFERENT channel, now pre-gated upstream by
    // CalCoincidence_factory instead -- see cal_nsigma_window's comment
    // above. Associations are then keyed on the
    // gated hits' rawHit ObjectIDs and cloned with the label-in-weight
    // convention of the tracker associations. One calo-specific wrinkle: a
    // SimCalorimeterHit can mix energy deposits from SEVERAL particles in one
    // cell (contributions), so "the" label is taken from the HIGHEST-ENERGY
    // contribution's particle -- the dominant depositor decides the cell.
    std::vector<std::unordered_set<uint64_t>> kept_calraw(m_calhit_in().size());
    for (size_t d = 0; d < m_calhit_in().size(); ++d) {
      auto& out = m_calhit_out().at(d);
      const auto* in = m_calhit_in().at(d);
      if (in == nullptr)
        continue;
      for (const auto& hit : *in) {
        const double terr  = double(hit.getTimeError());
        // Guard sigma<=0 (older digi chains hardcoded hit timeError to 0 --
        // see the CalorimeterHitReco note above): fall back to 1 ns, a
        // calo-scale width, instead of silently gating the whole family away.
        const double sigma = terr > 0.0 ? terr : 1.0;
        const double win   = double(cal_nsigma_window()) * sigma;
        const double dt    = hit.getTime() - t0;
        if (dt >= -win && dt <= win + double(cal_late_nsigma_window_asym()) * sigma) {
          out->push_back(hit.clone(true));
          if (hit.getRawHit().isAvailable())
            kept_calraw[d].insert(rawHitKey(hit.getRawHit()));
        }
      }
    }
    for (size_t i = 0; i < m_calasso_out().size(); ++i) {
      if (i >= m_calasso_in().size() || m_calasso_in().at(i) == nullptr || i >= kept_calraw.size())
        continue;
      const auto& keep = kept_calraw[i];
      for (const auto& a : *m_calasso_in().at(i)) {
        if (keep.count(rawHitKey(a.getRawHit())) == 0)
          continue;
        auto clone = a.clone(true);
        const auto sim = a.getSimHit();
        if (sim.isAvailable()) {
          float best_e = -1.0F;
          int32_t status = 0;
          bool found = false;
          for (const auto& c : sim.getContributions()) {
            if (c.getEnergy() > best_e && c.getParticle().isAvailable()) {
              best_e = c.getEnergy();
              status = c.getParticle().getGeneratorStatus();
              found  = true;
            }
          }
          if (found)
            clone.setWeight(static_cast<float>(status));
        }
        m_calasso_out().at(i)->push_back(clone);
      }
    }

    if (debug_gating()) {
      for (size_t d = 0; d < m_trk_in_names.size(); ++d) {
        const auto* tin = d < m_trk_in().size() ? m_trk_in().at(d) : nullptr;
        const auto* rin = d < m_raw_in().size() ? m_raw_in().at(d) : nullptr;
        const auto* lin = d < m_link_in().size() ? m_link_in().at(d) : nullptr;
        const auto* ain = d < m_asso_in().size() ? m_asso_in().at(d) : nullptr;
        jout << "[eb:gating] cand=" << child_idx << " det=" << m_trk_in_names[d]
             << " trk_in=" << (tin ? std::to_string(tin->size()) : "null")
             << " trk_out=" << m_trk_out().at(d)->size()
             << " kept_raw=" << kept_raw[d].size()
             << " raw_in=" << (rin ? std::to_string(rin->size()) : "null")
             << " raw_out=" << m_raw_out().at(d)->size()
             << " link_in=" << (lin ? std::to_string(lin->size()) : "null")
             << " link_out=" << m_link_out().at(d)->size()
             << " asso_in=" << (ain ? std::to_string(ain->size()) : "null")
             << " asso_out=" << m_asso_out().at(d)->size() << jendl;
      }
      for (size_t d = 0; d < m_calhit_in_names.size(); ++d) {
        const auto* hin = d < m_calhit_in().size() ? m_calhit_in().at(d) : nullptr;
        const auto* ain = d < m_calasso_in().size() ? m_calasso_in().at(d) : nullptr;
        jout << "[eb:gating] cand=" << child_idx << " det=" << m_calhit_in_names[d]
             << " hit_in=" << (hin ? std::to_string(hin->size()) : "null")
             << " hit_out=" << m_calhit_out().at(d)->size()
             << " kept_raw=" << kept_calraw[d].size()
             << " asso_in=" << (ain ? std::to_string(ain->size()) : "null")
             << " asso_out=" << m_calasso_out().at(d)->size() << jendl;
      }
    }

    // -- TimeAlign/TimeCoinc debug pass-through (full, ungated) ----------------
    copyAll(m_trk_in(), m_trk_align_out());
    copyAll(m_clu_in(), m_clu_align_out());
    copyAll(m_slow_align_in(), m_slow_align_out());

    // -- candidate-list debug pass-through -------------------------------------
    m_candidates_filtered_out()->setSubsetCollection(true);
    for (const auto& c : *m_candidates_in())
      m_candidates_filtered_out()->push_back(c);

    m_candidates_raw_out()->setSubsetCollection(true);
    if (m_candidates_raw_in() != nullptr)
      for (const auto& c : *m_candidates_raw_in())
        m_candidates_raw_out()->push_back(c);

    // -- MC particles: truth bookkeeping, gated around t0 (see mc_time_window) -
    // Moved ahead of the calo block below: it needs this loop's old-index ->
    // local-clone map to remap MCRecoClusterParticleAssociation.sim.
    const float mc_win = mc_time_window();
    std::unordered_map<uint32_t, edm4hep::MutableMCParticle> mcp_map;
    for (const auto& mcp : *m_mcparticles_in()) {
      if (mc_win >= 0.0f && std::abs(mcp.getTime() - t0) > mc_win)
        continue;
      auto newMcp = mcp.clone(true);
      m_mcparticles_out()->push_back(newMcp);
      mcp_map.emplace(static_cast<uint32_t>(mcp.getObjectID().index), newMcp);
    }

    // -- calo clusters + matched cluster->particle associations ----------------
    // Owned clones, NOT subset references: unlike tracks (recomputed fresh
    // per candidate by the ACTS chain, so their child-level output is already
    // self-contained), calo clusters are only computed once at the
    // Timeslice/frame level. A subset output (the previous behavior) stores
    // just an ObjectID back-reference into that frame-level collection, so
    // .energy/.position/.time/etc. never materialize in the child
    // PhysicsEvent tree -- downstream per-candidate calo efficiency/purity
    // reads the child tree alone and sees nothing. Clone instead, exactly
    // like MCParticles above. cin is already candidate-t0 gated (and energy-
    // floored) by CalCoincidence_factory -- see m_clu_in_names above -- so no
    // re-gating here; previously EVERY cluster in the frame was pushed into
    // EVERY candidate unconditionally, falsely attributing frame-wide clusters
    // to every candidate regardless of timing, which the (now upstream) gate fixes.
    for (size_t cd = 0; cd < m_clu_in().size(); ++cd) {
      auto& cout = m_clu_out().at(cd);
      auto& aout = m_cluasso_out().at(cd);
      const auto* cin = m_clu_in().at(cd);
      if (cin == nullptr) continue; // optional: e.g. FEMC not yet registered
      std::unordered_map<uint32_t, edm4eic::MutableCluster> included;
      for (const auto& clu : *cin) {
        auto newClu = clu.clone(false);
        cout->push_back(newClu);
        included.emplace(static_cast<uint32_t>(clu.getObjectID().index), newClu);
      }
      const auto* ain = (cd < m_cluasso_in().size()) ? m_cluasso_in().at(cd) : nullptr;
      if (ain == nullptr) continue;
      for (const auto& as : *ain) {
        const auto cit = included.find(static_cast<uint32_t>(as.getRec().getObjectID().index));
        if (cit == included.end()) continue; // cluster gated out above
        const auto mit = mcp_map.find(static_cast<uint32_t>(as.getSim().getObjectID().index));
        if (mit == mcp_map.end()) continue; // truth particle gated out (mc_time_window)
        auto newAs = as.clone(false);
        newAs.setRec(cit->second);
        newAs.setSim(mit->second);
        aout->push_back(newAs);
      }
    }

    // -- child headers (EventHeader / TimesliceHeader / EventBuilderInfo) ------
    // Deterministic numbering: with nthreads > 1 parents can reach the unfolder
    // out of order, so a running counter is not reproducible run-to-run. Encode
    // (frame, candidate) instead. Dense background frames have been observed at
    // a few hundred candidates, so guard the assumed <1000 margin loudly: a
    // 1000th candidate would silently collide with the next frame's numbering.
    if (child_idx >= 1000 && !m_warned_candidate_overflow) {
      m_warned_candidate_overflow = true;
      jerr << "EventUnfolder: frame " << parent_nr << " has >= 1000 candidates — "
           << "child eventNumber encoding (frame*1000+idx) is no longer unique; "
           << "tighten the trigger (eventbuilder:topology:threshold / "
              "eventbuilder:trigger:p_fake)" << jendl;
    }
    const uint64_t child_event_nr = parent_nr * 1000 + static_cast<uint64_t>(child_idx);
    child.SetEventNumber(child_event_nr);
    child.SetRunNumber(parent.GetRunNumber());

    // EventHeader: the child's own header, STANDARD edm4hep semantics only.
    edm4hep::MutableEventHeader hdr;
    hdr.setRunNumber(parent.GetRunNumber());
    hdr.setEventNumber(child_event_nr);
    hdr.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    hdr.setWeight(1.0); // real event weight — physics payload lives in EventBuilderInfo
    m_child_header_out()->push_back(hdr);

    // EventBuilderInfo: the candidate's full weights vector (see the layout
    // in EventBuilder_factory::emitCandidate — [0..17] incl. the measured
    // frame background and the applied configuration, making every output
    // file self-describing) copied verbatim, plus ONE appended entry:
    //   [last] = candidates in this frame (pileup context)
    edm4hep::MutableEventHeader info;
    info.setRunNumber(parent.GetRunNumber());
    info.setEventNumber(child_event_nr);
    info.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    for (const auto w : cand.getWeights())
      info.addToWeights(w);
    info.addToWeights(static_cast<double>(n_cand));
    m_info_out()->push_back(info);

    // Copy the frame-level GNN prefilter scores onto this child (owned clone,
    // cross-level). One frame → one score vector, replicated per child so it
    // is trivially usable in flat per-event offline cuts. Empty when the
    // prefilter did not run (no model) — the collection is still created.
    if (m_scores_in() != nullptr && m_scores_in()->size() > 0) {
      edm4hep::MutableEventHeader sc;
      sc.setRunNumber(parent.GetRunNumber());
      sc.setEventNumber(child_event_nr);
      for (const auto w : m_scores_in()->at(0).getWeights())
        sc.addToWeights(w);
      m_scores_out()->push_back(sc);
    }

    // TimesliceHeader: the parent time-frame's header, for provenance.
    if (m_event_header_in() != nullptr && m_event_header_in()->size() > 0) {
      m_timeslice_header_out()->setSubsetCollection(true);
      m_timeslice_header_out()->push_back(m_event_header_in()->at(0));
    }

    if (child_idx == n_cand - 1)
      return Result::NextChildNextParent;
    return Result::NextChildKeepParent;
  }

  // Copy every element of each input collection into the matching output subset
  // collection (guarding against missing optional inputs).
  template <typename InVec, typename OutVec>
  void copyAll(const InVec& in, OutVec& out) {
    for (size_t i = 0; i < out.size(); ++i) {
      out.at(i)->setSubsetCollection(true);
      if (i >= in.size() || in.at(i) == nullptr)
        continue;
      for (const auto& elem : *in.at(i))
        out.at(i)->push_back(elem);
    }
  }

  // Identity key for a raw hit (collectionID + index), for the kept-hit
  // lookup below. Raw hits carry no calibrated time, so they can't be gated
  // against t0 directly the way RecHits are -- this keys off which raw hit
  // each surviving (already t0-gated) RecHit points back to instead.
  // Templated: tracker uses edm4eic::RawTrackerHit, calo uses
  // edm4hep::RawCalorimeterHit; any podio object with getObjectID works.
  template <typename ObjT>
  static uint64_t rawHitKey(const ObjT& h) {
    const auto id = h.getObjectID();
    return (static_cast<uint64_t>(static_cast<uint32_t>(id.collectionID)) << 32) |
           static_cast<uint32_t>(id.index);
  }

  // Like copyAll, but only for elements whose raw hit (found via rawHitOf)
  // is one of THIS candidate's own gated RecHits (kept.at(i)) -- undoes the
  // previous unconditional carry-through, which silently gave every
  // candidate unfolded from one frame the frame's FULL raw-hit/truth
  // population regardless of its own t0. Only meaningful for per-candidate
  // consumers (e.g. GNN segmentation truth); the raw hits/links/associations
  // themselves carry no r/c correction, so this is the only gating available
  // at this level -- see the class comment above m_raw_in_names for why
  // cross-frame buffering doesn't apply here either.
  //
  // OWNED CLONES, not a subset: same trap as the calo clusters below -- a
  // subset output of a Timeslice-level collection materializes in the child
  // tree as a bare <name>_objIdx reference stub (the referenced frame-level
  // collection is never written), so nothing is readable offline. clone(true)
  // keeps the relations (rawHit/simHit/from/to), which point at parent-frame
  // collections that outlive this child by construction (see the BUG-002
  // note above for the lifetime argument).
  template <typename InVec, typename OutVec, typename RawHitOf>
  void copyGatedByRawHit(const InVec& in, OutVec& out,
                         const std::vector<std::unordered_set<uint64_t>>& kept,
                         RawHitOf&& rawHitOf) {
    for (size_t i = 0; i < out.size(); ++i) {
      if (i >= in.size() || in.at(i) == nullptr || i >= kept.size())
        continue;
      const auto& keep = kept[i];
      for (const auto& elem : *in.at(i))
        if (keep.count(rawHitKey(rawHitOf(elem))) > 0)
          out.at(i)->push_back(elem.clone(true));
    }
  }
};


extern "C" {
void InitPlugin(JApplication* app) {

  std::vector<std::string> m_simtrackerhit_collection_names_aligned = {
      "TOFBarrelTimeAlignRecHits",
      "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",
      "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits",
      "ForwardMPGDEndcapTimeAlignRecHits",
      "SiBarrelVertexTimeAlignRecHits",
      "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits",
      "B0TrackerTimeAlignRecHits"
  };
  // "TaggerTrackerTimeAlignRecHits",
  // "DIRCBarTimeAlignRecHits",
  // "DRICHTimeAlignRecHits",
  // "ForwardOffMTrackerTimeAlignRecHits",
  // "ForwardRomanPotTimeAlignRecHits",
  // "LumiSpecTrackerTimeAlignRecHits",
  // "RICHEndcapNTimeAlignRecHits"

  std::vector<std::string> m_simtrackerhit_collection_names = {
      "TOFBarrelRecHitFrame",
      "TOFEndcapRecHitFrame",
      "MPGDBarrelRecHitFrame",
      "OuterMPGDBarrelRecHitFrame",
      "BackwardMPGDEndcapRecHitFrame",
      "ForwardMPGDEndcapRecHitFrame",
      "SiBarrelVertexRecHitFrame",
      "SiBarrelTrackerRecHitFrame",
      "SiEndcapTrackerRecHitFrame",
      "B0TrackerRecHitFrame"
  };
  // "TaggerTrackerRecHitFrame",
  // "DIRCBarRecHitFrame",
  // "DRICHRecHitFrame",
  // "ForwardOffMTrackerRecHitFrame",
  // "ForwardRomanPotRecHitFrame",
  // "LumiSpecRecHitFrame",
  // "RICHEndcapNRecHitFrame"

  std::vector<std::string> m_simcalocluster_collection_names_aligned = {
      "B0ECalTimeAlignClusters",
      "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters",
      "EcalEndcapPTimeAlignClusters",
      "HcalBarrelTimeAlignClusters",
      "HcalEndcapPInsertTimeAlignClusters",
      "LFHCALTimeAlignClusters",
      "EcalFarForwardZDCTimeAlignClusters",
      "HcalFarForwardZDCTimeAlignClusters"
  };
  // "EcalFarForwardZDCTimeAlignClusters",
  // "EcalLumiSpecTimeAlignClusters",
  // "HcalBarrelTimeAlignClusters",
  // "HcalEndcapNTimeAlignClusters",
  // "HcalEndcapPInsertTimeAlignClusters",
  // "HcalFarForwardZDCTimeAlignClusters",
  // "LFHCALTimeAlignClusters"

  std::vector<std::string> m_simcalocluster_collection_names = {
      "B0ECalClusterFrame",
      // ScFi clusters, NOT the merged "EcalBarrelClusterFrame": the merger
      // (EnergyPositionClusterMerger, ScFi x Imaging pair matching) is
      // starved at Timeslice level — measured 5 barrel clusters in 19 frames
      // vs ~50 per endcap, and a 0.65 GeV electron shower produced none —
      // so barrel-only events silently lost their E_calo tag. The ScFi side
      // carries the EM energy and clusters reliably; the TimeAlign OUTPUT
      // name (EcalBarrelTimeAlignClusters) is unchanged, so every consumer
      // (E_calo weights[24], tag bitmask, unfolder) picks this up as-is.
      "EcalBarrelScFiClusterFrame",
      "EcalEndcapNClusterFrame",
      "EcalEndcapPClusterFrame",
      "HcalBarrelClusterFrame",
      "HcalEndcapPInsertClusterFrame",
      "LFHCALClusterFrame",
      "EcalFarForwardZDCClusterFrame",
      "HcalFarForwardZDCClusterFrame"
  };
  // "EcalLumiSpecClusterFrame",
  // "HcalEndcapNClusterFrame",
  // "EcalBarrelImagingClusterFrame",
  // "EcalBarrelScFiClusterFrame",
  // "EcalEndcapNImagingClusterFrame",
  // "EcalEndcapPImagingClusterFrame",
  // "EcalFarForwardZDCImagingClusterFrame",
  // "EcalLumiSpecImagingClusterFrame"

  InitJANAPlugin(app);

  // Process-wide services: cross-frame Si/B0 hit buffer (thread-safe, keyed by
  // frame number) and the shared ONNX runtime used by the prefilter.
  app->ProvideService(std::make_shared<TimesliceBuffer_service>());
  app->ProvideService(std::make_shared<ONNXRuntime_service>(app));

  app->Add(new JOmniFactoryGeneratorT<TrkTimeAlignment_factory>(
      JOmniFactoryGeneratorT<TrkTimeAlignment_factory>::TypedWiring{
          .m_tag                 = "align",
          .m_default_input_tags  = m_simtrackerhit_collection_names,
          .m_default_output_tags = m_simtrackerhit_collection_names_aligned,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<CalTimeAlignment_factory>(
      JOmniFactoryGeneratorT<CalTimeAlignment_factory>::TypedWiring{
          .m_tag                 = "CalTimeAlignment",
          .m_default_input_tags  = m_simcalocluster_collection_names,
          .m_default_output_tags = m_simcalocluster_collection_names_aligned,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // EventBuilder_factory: finds candidate physics events in a time-frame and
  // assigns each a (t0, t0sigma). Runs at the Timeslice (time-frame) level.
  // Only fast detectors (TOF, MPGD, indices 0-5) participate in event finding
  // (kNumTriggerDet=6 in EventBuilder_factory.h). Si/B0 are wired as optional
  // inputs for completeness but are never read in the finding loop.
  // Cross-frame slow-hit recovery for RecHits is handled by the EventUnfolder
  // downstream. Raw hits use current frame only (no cross-frame buffer).
  // NOTE: this flat tag list is distributed POSITIONALLY over the factory's
  // declared inputs (10-slot tracker variadic, then MCParticles, then the
  // 9-slot calo-cluster variadic for the E_calo score, weights[24]) — adding
  // an input to EventBuilder_factory.h without extending this list shifts
  // every later binding and silently nulls optional inputs (observed as
  // "all candidates fake + E_calo always 0").
  std::vector<std::string> m_eventbuilder_input_tags = {
      "TOFBarrelTimeAlignRecHits",
      "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",
      "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits",
      "ForwardMPGDEndcapTimeAlignRecHits",
      "SiBarrelVertexTimeAlignRecHits",
      "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits",
      "B0TrackerTimeAlignRecHits",
      "MCParticles",
      "B0ECalTimeAlignClusters",
      "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters",
      "EcalEndcapPTimeAlignClusters",
      "HcalBarrelTimeAlignClusters",
      "HcalEndcapPInsertTimeAlignClusters",
      "LFHCALTimeAlignClusters",
      "EcalFarForwardZDCTimeAlignClusters",
      "HcalFarForwardZDCTimeAlignClusters",
      // placeholder, not yet wired (optional -> nullptr); present so the two
      // variadic inputs stay the SAME SIZE — JOmniFactory splits this flat
      // list evenly across variadics and throws otherwise.
      "EcalLumiSpecTimeAlignClusters"
  };
  app->Add(new JOmniFactoryGeneratorT<EventBuilder_factory>(
      JOmniFactoryGeneratorT<EventBuilder_factory>::TypedWiring{
          .m_tag                 = "trigger",
          .m_default_input_tags  = m_eventbuilder_input_tags,
          .m_default_output_tags = {"EventCandidates"},
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // TimeCoincidence_factory: candidate-t0 gate for slow (Si/B0) RecHits.
  // Gates the current frame's *TimeAlignRecHits on each candidate's t0 window
  // (nsigma × sigma_hit, sigma_hit = each hit's own getTimeError()) and
  // produces a subset collection for EventUnfolder. This replaced a fixed
  // 2000 ns window, which was 200-600x wider than the real per-hit
  // resolution (see TimeCoincidence_factory.h).
  // Cross-frame boundary recovery happens downstream: EventUnfolder merges the
  // adjacent frames' slow hits from the TimesliceBuffer_service.
  // Input: EventCandidates + 4 slow *TimeAlignRecHits (current frame, Si/B0).
  // Output: 4 *TimeCoincRecHits (subset) → EventUnfolder slow RecHit inputs.
  const std::vector<std::string> coinc_inputs = {
      "EventCandidates",
      "SiBarrelVertexTimeAlignRecHits", "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits", "B0TrackerTimeAlignRecHits"
  };
  const std::vector<std::string> coinc_outputs = {
      "SiBarrelVertexTimeCoincRecHits", "SiBarrelTrackerTimeCoincRecHits",
      "SiEndcapTrackerTimeCoincRecHits", "B0TrackerTimeCoincRecHits"
  };
  app->Add(new JOmniFactoryGeneratorT<TimeCoincidence_factory>(
      JOmniFactoryGeneratorT<TimeCoincidence_factory>::TypedWiring{
          .m_tag                 = "coincidence",
          .m_default_input_tags  = coinc_inputs,
          .m_default_output_tags = coinc_outputs,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // TrkCoincidence_factory: the SAME candidate-t0 gate as TimeCoincidence_
  // factory above, applied to the FAST (TOF+MPGD) hits instead of slow
  // (Si/B0) — discards the frame-wide fast-hit population down to just
  // what falls within some candidate's window. Exists so downstream
  // consumers (EventPrefilter_factory) see a minimized, coincidence-cleaned
  // hit set rather than the whole frame.
  const std::vector<std::string> trkcoinc_inputs = {
      "EventCandidates",
      "TOFBarrelTimeAlignRecHits",          "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",         "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits", "ForwardMPGDEndcapTimeAlignRecHits"
  };
  const std::vector<std::string> trkcoinc_outputs = {
      "TOFBarrelTrkCoincRecHits",          "TOFEndcapTrkCoincRecHits",
      "MPGDBarrelTrkCoincRecHits",         "OuterMPGDBarrelTrkCoincRecHits",
      "BackwardMPGDEndcapTrkCoincRecHits", "ForwardMPGDEndcapTrkCoincRecHits"
  };
  app->Add(new JOmniFactoryGeneratorT<TrkCoincidence_factory>(
      JOmniFactoryGeneratorT<TrkCoincidence_factory>::TypedWiring{
          .m_tag                 = "trkcoincidence",
          .m_default_input_tags  = trkcoinc_inputs,
          .m_default_output_tags = trkcoinc_outputs,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // CalCoincidence_factory: same candidate-t0 gate, applied to calo clusters,
  // plus an (off-by-default) energy floor. Promotes what E_calo currently
  // computes inline as a scalar sum into a real, reusable gated collection.
  // 10 names (matching EventBuilder_factory's own m_calo_collection_names):
  // EcalLumiSpec is anticipated-but-not-yet-wired, resolves to nullptr.
  const std::vector<std::string> calcoinc_inputs = {
      "EventCandidates",
      "B0ECalTimeAlignClusters",            "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters",       "EcalEndcapPTimeAlignClusters",
      "HcalBarrelTimeAlignClusters",        "HcalEndcapPInsertTimeAlignClusters",
      "LFHCALTimeAlignClusters",            "EcalFarForwardZDCTimeAlignClusters",
      "HcalFarForwardZDCTimeAlignClusters", "EcalLumiSpecTimeAlignClusters"
  };
  const std::vector<std::string> calcoinc_outputs = {
      "B0ECalCoincClusters",           "EcalBarrelCoincClusters",
      "EcalEndcapNCoincClusters",      "EcalEndcapPCoincClusters",
      "HcalBarrelCoincClusters",       "HcalEndcapPInsertCoincClusters",
      "LFHCALCoincClusters",           "EcalFarForwardZDCCoincClusters",
      "HcalFarForwardZDCCoincClusters", "EcalLumiSpecCoincClusters"
  };
  app->Add(new JOmniFactoryGeneratorT<CalCoincidence_factory>(
      JOmniFactoryGeneratorT<CalCoincidence_factory>::TypedWiring{
          .m_tag                 = "calcoincidence",
          .m_default_input_tags  = calcoinc_inputs,
          .m_default_output_tags = calcoinc_outputs,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // EventPrefilter_factory: optional GNN background rejection between EventBuilder
  // and EventUnfolder. When eventbuilder:onnx_model is empty (the default), all
  // EventCandidates are forwarded unchanged and no ONNX session is created.
  // When a model path is set, only candidates with
  // probs[BKG=0] < eventbuilder:score_threshold survive.
  //
  // Hit/cluster features per candidate frame: [x/4000, y/4000, z/5000, t/50,
  // eDep-or-energy/10, det_id/6]. Kept in sync with EventPrefilter_factory.h's
  // own hardcoded m_trk_hit_names/m_calo_cluster_names (that factory's
  // VariadicPodioInputs bind their names directly at construction, matching
  // EventBuilder_factory's m_trk_in/m_clu_in pattern) — this wiring-level
  // list is redundant with those but kept matching for safety, same as
  // EventBuilder_factory's own dual-declaration below.
  // ALL hit/cluster inputs are the COINCIDENCE-GATED collections (Trk/Time/
  // CalCoincidence output), not the raw frame-wide *TimeAlign* mirrors —
  // "as little as possible, after cleaning coincidence," not the whole frame.
  const std::vector<std::string> prefilter_inputs = {
      "EventCandidates",
      "TOFBarrelTrkCoincRecHits",          "TOFEndcapTrkCoincRecHits",
      "MPGDBarrelTrkCoincRecHits",         "OuterMPGDBarrelTrkCoincRecHits",
      "BackwardMPGDEndcapTrkCoincRecHits", "ForwardMPGDEndcapTrkCoincRecHits",
      "SiBarrelVertexTimeCoincRecHits",    "SiBarrelTrackerTimeCoincRecHits",
      "SiEndcapTrackerTimeCoincRecHits",   "B0TrackerTimeCoincRecHits",
      "B0ECalCoincClusters",               "EcalBarrelCoincClusters",
      "EcalEndcapNCoincClusters",          "EcalEndcapPCoincClusters",
      "HcalBarrelCoincClusters",           "HcalEndcapPInsertCoincClusters",
      "LFHCALCoincClusters",               "EcalFarForwardZDCCoincClusters",
      "HcalFarForwardZDCCoincClusters",    "EcalLumiSpecCoincClusters"
  };
  app->Add(new JOmniFactoryGeneratorT<EventPrefilter_factory>(
      JOmniFactoryGeneratorT<EventPrefilter_factory>::TypedWiring{
          .m_tag                 = "prefilter",
          .m_default_input_tags  = prefilter_inputs,
          .m_default_output_tags = {"EventCandidatesFiltered", "PrefilterScores"},
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // Unfolder: materialises each surviving candidate as a PhysicsEvent.
  app->Add(new EventUnfolder());

  // Detector plugins (BTOF, MPGD, calorimeters, etc.) are NOT loaded here.
  // eventbuilder depends on their digitization output (per-hit getTimeError(),
  // see EventBuilder_factory.h's derivation comment), including a
  // Timeslice-level, "Frame"-suffixed mirror of each detector's chain (see
  // e.g. src/detectors/BTOF/BTOF.cc) that TimeAlignment_factory reads. That
  // mirror now lives in each detector's own stock InitPlugin, alongside its
  // normal PhysicsEvent-level registration, so eventbuilder does not
  // initialize any of it itself: the standard eicrecon executable already
  // loads every default plugin unconditionally (eicrecon_cli.cc's
  // AddAvailablePluginsToOptionParams), and calling a detector's InitPlugin
  // a second time here would only double-initialize it.
}
} // extern "C"
