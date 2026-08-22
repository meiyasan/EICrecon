// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
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
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/CalorimeterHitCollection.h>
#include <edm4eic/MCRecoCalorimeterHitLinkCollection.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/eventbuilder/TimeAlignment_factory.h"
#include "factories/eventbuilder/TimeCoincidence_factory.h"
#include "factories/eventbuilder/EventBuilder_factory.h"
#include "factories/eventbuilder/EventPrefilter_factory.h"
#include "services/eventbuilder/TimesliceBuffer_service.h"
#include "services/onnx/ONNXRuntime_service.h"


// ---------------------------------------------------------------------------
// EventUnfolder: creates one PhysicsEvent per event candidate.
//
// EventBuilder_factory finds candidates and writes them to "EventCandidates"
// (an edm4hep::EventHeaderCollection with t0 and t0sigma). This class reads
// that list and creates one child PhysicsEvent for each candidate. It does
// no event-finding or coincidence logic itself: only a JEventUnfolder can
// create child events in JANA, so this class does only that step.
// ---------------------------------------------------------------------------
struct EventUnfolder : public JEventUnfolder {

  // Acceptance window for tracker and TOF hits: N sigma around t0, where
  // sigma is the hit's own edm4eic::TrackerHit::getTimeError(). Each
  // detector's own digitization sets this value, so the window always
  // matches that hit's true timing resolution.
  //
  // The window is symmetric. Measured residuals (`make rescheck`, see
  // the external resolution-check tooling) show no early/late bias for tracker
  // hits. Calo hits do show a real late-side bias; see
  // cal_late_nsigma_window_asym below for that separate, asymmetric window.
  Parameter<float> trk_nsigma_window{this, "eventbuilder:unfolder:trk_nsigma_window", 3.0,
                                 "N-sigma half-width for per-detector hit acceptance"};

  // Calo HITS use their own window (cal_nsigma_window and
  // cal_late_nsigma_window_asym), not the tracker/TOF window above. Reason:
  // calo timing (shower development plus light-collection spread) is a
  // different physical process from tracker crossing time, and it has a
  // real late-side bias that the tracker case does not.
  //
  // Calo CLUSTERS do not use these two parameters. CalTimeCoincidence gates
  // clusters upstream (its own nsigma/late_nsigma_window_asym, plus a
  // minimum cluster energy) before EventUnfolder sees them; see
  // m_clu_in_names below. These two parameters apply only to the raw
  // CalorimeterHit gate further down, which has no coincidence factory of
  // its own.
  //
  // cal_late_nsigma_window_asym's default is not yet tuned; retune with an
  // efficiency/purity scan.
  Parameter<float> cal_late_nsigma_window_asym{this, "eventbuilder:unfolder:cal_late_nsigma_window_asym", 0.5,
      "extra LATE acceptance, in units of that hit's own sigma "
      "(getTimeError()), added to the + side of the calo HIT gate only "
      "(shower-development/light-collection delay; separate from "
      "CalTimeCoincidence's cluster-level version -- the tracker gate has no "
      "late-side widening at all, see the trk_nsigma_window comment above). "
      "UNMEASURED PLACEHOLDER — retune via an efficiency/purity scan."};
  Parameter<float> cal_nsigma_window{this, "eventbuilder:unfolder:cal_nsigma_window", 5.0,
      "N-sigma half-width for the calo HIT gate (separate from trk_nsigma_window "
      "above, and from CalTimeCoincidence's cluster-level nsigma)"};

  // MCParticles carried into each child. Cloning the full frame MCParticle
  // list into every candidate multiplies file size by the candidate count.
  // Default: keep only particles produced within mc_time_window of the
  // candidate's t0. This keeps the matched collision's products and any
  // real coincident background, and drops almost everything for a fake
  // candidate. Hit-to-MC truth links are unaffected: they point to the
  // original Timeslice-level collection, not to these clones. Set to -1 to
  // keep every particle, as before this parameter existed.
  Parameter<float> mc_time_window{this, "eventbuilder:truth:mc_time_window", 50.0,
      "half-width (ns) around t0 for MCParticles carried into the child; -1 = all"};

  // Scoped, temporary mitigation for a real EVGEN-cache bug (signal.sh's
  // per-class evgen directory is not beam-scoped, so a class with no evgen
  // for the requested beam silently falls back to whatever other beam's
  // files happen to be cached -- confirmed on live 10x100 production:
  // ddvcs/jpsiphoto/omega/upi0 only ever had 18x275 evgen, so their
  // contribution to "10x100" frames is actually 18x275-kinematics events
  // run through the 10x100 detector/beampipe geometry). Not a code fix for
  // that root cause (deferred) -- this only stops the exposed TRUTH LABEL
  // for these classes' hits/links from claiming a specific, wrong class.
  // Comma-separated "lo:hi" generatorStatus bands (inclusive); a leading
  // particle whose effective status falls in one is truthWeight()'d as
  // status 0 (unknown) instead of its real class band. Hits/particles
  // still simulate and reconstruct normally -- only the label is scrubbed.
  Parameter<std::string> exclude_status_ranges{
      this, "eventbuilder:truth:exclude_status_ranges", "",
      "comma-separated lo:hi generatorStatus bands scrubbed to status 0 in "
      "truth labels (mitigates cross-beam EVGEN contamination for specific "
      "classes; empty = no exclusion)"};

  // Cross-frame Si/B0 recovery. EventUnfolder fetches slow-detector hits
  // from adjacent frames through TimesliceBuffer_service and gates them per
  // candidate, the same way as current-frame hits. forward_frames adds no
  // output delay: it only waits for frames that are already in flight (see
  // the nthreads/max_inflight requirement below).
  Parameter<int32_t> backward_frames{this, "eventbuilder:unfolder:backward_frames", 1,
      "Past frames whose Si/B0 hits are included when gating each candidate. 0 = off."};
  // Default is 1. Measured hit-time residuals around t0 have a late tail,
  // so hits spill more often into the NEXT frame than into the previous
  // one; forward recovery matters at least as much as backward recovery.
  // It needs a concurrent pipeline (nthreads >= 2). When that is not
  // available, checkConfig() falls back to 0 and prints a warning.
  Parameter<int32_t> forward_frames{this, "eventbuilder:unfolder:forward_frames", 1,
      "Future frames whose Si/B0 hits are included when gating each candidate. 0 = off. "
      "Requires nthreads >= 2 and jana:max_inflight_timeslices > forward_frames; "
      "degrades to 0 with a warning when unsupported."};

  // Diagnostic for raw-hit truth gating: prints each detector's input,
  // kept, and output sizes per candidate for the tracker and raw/link/asso
  // chains. Use this to tell "the gate kept nothing" apart from "the gate
  // worked correctly": both look the same in the output file otherwise,
  // since an empty gated output causes no error.
  Parameter<bool> debug_gating{this, "eventbuilder:unfolder:debug_gating", false,
      "Print per-candidate input/kept/output sizes for the raw-hit truth gating."};

  std::shared_ptr<TimesliceBuffer_service> m_ts_buffer;
  bool m_warned_candidate_overflow = false;
  bool m_checked_config = false;
  std::vector<std::pair<int, int>> m_excluded_ranges; // parsed exclude_status_ranges, see checkConfig()
  bool isExcludedStatus(int status) const {
    for (const auto& [lo, hi] : m_excluded_ranges)
      if (status >= lo && status <= hi)
        return true;
    return false;
  }

  // Adjacent-frame slow-detector hits, fetched once per parent frame
  // (Unfold() runs sequentially, one parent at a time).
  uint64_t m_adjacent_parent = std::numeric_limits<uint64_t>::max();
  std::vector<TimesliceBuffer_service::FrameHits> m_adjacent;

  // -- tracker RecHits (10 detectors) ------------------------------------------
  // Fast detectors (TOF, MPGD): read from TrkTimeCoincidence output
  // (*TrkCoincRecHits, current frame). Already gated to candidate t0.
  // Slow detectors (Si, B0): also read from TrkTimeCoincidence output
  // (*TimeCoincRecHits, current frame, tag "coincidence"). Unfold() adds
  // adjacent-frame hits from TimesliceBuffer_service (backward_frames and
  // forward_frames) and gates those inline with the same N-sigma test, since
  // no coincidence factory runs on adjacent frames. Measurement2D hits below
  // use the same plain N-sigma gate, for the same reason.
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
  // Current frame only, for every detector. Cross-frame recovery for RecHits
  // happens upstream, in TrkTimeAlignment_factory. Raw hits carry no timing
  // correction, so cross-frame buffering does not apply to them.
  //
  // TOF uses the plain "RawHitFrame" chain here, not the "Shared"
  // (charge-sharing) chain that BTOF.cc/ECTOF.cc also define. The gate below
  // matches each raw hit to a kept RecHit through the getRawHit() relation,
  // and the kept TOF RecHits (TOFBarrel/EndcapTimeAlign) come from the plain
  // chain. Both sides must be the same chain, or no match succeeds and the
  // TOF truth output is silently empty.
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

  // RawHitAssociations are retired for both families: RawHitLink is the
  // single raw-hit truth carrier, and eicrecon's own ACTS truth-matching
  // chain (ActsToTracks -> CentralCKFTrackAssociations) consumes the links
  // directly. A link whose particle is not part of this child event gets
  // its `to` relation cleared (see the link loop), so truth matching can
  // never resolve to a particle outside the child's MCParticles.

  // -- TOF Shared-chain truth links (charge-sharing digitization) --------------
  // TOF has two digitization chains: the plain one (one raw hit per Geant4
  // deposit -- the *RawHitLinks above, the GNN's per-hit truth) and the
  // Shared one (SiliconChargeSharing splits each deposit across neighbouring
  // AC-LGAD pads; LGADHitClustering then builds the Measurement2D central
  // tracking consumes). The child's TOF measurements reference the SHARED
  // chain's raw hits, so these links are what give TOF hits ACTS truth
  // matching (CentralTrackingRawHitLinks lists them). Gated by the raw hits
  // behind this candidate's own gated Measurement2D, labelled and re-pointed
  // exactly like the tracker links. The frame inputs are optional: a
  // geometry without the charge-sharing chain simply leaves these empty and
  // the plain chain remains the only TOF truth (Shared is the default
  // whenever the frame provides it). Aligned 1:1 with m_meas_*_names
  // (index 0 = barrel, 1 = endcap).
  std::vector<std::string> m_tof_shared_link_in_names = {
      "TOFBarrelSharedRawHitLinkFrame", "TOFEndcapSharedRawHitLinkFrame"};
  std::vector<std::string> m_tof_shared_link_out_names = {
      "TOFBarrelSharedRawHitLinks", "TOFEndcapSharedRawHitLinks"};
  // Child-local copies of the Shared-chain SimTrackerHits (split deposits,
  // a different granularity than the plain *SimHits -- kept in their own
  // collections so nobody double-counts energy by summing both).
  std::vector<std::string> m_tof_shared_simhit_out_names = {
      "TOFBarrelSharedSimHits", "TOFEndcapSharedSimHits"};

  // -- calo hits + per-hit truth (11 detectors) --------------------------------
  // Same pattern as the tracker block above, one level below the clusters:
  // per-candidate, time-gated CalorimeterHits (the GNN's calo input
  // features) and their *RawHitLinks (per-hit provenance truth, below).
  // Every detector plugin provides the frame-level mirror chain
  // (…RecHitFrame / …RawHitLinkFrame, Timeslice level; see e.g. BEMC.cc).
  // EcalLumiSpec has no such mirror and is left out here.
  // Unlike the trackers, no calo *RawHitAssociations are emitted: nothing
  // consumes them in the child event (the ACTS truth chain is tracker-only,
  // and calo clusters + their particle associations pass through from the
  // frame), so cloning them was dead weight.
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
  // Calo truth links, the calo counterpart of the tracker *RawHitLinks and
  // the single per-hit truth carrier for the calo family. Same
  // label-in-weight convention, same reason: the link's `to`
  // (SimCalorimeterHit) belongs to the parent frame and is not written into
  // a child, so the label is baked into the weight rather than chased.
  //
  // edm4eic::MCRecoCalorimeterHitLink is
  // podio::LinkCollection<edm4hep::RawCalorimeterHit, edm4hep::SimCalorimeterHit>,
  // produced by CalorimeterHitDigi for all 11 detectors alongside the
  // Association form. The associations stay as plain clones (weight 1.0) for
  // anything downstream that still consumes them.
  std::vector<std::string> m_calink_in_names = {
      "EcalBarrelScFiRawHitLinkFrame",   "EcalBarrelImagingRawHitLinkFrame",
      "EcalEndcapNRawHitLinkFrame",      "EcalEndcapPRawHitLinkFrame",
      "B0ECalRawHitLinkFrame",           "HcalBarrelRawHitLinkFrame",
      "HcalEndcapNRawHitLinkFrame",      "HcalEndcapPInsertRawHitLinkFrame",
      "LFHCALRawHitLinkFrame",           "EcalFarForwardZDCRawHitLinkFrame",
      "HcalFarForwardZDCRawHitLinkFrame"};
  std::vector<std::string> m_calink_out_names = {
      "EcalBarrelScFiRawHitLinks",   "EcalBarrelImagingRawHitLinks",
      "EcalEndcapNRawHitLinks",      "EcalEndcapPRawHitLinks",
      "B0ECalRawHitLinks",           "HcalBarrelRawHitLinks",
      "HcalEndcapNRawHitLinks",      "HcalEndcapPInsertRawHitLinks",
      "LFHCALRawHitLinks",           "EcalFarForwardZDCRawHitLinks",
      "HcalFarForwardZDCRawHitLinks"};

  // -- TOF Measurement2D (LGAD cluster hits): the ACTS measurement inputs ------
  // Central tracking reads TOF as Measurement2D ("TOFBarrelClusterHits" and
  // "TOFEndcapClusterHits", from LGADHitClustering), not as TrackerHit. This
  // copies the frame-level hits into each child, gated around t0, so track
  // finding on the unfolded candidate uses them directly instead of
  // re-running digitization (which aborts on missing sim inputs; see the
  // README).
  std::vector<std::string> m_meas_in_names = {
      "TOFBarrelClusterHitFrame", "TOFEndcapClusterHitFrame"};
  std::vector<std::string> m_meas_out_names = {
      "TOFBarrelClusterHits", "TOFEndcapClusterHits"};

  // -- calo clusters + cluster->particle associations (4 detectors) ------------
  // Read from CalTimeCoincidence output (*CoincClusters). Already gated to
  // candidate t0 and to a minimum cluster energy. EcalLumiSpec stays
  // excluded, as before.
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
  // Inputs live in the parent frame (Timeslice level); outputs go to the
  // child (PhysicsEvent). Every input below sets .level = Timeslice so JANA
  // reads it from the parent, not from the (empty) child.
  PodioInput<edm4hep::EventHeader> m_candidates_in{
      this, {.name = "EventCandidatesFiltered", .level = JEventLevel::Timeslice}};

  // Debug copies of both candidate-list stages (same level issue and fix as
  // the TimeAlign/TimeCoinc pass-through below). EventCandidates is the list
  // before the prefilter; EventCandidatesFiltered is the list
  // m_candidates_in reads. Both are re-exported here as PhysicsEvent-level
  // subsets, only so they appear in the output ROOT file for inspection.
  PodioInput<edm4hep::EventHeader> m_candidates_raw_in{
      this, {.name = "EventCandidates", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_candidates_raw_out{this, "EventCandidates"};
  PodioOutput<edm4hep::EventHeader> m_candidates_filtered_out{this, "EventCandidatesFiltered"};

  // Three headers per child:
  //   EventHeader      — the child's own header. Standard edm4hep fields
  //                      only (eventNumber, runNumber, timeStamp = t0 in ps,
  //                      weight = 1). Any standard tool can read this
  //                      header safely: its weights list carries no
  //                      eventbuilder data.
  //   TimesliceHeader  — a copy of the parent frame's EventHeader, kept for
  //                      provenance.
  //   EventBuilderInfo — the eventbuilder data (t0, t0sigma, significance,
  //                      truth labels, and more; see the layout comment
  //                      where this is filled, below). Same type as
  //                      EventHeader, but a different name, so no standard
  //                      tool mistakes its weights list for MC generator
  //                      weights.
  PodioInput<edm4hep::EventHeader> m_event_header_in{
      this, {.name = "EventHeader", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_timeslice_header_out{this, "TimesliceHeader"};
  PodioOutput<edm4hep::EventHeader> m_child_header_out{this, "EventHeader"};
  PodioOutput<edm4hep::EventHeader> m_info_out{this, "EventBuilderInfo"};
  // Frame-constant configuration record (one entry per frame at Timeslice
  // level; see EventBuilder_factory::emitFrameInfo for the layout). Copied
  // onto each child with eventNumber = frame*1000, so offline readers join
  // it to EventBuilderInfo (frame*1000 + candidate) via eventNumber/1000.
  PodioInput<edm4hep::EventHeader> m_frame_info_in{
      this, {.name = "EventBuilderFrameInfo", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_frame_info_out{this, "EventBuilderFrameInfo"};
  // Frame-level GNN prefilter scores (weights[0..6]: 7-class softmax, class
  // 0 = background). EventPrefilter_factory produces these when a model
  // runs. Each child PhysicsEvent gets a copy, so the prefilter cut or
  // class tag can be applied or changed offline, with no need to reload the
  // ONNX model. Empty when the prefilter is off or has no model.
  PodioInput<edm4hep::EventHeader> m_scores_in{
      this, {.name = "PrefilterScores", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_scores_out{this, "PrefilterScores"};

  // Same pass-through as PrefilterScores above, for the two raw model
  // outputs EventPrefilter_factory stores only when
  // eventbuilder:prefilter:store_aux_outputs is on (see that file). Neither
  // lambda_k nor seg_logits is trained/meaningful yet; both stay empty
  // unless that parameter is set and the loaded model declares them.
  PodioInput<edm4hep::EventHeader> m_lambda_k_in{
      this, {.name = "PrefilterLambdaK", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_lambda_k_out{this, "PrefilterLambdaK"};
  PodioInput<edm4hep::EventHeader> m_seg_logits_in{
      this, {.name = "PrefilterSegLogits", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_seg_logits_out{this, "PrefilterSegLogits"};

  PodioInput<edm4hep::MCParticle> m_mcparticles_in{
      this, {.name = "MCParticles", .level = JEventLevel::Timeslice}};
  PodioOutput<edm4hep::MCParticle> m_mcparticles_out{this, "MCParticles"};

  VariadicPodioInput<edm4eic::TrackerHit> m_trk_in{
      this, {.names = m_trk_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::TrackerHit> m_trk_out{this, m_trk_out_names};

  // Declared BEFORE m_raw_in/m_link_in below. JANA fills a
  // component's inputs in declaration order. For an optional input, JANA
  // only checks whether the databundle already exists; it does not force
  // the upstream factory to run (see PodioInput<T>::Populate() in
  // JHasInputs.h). TOF's raw-hit truth comes from the "Shared"
  // charge-sharing chain (TOFBarrel/EndcapShared RawHit*Frame in
  // BTOF.cc/ECTOF.cc), the same chain that produces the Measurement2D input
  // below (TOFBarrel/EndcapClusterHitFrame). Declaring m_meas_in first
  // forces that chain to run, so its raw-hit outputs already exist when
  // m_raw_in/m_link_in look for them further down. With the
  // declaration order reversed, every TOF RawHit/RawHitLink
  // comes out empty, with no error.
  //
  // B0ECAL, BEMC, and EEMC are registered at Timeslice level. FEMC
  // (EcalEndcapP) has no plugin yet, so its optional inputs stay nullptr
  // (handled below).
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
  // TOF Shared-chain truth: see m_tof_shared_link_*_names.
  VariadicPodioInput<podio::Link<edm4eic::RawTrackerHit, edm4hep::SimTrackerHit>>
      m_tof_shared_link_in{
          this,
          {.names = m_tof_shared_link_in_names, .level = JEventLevel::Timeslice,
           .is_optional = true}};
  VariadicPodioOutput<podio::Link<edm4eic::RawTrackerHit, edm4hep::SimTrackerHit>>
      m_tof_shared_link_out{this, m_tof_shared_link_out_names};
  VariadicPodioOutput<edm4hep::SimTrackerHit> m_tof_shared_simhit_out{
      this, m_tof_shared_simhit_out_names};

  // Optional per-candidate SimTrackerHit copies (eventbuilder:truth:simhits).
  //
  // Off by default. The truth LABEL is always available without them: the
  // unfolder bakes generatorStatus(+slot) into each RawHitLink's weight
  // while the parent frame is still in memory, so a child event needs no
  // pointer chase to know what produced a hit.
  //
  // Turn this on when you need the FULL truth chain offline
  // (link -> SimTrackerHit -> MCParticle): exact particle kinematics,
  // energy deposits, or walking Geant4 secondaries back to their primary.
  // Without it the link's `to` side dangles, because the parent frame's
  // SimTrackerHit collection is never written into the child tree.
  //
  // Cost: sim hits are duplicated into every candidate whose window contains
  // them, and candidate windows overlap — expect a large output-size
  // increase on the collection type that already dominates volume.
  std::vector<std::string> m_simhit_out_names = {
      "TOFBarrelSimHits",       "TOFEndcapSimHits",          "MPGDBarrelSimHits",
      "OuterMPGDBarrelSimHits", "BackwardMPGDEndcapSimHits", "ForwardMPGDEndcapSimHits",
      "SiBarrelVertexSimHits",  "SiBarrelSimHits",           "SiEndcapTrackerSimHits",
      "B0TrackerSimHits"};
  VariadicPodioOutput<edm4hep::SimTrackerHit> m_simhit_out{this, m_simhit_out_names};

  // On by default. It used to be off because the only way to make the chain
  // resolve was mc_time_window=-1, which pulled the whole frame's particle
  // list into every candidate (1121 MB of MCParticles in a 1161 MB file).
  // Sim hits are now written only for particles that pass the gate, so the
  // chain resolves fully at the default window: 113 MB against 104 MB with
  // this off. That is cheap enough to be the default.
  bool m_simhits = true;
  Parameter<bool> simhits{
      this, "eventbuilder:truth:simhits", m_simhits,
      "1 (default) = also write each candidate's SimTrackerHits (*SimHits) and "
      "re-point the RawHitLinks at those child-local copies, so the full "
      "link -> SimTrackerHit -> MCParticle chain resolves offline. Only sim hits "
      "whose particle passes the mc_time_window gate are written, so everything "
      "written resolves. 0 = links carry the truth label in their weight only "
      "(smaller output; the per-hit label is identical either way)."};

  VariadicPodioInput<edm4eic::CalorimeterHit> m_calhit_in{
      this, {.names = m_calhit_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::CalorimeterHit> m_calhit_out{this, m_calhit_out_names};

  VariadicPodioInput<edm4eic::MCRecoCalorimeterHitLink> m_calink_in{
      this, {.names = m_calink_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::MCRecoCalorimeterHitLink> m_calink_out{
      this, m_calink_out_names};

  VariadicPodioInput<edm4eic::Cluster> m_clu_in{
      this, {.names = m_clu_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::Cluster> m_clu_out{this, m_clu_out_names};

  VariadicPodioInput<edm4eic::MCRecoClusterParticleAssociation> m_cluasso_in{
      this, {.names = m_cluasso_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::MCRecoClusterParticleAssociation> m_cluasso_out{
      this, m_cluasso_out_names};

  // -- TimeAlign/TimeCoinc debug pass-through -----------------------------------
  // m_trk_in/m_clu_in above exist only at Timeslice level. The PODIO writer
  // runs at PhysicsEvent level and cannot see them there: it lists
  // collections through event->GetAllCollectionNames(), which checks only
  // the current JEvent's own factories, not the parent's. This re-exports
  // them under their original names, as full (ungated) copies, so they
  // appear in the output ROOT file. Every PhysicsEvent child of the same
  // parent gets an identical copy. This is for debugging only, not for
  // production output.
  VariadicPodioOutput<edm4eic::TrackerHit> m_trk_align_out{this, m_trk_in_names};
  VariadicPodioOutput<edm4eic::Cluster> m_clu_align_out{this, m_clu_in_names};

  // Ungated slow-detector *TimeAlignRecHits, from before the coincidence
  // gate. TrkTimeAlignment_factory produces these; only
  // TimeCoincidence_factory normally reads them. m_trk_in above reads only
  // the gated *TimeCoincRecHits version for Si/B0 (indices 6-9). This copy
  // exists so the coincidence gate itself can be checked, by comparing
  // against *TimeCoincRecHits.
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

  // Gate half-width for one slow-detector (Si/B0) hit, from that hit's own
  // getTimeError(). (Not const: Parameter::operator() is non-const.)
  double slowHitWindow(double hit_sigma) {
    return double(trk_nsigma_window()) * hit_sigma;
  }

  // Fast (trigger) detectors are indices 0-5. Slow (Si/B0) are indices 6-9.
  static constexpr size_t kNumFastDet = 6;

  // Checks the forward_frames setting once. Waiting for a future frame only
  // works if that frame can be processed while this one is still in the
  // (sequential) unfold stage. When the pipeline cannot do that
  // (single-threaded), this falls back to 0 and prints a warning instead of
  // aborting: forward recovery improves physics but is not required.
  int32_t m_fwd_effective = 0;
  void checkConfig() {
    if (m_checked_config)
      return;
    m_checked_config = true;
    // Parse exclude_status_ranges ("lo:hi,lo:hi,...") once. Malformed
    // entries are skipped with a warning rather than aborting the job --
    // this is a scoped mitigation flag, not a hard config requirement.
    {
      std::stringstream ss(exclude_status_ranges());
      std::string tok;
      while (std::getline(ss, tok, ',')) {
        if (tok.empty())
          continue;
        const auto colon = tok.find(':');
        if (colon == std::string::npos) {
          jout << "[eventbuilder] exclude_status_ranges: skipping malformed entry '" << tok
               << "' (expected lo:hi)" << jendl;
          continue;
        }
        try {
          const int lo = std::stoi(tok.substr(0, colon));
          const int hi = std::stoi(tok.substr(colon + 1));
          m_excluded_ranges.emplace_back(lo, hi);
          jout << "[eventbuilder] truth label exclusion active: status " << lo << "-" << hi
               << jendl;
        } catch (const std::exception&) {
          jout << "[eventbuilder] exclude_status_ranges: skipping malformed entry '" << tok << "'"
               << jendl;
        }
      }
    }
    const int32_t n_fwd = forward_frames();
    m_fwd_effective = n_fwd;
    if (n_fwd <= 0)
      return;
    auto* app = GetApplication();
    const int nthreads = app->GetNThreads();
    // jana:max_inflight_timeslices defaults to nthreads, but JANA does not
    // register it unless the user sets it. Reading an unset parameter
    // throws, so check first and use JANA's default if it is not set.
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
    const double t0   = cand.getWeights(0); // t0, in ns
    const double t0sigma  = cand.getWeights(1); // t0 sigma, in ns
    const double phys = cand.getWeights(2); // flag: real-collision count; see EventBuilder_factory.h

    // -- true-pile-up collision-instance separation ---------------------------
    // eventbuilder:trigger:store_coincident_list (on by default) appends
    // the per-candidate trigger_classes list of (time, stream) pairs for
    // every real collision within dt of t0, at weights[25+]. When the
    // parameter is off, or the candidate predates this feature,
    // coincident_times stays empty and coincidentSlot() always returns -1:
    // unchanged behavior.
    //
    // EventBuilderInfo (below) copies cand.getWeights() verbatim, so this
    // list already reaches every child. A consumer can match a kept
    // MCParticle's own time against it offline. What is added below is the
    // same match, baked into each kept hit's tracker or calo RawHitLink
    // weight (weight = generatorStatus + (slot+1)*1e6), using the same
    // label-in-weight convention as the rest of that code.
    std::vector<double> coincident_times;
    std::vector<int> coincident_stream;
    {
      const auto& w = cand.getWeights();
      if (w.size() > 25) {
        const int n = static_cast<int>(w[25]);
        for (int k = 0; k < n && 27 + 2 * k < static_cast<int>(w.size()); ++k) {
          coincident_times.push_back(w[26 + 2 * k]);
          coincident_stream.push_back(static_cast<int>(w[27 + 2 * k]));
        }
      }
    }
    // Nearest-TIME slot assignment is the wrong primary signal: during
    // real pileup, distinct collisions routinely land within dt of each
    // other (that is what makes them "coincident" at all), so "closest in
    // time" frequently picks a different collision than the one that
    // actually produced this hit's particle -- confirmed empirically,
    // ~12% of labelled links in true-pileup candidates decoded to a slot
    // whose trigger_classes entry did not match the hit's own class.
    // truthWeight already has `status` in hand at the point it needs a
    // slot, and coincident_stream already carries each candidate
    // collision's class -- match on THAT first (never assign a slot whose
    // class disagrees with the hit's own), and use nearest time only to
    // break a tie between two same-class collisions genuinely coincident
    // in this candidate (e.g. two ncdisq1 collisions overlapping). No
    // stream match at all means this particle's own collision was not
    // independently recognized as coincident (e.g. carried in only via
    // the wider mc_time_window gate, past dt) -- -1 (no slot) is the
    // honest answer, not a guess.
    auto coincidentSlot = [&](double t, int status) -> int {
      if (coincident_times.empty())
        return -1;
      const int stream = (status == 1) ? 0 : status / 1000;
      int best = -1;
      double best_d = 0.0;
      for (size_t k = 0; k < coincident_stream.size(); ++k) {
        if (coincident_stream[k] != stream)
          continue;
        const double d = std::abs(t - coincident_times[k]);
        if (best < 0 || d < best_d) {
          best_d = d;
          best   = int(k);
        }
      }
      return best;
    };

    // The per-hit truth label written into every RawHitLink weight:
    //
    //   weight = (slot+1)*1e6 + (generatorStatus & 0xFFFF)
    //
    // Two fields, nested, never overlapping. The low half is the generator
    // status (class band + native code, <= 35999 by construction); the high
    // half is which coincident collision this hit belongs to, offset by one
    // so a high half of exactly 0 unambiguously means "no matching
    // collision in this candidate's trigger_classes" (coincidentSlot
    // returned -1), never "slot 0". Decode: high = floor(weight/1e6); slot
    // = high - 1 if high > 0, else "no slot".
    //
    // The mask is a guard, not a transformation: DDG4 stores the generator
    // status in an unsigned short, so it cannot exceed 65535 anyway. Stating
    // it here means the packing declares the invariant it depends on rather
    // than inheriting it by luck -- if a future scheme widens the status, the
    // label truncates visibly instead of silently colliding with the slot.
    //
    // slot stays below 16: float32 holds integers exactly only to 2^24, and
    // 16*1e6 + 65535 is the last value that fits.
    checkConfig(); // idempotent; guarantees m_excluded_ranges is parsed before use below
    auto truthWeight = [&](int status, double ptime) -> float {
      // Slot is matched against the hit's TRUE class (coincidentSlot's own
      // stream match), before any scrubbing below -- an excluded particle
      // still really did come from a specific coincident collision, and
      // zeroing status first would make it look native-primary (stream 0)
      // instead of correctly finding no match.
      const int slot = coincidentSlot(ptime, status);
      // Scrub the label for a scoped, temporary EVGEN-contamination
      // mitigation -- see exclude_status_ranges's declaration comment.
      if (isExcludedStatus(status))
        status = 0;
      float w = static_cast<float>(status & 0xFFFF);
      // Packed as slot+1, not slot: 0 must mean "no slot" unambiguously.
      // With true pileup, coincidentSlot can now legitimately return -1 for
      // ONE hit in a candidate while another hit in the SAME candidate
      // matches slot 0 -- unlike before this session's fix, when -1 only
      // ever happened candidate-wide (empty trigger_classes), so slot 0
      // and "no slot" never needed to coexist in one candidate's hits.
      // Packing slot 0 as bare 0*1e6 would make it indistinguishable from
      // "no slot" on a per-hit basis; +1 removes the collision.
      assert(slot < 16 && "coincident slot exceeds float32-exact packing range");
      if (slot >= 0)
        w += static_cast<float>(slot + 1) * 1e6f;
      return w;
    };

    // -- adjacent-frame Si/B0 hits (cross-frame boundary recovery) -------------
    // Fetched once per parent from TimesliceBuffer_service. Keyed by
    // absolute frame number, so this works regardless of the order in which
    // frames were processed upstream.
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
    // Fast detectors: a subset of the parent's current-frame collections,
    // windowed per-hit by that hit's own getTimeError(). Slow detectors:
    // owned clones, merging the current frame with the adjacent frames
    // above (hits from other frames cannot be subset references), windowed
    // the same way.
    //
    // Per-detector set of the RawTrackerHits behind this candidate's own
    // gated RecHits below. Used to gate the raw hits/links/associations
    // further down to just this candidate, instead of the whole frame; see
    // copyGatedByRawHit's comment.
    std::vector<std::unordered_set<uint64_t>> kept_raw(m_trk_in().size());

    // Adjacent-frame hits kept by this candidate, per detector, paired with
    // their buffered source. Their raw hit and truth link are rebuilt from
    // the buffer after the current-frame copies below, because their
    // originals died with the frame that produced them.
    std::vector<std::vector<std::pair<edm4eic::MutableTrackerHit,
                                      const TimesliceBuffer_service::SlowHit*>>>
        adj_rebuild(m_trk_in().size());

    for (size_t det = 0; det < m_trk_in().size(); ++det) {
      auto& out = m_trk_out().at(det);
      const auto* in = m_trk_in().at(det);

      if (det < kNumFastDet) {
        // Owned clones, not a subset: fast-detector inputs are Timeslice-level
        // collections that are never written to the child tree. A subset
        // output would only store an unreadable <name>_objIdx reference to
        // them.
        if (in == nullptr)
          continue;
        // in is already gated to candidate t0 by TrkTimeCoincidence (see
        // m_trk_in_names above); clone through, with no re-gating here.
        for (const auto& hit : *in) {
          out->push_back(hit.clone(true));
          kept_raw[det].insert(rawHitKey(hit.getRawHit()));
        }
        continue;
      }

      // (hit, source) pairs: source is null for current-frame hits, whose
      // raw hit / link / association copies are handled by the gated
      // whole-collection copies below. Adjacent-frame hits carry a pointer
      // to their buffered SlowHit instead, and get their own raw hit and
      // link rebuilt from it further down -- see adj_rebuild.
      std::vector<std::pair<edm4eic::MutableTrackerHit,
                            const TimesliceBuffer_service::SlowHit*>> merged;
      if (in != nullptr)
        for (const auto& hit : *in) {
          const double sigma = double(hit.getTimeError());
          const double w = slowHitWindow(sigma);
          const double dt = hit.getTime() - t0;
          if (dt >= -w && dt <= w) {
            merged.emplace_back(hit.clone(), nullptr);
            kept_raw[det].insert(rawHitKey(hit.getRawHit()));
          }
        }
      // Adjacent-frame hits are cloned WITHOUT their relations
      // (cloneRelations=false), unlike the current-frame clone above. A
      // default hit.clone() keeps the relation pointing at the original
      // RawTrackerHit, which here lives in one of the m_adjacent frames.
      // m_adjacent is cleared and replaced (see above) as soon as the
      // unfolder moves to the next parent timeslice. If this child's
      // TrackerHit is still unwritten at that point, prepareForWrite()
      // later dereferences an already-destroyed RawTrackerHit and crashes.
      // The current-frame clone above does not have this problem: its
      // source collection (m_trk_in(), owned by the parent event) always
      // outlives this child.
      const size_t si = det - kNumFastDet;
      for (const auto& frame : m_adjacent)
        if (si < frame.size())
          for (const auto& sh : frame[si]) {
            const double sigma = double(sh.hit.getTimeError());
            const double w = slowHitWindow(sigma);
            const double dt = sh.hit.getTime() - t0;
            if (dt >= -w && dt <= w)
              merged.emplace_back(sh.hit.clone(false), &sh);
          }
      std::sort(merged.begin(), merged.end(),
                [](const auto& a, const auto& b) {
                  return a.first.getTime() < b.first.getTime();
                });
      for (auto& [hit, src] : merged) {
        out->push_back(hit);
        if (src != nullptr)
          adj_rebuild[det].emplace_back(hit, src);
      }
    }

    // -- TOF Measurement2D (LGAD cluster hits): gated like the fast hits -------
    // Per-measurement time error is sqrt(covariance.zz) (loc0/loc1/time
    // ordering). Falls back to the 25 ps LGAD resolution when unset.
    // kept_shared collects the Shared-chain RawTrackerHits behind each kept
    // measurement's own hits; it gates the TOF Shared truth links below the
    // same way kept_raw gates the tracker links.
    std::vector<std::unordered_set<uint64_t>> kept_shared(m_meas_in().size());
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
        // Clone, not a subset: same reason as the fast detectors above.
        if (dt >= -win && dt <= win) {
          out->push_back(m.clone(true));
          for (const auto& h : m.getHits())
            if (h.getRawHit().isAvailable())
              kept_shared[d].insert(rawHitKey(h.getRawHit()));
        }
      }
    }

    // The MCParticle gate, computed here because the sim-hit copies below
    // have to obey it too. Same value the MCParticle clone loop uses further
    // down; see its comment for why dt is a floor.
    float mc_win = mc_time_window();
    if (mc_win >= 0.0f && m_frame_info_in() != nullptr && m_frame_info_in()->size() > 0) {
      const auto& fi = m_frame_info_in()->at(0).getWeights();
      if (fi.size() > 3 && static_cast<float>(fi[3]) > mc_win)
        mc_win = static_cast<float>(fi[3]);
    }
    // True when a particle belongs to THIS candidate, by the same rule the
    // MCParticle clone loop applies.
    auto particleInGate = [&](double ptime) {
      return mc_win < 0.0f || std::abs(ptime - t0) <= mc_win;
    };

    // Child-local SimTrackerHit copies (eventbuilder:truth:simhits).
    // Keyed by the parent hit's (collectionID, index) so one sim hit backing
    // several kept raw hits is written once. The sim -> MCParticle relation
    // cannot be set here: the MCParticle clones (mcp_map) are made further
    // below, so each new hit records the parent particle index and is fixed
    // up right after that loop.
    //
    // A sim hit is written only when its MCParticle passes the gate, so the
    // child stays self-consistent: everything written belongs to this
    // candidate, and every sim hit resolves to a particle that is also
    // present. The alternative -- writing every sim hit and letting the gated
    // ones keep an unset relation -- produced collections that were mostly
    // unresolvable (about 90% for TOF/MPGD, whose particles are typically
    // produced hundreds of ns from t0), and widening the gate to fix that
    // pulled the whole frame's particle list into every candidate: 1.1 GB of
    // MCParticles per file to make 2.4% of them reachable.
    //
    // Ground truth does not depend on any of this. The per-hit label lives in
    // the RawHitLink weight and is written for every kept hit whether or not
    // its particle belongs to this candidate.
    std::unordered_map<uint64_t, edm4hep::MutableSimTrackerHit> simhit_map;
    std::vector<std::pair<edm4hep::MutableSimTrackerHit, uint32_t>> simhit_fixups;

    // -- raw hits / links / hit associations: gated to THIS candidate ----------
    // Gated per candidate, not copied whole-frame: see copyGatedByRawHit's
    // comment for why an ungated copy is wrong for a per-candidate consumer
    // (for example, GNN segmentation truth).
    copyGatedByRawHit(m_raw_in(), m_raw_out(), kept_raw,
                       [](const auto& h) -> const edm4eic::RawTrackerHit& { return h; });
    // Tracker RawHitAssociations are cloned in their own loop AFTER the
    // links loop below (they need its simhit_map) — see that loop's comment.
    // Links use an explicit loop instead of copyGatedByRawHit, because the
    // clone's weight is overwritten with the sim particle's
    // generatorStatus. Reason: the link's `to` (SimTrackerHit) relation
    // points into the parent frame's SimTrackerHit collection, which is
    // never written to the child tree. The link->simHit->MCParticle chain
    // that turns a hit into a signal/background label cannot be resolved
    // offline, so the label is baked into the weight here instead, while
    // the whole chain is still in memory. (Formerly this lived on the
    // tracker RawHitAssociations, now retired — RawHitLink is the single
    // truth-label carrier for the tracker family.)
    //
    // generatorStatus encodes provenance (assigned at generation time):
    //   < 10000   machine background / native (bkg sub-bands at 2000-6999)
    //   >= 10000  physics classes, one 1000-wide band per class (base =
    //             10000 + class_index*1000, max 35999). A wider,
    //             10000-wide-per-class layout with 1000-wide same-class
    //             pile-up instance sub-bands is designed (rewrite_status.py
    //             in the production repo) but not live: SBM's merged output
    //             format blocks that rewrite today, so every real status is
    //             instance 1 of the 1000-wide layout above.
    // The original weight is always 1.0 (from SiliconTrackerDigi), so
    // nothing is lost by overwriting it. Offline, join a link to a RecHit
    // by rawHit (`from`) ObjectID (collectionID, index) equality; the
    // clone keeps that relation.
    for (size_t i = 0; i < m_link_out().size(); ++i) {
      if (i >= m_link_in().size() || m_link_in().at(i) == nullptr || i >= kept_raw.size())
        continue;
      const auto& keep = kept_raw[i];
      for (const auto& l : *m_link_in().at(i)) {
        if (keep.count(rawHitKey(l.getFrom())) == 0)
          continue;
        auto clone = l.clone(true);
        const auto sim = l.getTo();
        // Give this child its own SimTrackerHit copy and re-point the link
        // at it, so the truth chain resolves offline AND ActsToTracks (which
        // walks link -> simHit -> particle) lands on the child's own
        // MCParticles. Deduplicated per detector — one sim hit can back
        // several kept raw hits. When no child-local copy exists (simhits
        // off, or the particle is outside the mc_time_window gate), the `to`
        // relation is CLEARED instead of left pointing at the parent frame:
        // a frame-level ref serializes as (name-hash, frame index) and would
        // silently mis-resolve against the child's same-named collections —
        // the exact corruption that once made 37.6% of CKF association
        // targets decode as neutral particles. The weight below carries the
        // per-hit label either way.
        bool repointed = false;
        if (simhits() && sim.isAvailable() && i < m_simhit_out().size() &&
            sim.getParticle().isAvailable() &&
            particleInGate(double(sim.getParticle().getTime()))) {
          const uint64_t skey = (static_cast<uint64_t>(sim.getObjectID().collectionID) << 32) |
                                static_cast<uint32_t>(sim.getObjectID().index);
          auto it = simhit_map.find(skey);
          if (it == simhit_map.end()) {
            auto newSim = sim.clone(false); // relations re-set below, not carried
            if (sim.getParticle().isAvailable())
              simhit_fixups.emplace_back(
                  newSim, static_cast<uint32_t>(sim.getParticle().getObjectID().index));
            m_simhit_out().at(i)->push_back(newSim);
            it = simhit_map.emplace(skey, newSim).first;
          }
          clone.setTo(it->second);
          repointed = true;
        }
        if (!repointed)
          clone.setTo(edm4hep::SimTrackerHit::makeEmpty());
        if (sim.isAvailable() && sim.getParticle().isAvailable()) {
          // Secondaries inherit their first labelled ancestor's
          // (status, time) -- see effectiveProvenance(). The weight is the
          // provenance label; the sim -> particle chain written above still
          // points at the actual (possibly status-0) particle.
          const auto [est, etime] = effectiveProvenance(sim.getParticle());
          clone.setWeight(truthWeight(est, etime));
        }
        m_link_out().at(i)->push_back(clone);
      }
    }

    // -- TOF Shared-chain truth links: gated to THIS candidate ------------------
    // Same treatment as the tracker links above, but gated by kept_shared
    // (the raw hits behind this candidate's own gated Measurement2D) and
    // deduplicating sim copies into the *SharedSimHits collections. These
    // links are what let ActsToTracks truth-match TOF measurements; the
    // plain-chain TOF links above stay the GNN's per-hit truth. See
    // m_tof_shared_link_in_names for the two-chain picture.
    std::unordered_map<uint64_t, edm4hep::MutableSimTrackerHit> shared_simhit_map;
    for (size_t i = 0; i < m_tof_shared_link_out().size(); ++i) {
      if (i >= m_tof_shared_link_in().size() || m_tof_shared_link_in().at(i) == nullptr ||
          i >= kept_shared.size())
        continue;
      const auto& keep = kept_shared[i];
      for (const auto& l : *m_tof_shared_link_in().at(i)) {
        if (keep.count(rawHitKey(l.getFrom())) == 0)
          continue;
        auto clone = l.clone(true);
        const auto sim = l.getTo();
        bool repointed = false;
        if (simhits() && sim.isAvailable() && i < m_tof_shared_simhit_out().size() &&
            sim.getParticle().isAvailable() &&
            particleInGate(double(sim.getParticle().getTime()))) {
          const uint64_t skey = (static_cast<uint64_t>(sim.getObjectID().collectionID) << 32) |
                                static_cast<uint32_t>(sim.getObjectID().index);
          auto it = shared_simhit_map.find(skey);
          if (it == shared_simhit_map.end()) {
            auto newSim = sim.clone(false); // relations re-set below, not carried
            if (sim.getParticle().isAvailable())
              simhit_fixups.emplace_back(
                  newSim, static_cast<uint32_t>(sim.getParticle().getObjectID().index));
            m_tof_shared_simhit_out().at(i)->push_back(newSim);
            it = shared_simhit_map.emplace(skey, newSim).first;
          }
          clone.setTo(it->second);
          repointed = true;
        }
        if (!repointed)
          clone.setTo(edm4hep::SimTrackerHit::makeEmpty());
        if (sim.isAvailable() && sim.getParticle().isAvailable()) {
          const auto [est, etime] = effectiveProvenance(sim.getParticle());
          clone.setWeight(truthWeight(est, etime));
        }
        m_tof_shared_link_out().at(i)->push_back(clone);
      }
    }

    // -- adjacent-frame hits: rebuild raw hit + truth link from the buffer ----
    // Cross-frame recovery hands the unfolder a bare RecHit: its rawHit
    // relation, and the link chain behind it, belong to a frame that is
    // already gone (see the cloneRelations=false comment above). Left alone,
    // these hits reach the child event with rawHit index -1 and no link, so
    // they carry no class label at all -- 20-30% of silicon hits at low
    // background, ~50-67% at high background, and the fraction moves with
    // forward_frames, which made it look like a background-correlated
    // detector dropout. TimesliceBuffer_service resolves the truth at
    // deposit time instead and carries it here, so the child gets a real
    // RawTrackerHit to point at plus a link with the usual
    // generatorStatus(+slot) weight. Offline decoding is then identical for
    // current-frame and adjacent-frame hits.
    //
    // With eventbuilder:truth:simhits on, the sim hit and its MCParticle
    // are buffered too, so the full link -> SimTrackerHit -> MCParticle chain
    // resolves for adjacent-frame hits as well. Both objects belong to a frame
    // that no longer exists, so the child gets its own copies: the sim hit
    // joins this detector's *SimHits collection and the particle is appended
    // to the child's MCParticles. Appending is safe because a neighbouring
    // frame's MCParticles are disjoint from this frame's, so they cannot
    // collide with the gated clones made further below.
    //
    // The particle is cloned without relations: its parents and daughters
    // live in the dead frame. Kinematics, PDG and generatorStatus survive;
    // lineage does not.
    // Keyed by (source frame, object index): podio indices restart every
    // frame, so keying on the index alone silently attaches a backward
    // neighbour's particle to a forward neighbour's hit.
    std::map<std::pair<uint64_t, uint64_t>, edm4hep::MutableSimTrackerHit> adj_sim;
    std::map<std::pair<uint64_t, uint32_t>, edm4hep::MutableMCParticle>    adj_mcp;
    for (size_t det = 0; det < adj_rebuild.size(); ++det) {
      if (det >= m_raw_out().size() || det >= m_link_out().size())
        continue;
      for (auto& [hit, src] : adj_rebuild[det]) {
        if (!src->has_raw)
          continue;
        auto raw = src->raw.clone(false);
        m_raw_out().at(det)->push_back(raw);
        hit.setRawHit(raw);
        if (!src->has_truth)
          continue;
        auto link = m_link_out().at(det)->create();
        link.setFrom(raw);
        link.setWeight(truthWeight(static_cast<int>(src->status), double(src->ptime)));

        // Same gate as the current-frame path: only particles belonging to
        // this candidate are written, so the child stays self-consistent.
        // The link weight above already carries this hit's label either way.
        if (!simhits() || !src->has_sim || det >= m_simhit_out().size() ||
            !src->has_particle || !particleInGate(double(src->ptime)))
          continue;
        const auto skey = std::make_pair(src->frame, src->sim_key);
        auto sit = adj_sim.find(skey);
        if (sit == adj_sim.end()) {
          auto newSim = src->sim.clone(false);
          if (src->has_particle) {
            const auto pkey = std::make_pair(src->frame, src->mcp_key);
            auto pit = adj_mcp.find(pkey);
            if (pit == adj_mcp.end()) {
              auto newMcp = src->particle.clone(false);
              m_mcparticles_out()->push_back(newMcp);
              pit = adj_mcp.emplace(pkey, newMcp).first;
            }
            newSim.setParticle(pit->second);
          }
          m_simhit_out().at(det)->push_back(newSim);
          sit = adj_sim.emplace(skey, newSim).first;
        }
        link.setTo(sit->second);
      }
    }

    // -- calo hits + per-hit truth: gated to THIS candidate --------------------
    // Mirrors the tracker block above, at the hit level. Hits are time-gated
    // inline, with cal_nsigma_window/cal_late_nsigma_window_asym, since no
    // coincidence factory exists for raw CalorimeterHit. Clusters (further
    // down) are a different channel, gated upstream by CalTimeCoincidence
    // instead; see cal_nsigma_window's comment above. Associations are keyed
    // on the gated hits' rawHit ObjectIDs and cloned with the same
    // label-in-weight convention as the tracker links.
    //
    // One calo-specific case: a SimCalorimeterHit can mix energy deposits
    // from several particles in one cell (contributions). The label is
    // taken from the highest-energy contribution's particle, since that
    // particle dominates the cell.
    std::vector<std::unordered_set<uint64_t>> kept_calraw(m_calhit_in().size());
    for (size_t d = 0; d < m_calhit_in().size(); ++d) {
      auto& out = m_calhit_out().at(d);
      const auto* in = m_calhit_in().at(d);
      if (in == nullptr)
        continue;
      for (const auto& hit : *in) {
        const double terr  = double(hit.getTimeError());
        // Guard sigma <= 0: fall back to 1 ns, a calo-scale width, instead
        // of gating the whole family away.
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
    // Calo truth links: THE per-hit calo label, same convention as the
    // tracker *RawHitLinks. A SimCalorimeterHit can mix deposits from several
    // particles in one cell, so the label is taken from the highest-energy
    // contribution's particle, which dominates the cell.
    for (size_t i = 0; i < m_calink_out().size(); ++i) {
      if (i >= m_calink_in().size() || m_calink_in().at(i) == nullptr || i >= kept_calraw.size())
        continue;
      const auto& keep = kept_calraw[i];
      for (const auto& l : *m_calink_in().at(i)) {
        if (keep.count(rawHitKey(l.getFrom())) == 0)
          continue;
        auto clone = l.clone(true);
        const auto sim = l.getTo();
        if (sim.isAvailable()) {
          float best_e = -1.0F;
          int32_t status = 0;
          double leading_t = 0.0;
          bool found = false;
          for (const auto& c : sim.getContributions()) {
            if (c.getEnergy() > best_e && c.getParticle().isAvailable()) {
              best_e = c.getEnergy();
              // Secondaries inherit their first labelled ancestor's
              // (status, time) -- see effectiveProvenance(). Status and time
              // always describe the same particle, so slot and label agree.
              std::tie(status, leading_t) = effectiveProvenance(c.getParticle());
              found  = true;
            }
          }
          if (found) {
            clone.setWeight(truthWeight(status, leading_t));
          }
        }
        m_calink_out().at(i)->push_back(clone);
      }
    }

    if (debug_gating()) {
      for (size_t d = 0; d < m_trk_in_names.size(); ++d) {
        const auto* tin = d < m_trk_in().size() ? m_trk_in().at(d) : nullptr;
        const auto* rin = d < m_raw_in().size() ? m_raw_in().at(d) : nullptr;
        const auto* lin = d < m_link_in().size() ? m_link_in().at(d) : nullptr;
        // kept_raw[d].size() is omitted: it always equals raw_out below,
        // so printing both would be redundant.
        jout << "[eb:gating] cand=" << child_idx << " det=" << m_trk_in_names[d]
             << " trk_in=" << (tin ? std::to_string(tin->size()) : "null")
             << " trk_out=" << m_trk_out().at(d)->size()
             << " raw_in=" << (rin ? std::to_string(rin->size()) : "null")
             << " raw_out=" << m_raw_out().at(d)->size()
             << " link_in=" << (lin ? std::to_string(lin->size()) : "null")
             << " link_out=" << m_link_out().at(d)->size() << jendl;
      }
      for (size_t d = 0; d < m_calhit_in_names.size(); ++d) {
        const auto* hin = d < m_calhit_in().size() ? m_calhit_in().at(d) : nullptr;
        const auto* lin = d < m_calink_in().size() ? m_calink_in().at(d) : nullptr;
        // kept_calraw[d].size() is omitted: it always equals hit_out, so
        // printing both would be redundant.
        jout << "[eb:gating] cand=" << child_idx << " det=" << m_calhit_in_names[d]
             << " hit_in=" << (hin ? std::to_string(hin->size()) : "null")
             << " hit_out=" << m_calhit_out().at(d)->size()
             << " link_in=" << (lin ? std::to_string(lin->size()) : "null")
             << " link_out=" << m_calink_out().at(d)->size() << jendl;
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
    // Runs before the calo block below, which needs this loop's
    // old-index-to-local-clone map to remap MCRecoClusterParticleAssociation.sim.
    // Floor the gate at this frame's own coincidence window dt
    // (EventBuilderFrameInfo weights[3]). flag/trigger_classes count every
    // real collision within +-dt of t0, so a narrower particle gate would
    // drop the products of a collision the candidate is credited with. dt is
    // derived per frame from the measured hit resolution, while
    // mc_time_window is a fixed number — without this floor the invariant
    // holds only by coincidence for the current detector mix.
    // mc_win / particleInGate are computed once, further up, because the
    // sim-hit copies obey the same gate.
    std::unordered_map<uint32_t, edm4hep::MutableMCParticle> mcp_map;
    for (const auto& mcp : *m_mcparticles_in()) {
      if (!particleInGate(double(mcp.getTime())))
        continue;
      auto newMcp = mcp.clone(true);
      m_mcparticles_out()->push_back(newMcp);
      mcp_map.emplace(static_cast<uint32_t>(mcp.getObjectID().index), newMcp);
    }

    // Now that this child's MCParticles exist, complete the sim-hit chain
    // (simhits only). A hit whose particle fell outside the
    // mc_time_window gate keeps an unset relation rather than a dangling one.
    for (auto& [newSim, parent_idx] : simhit_fixups) {
      auto mit = mcp_map.find(parent_idx);
      if (mit != mcp_map.end())
        newSim.setParticle(mit->second);
    }

    // -- calo clusters + matched cluster->particle associations ----------------
    // Owned clones, not subset references. Unlike tracks, which are
    // recomputed per candidate by the ACTS chain, calo clusters are computed
    // once at the Timeslice/frame level. A subset output would store only
    // an ObjectID reference into that frame-level collection, so
    // .energy/.position/.time and so on would never appear in the child
    // PhysicsEvent tree, and downstream calo efficiency/purity checks (which
    // read the child tree alone) would see nothing. This clones instead,
    // the same way as MCParticles above. cin is already gated to candidate
    // t0 (and to a minimum energy) by CalTimeCoincidence; see
    // m_clu_in_names above. No re-gating happens here.
    for (size_t cd = 0; cd < m_clu_in().size(); ++cd) {
      auto& cout = m_clu_out().at(cd);
      auto& aout = m_cluasso_out().at(cd);
      const auto* cin = m_clu_in().at(cd);
      if (cin == nullptr) continue; // optional input; e.g. FEMC has no plugin yet
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
    // Deterministic numbering: with nthreads > 1, parent frames can reach the
    // unfolder out of order, so a running counter would not be reproducible
    // between runs. This encodes (frame, candidate) instead, and assumes
    // fewer than 1000 candidates per frame. Warn loudly if that assumption
    // ever breaks, since a 1000th candidate would collide with the next
    // frame's numbering with no other sign of error.
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
    hdr.setWeight(1.0); // standard event weight; eventbuilder data is in EventBuilderInfo instead
    m_child_header_out()->push_back(hdr);

    // EventBuilderInfo: copies the candidate's full weights list verbatim
    // (see the layout comment above EventBuilder_factory::emitCandidate).
    // Nothing is appended: the old trailing candidate-count entry was
    // retired (redundant with the eventNumber = frame*1000 + candidate
    // encoding every consumer already derives it from).
    edm4hep::MutableEventHeader info;
    info.setRunNumber(parent.GetRunNumber());
    info.setEventNumber(child_event_nr);
    info.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    for (const auto w : cand.getWeights())
      info.addToWeights(w);
    m_info_out()->push_back(info);

    // EventBuilderFrameInfo: the parent frame's frame-constant record (see
    // EventBuilder_factory::emitFrameInfo for the layout), copied onto this
    // child with eventNumber = frame*1000 so offline readers recover the
    // join key as eventNumber/1000 on both sides.
    if (m_frame_info_in() != nullptr && m_frame_info_in()->size() > 0) {
      edm4hep::MutableEventHeader fi;
      fi.setRunNumber(parent.GetRunNumber());
      fi.setEventNumber(parent_nr * 1000);
      for (const auto w : m_frame_info_in()->at(0).getWeights())
        fi.addToWeights(w);
      m_frame_info_out()->push_back(fi);
    }

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

    // Same copy, for the two raw aux outputs (see m_lambda_k_in's comment
    // above); empty unless eventbuilder:prefilter:store_aux_outputs is on.
    if (m_lambda_k_in() != nullptr && m_lambda_k_in()->size() > 0) {
      edm4hep::MutableEventHeader lk;
      lk.setRunNumber(parent.GetRunNumber());
      lk.setEventNumber(child_event_nr);
      for (const auto w : m_lambda_k_in()->at(0).getWeights())
        lk.addToWeights(w);
      m_lambda_k_out()->push_back(lk);
    }
    if (m_seg_logits_in() != nullptr && m_seg_logits_in()->size() > 0) {
      edm4hep::MutableEventHeader sl;
      sl.setRunNumber(parent.GetRunNumber());
      sl.setEventNumber(child_event_nr);
      for (const auto w : m_seg_logits_in()->at(0).getWeights())
        sl.addToWeights(w);
      m_seg_logits_out()->push_back(sl);
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

  // Identity key for a raw hit (collectionID plus index), for the kept-hit
  // lookup below. Raw hits carry no calibrated time, so they cannot be
  // gated against t0 directly, the way RecHits are. This keys off which raw
  // hit each surviving, already-gated RecHit points back to instead. Works
  // for any podio object with getObjectID: tracker uses
  // edm4eic::RawTrackerHit, calo uses edm4hep::RawCalorimeterHit.
  template <typename ObjT>
  static uint64_t rawHitKey(const ObjT& h) {
    const auto id = h.getObjectID();
    return (static_cast<uint64_t>(static_cast<uint32_t>(id.collectionID)) << 32) |
           static_cast<uint32_t>(id.index);
  }

  // Like copyAll, but keeps only elements whose raw hit (found through
  // rawHitOf) is one of this candidate's own gated RecHits (kept.at(i)).
  // Without this gate, every candidate unfolded from one frame would get
  // the frame's full raw-hit/truth population, regardless of its own t0.
  // This matters for any per-candidate consumer, for example GNN
  // segmentation truth. Raw hits, links, and associations carry no timing
  // correction of their own, so this rawHit match is the only gate
  // available at this level.
  //
  // Owned clones, not a subset: a subset output of a Timeslice-level
  // collection would only store a bare <name>_objIdx reference, since the
  // referenced frame-level collection is never written to the child tree.
  // clone(true) keeps the relations (rawHit/simHit/from/to), which point at
  // parent-frame collections that always outlive this child.
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
      "B0TrackerRecHitFrame",
      // Trailing, non-variadic inputs of TrkTimeAlignment_factory: the truth
      // links of the four slow detectors, in slow-index order. They let the
      // factory resolve each buffered Si/B0 hit's generatorStatus while its
      // producing frame is still alive, so cross-frame-recovered hits reach
      // the child event with a real label instead of a dangling rawHit.
      // Must stay last and stay exactly kNumSlowLink long.
      "SiBarrelVertexRawHitLinkFrame",
      "SiBarrelRawHitLinkFrame",
      "SiEndcapTrackerRawHitLinkFrame",
      "B0TrackerRawHitLinkFrame"
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
      // Use the ScFi clusters, not the merged "EcalBarrelClusterFrame": the
      // merger (EnergyPositionClusterMerger, ScFi x Imaging pair matching)
      // produces very few clusters at Timeslice level, so barrel events
      // would lose their E_calo tag. The ScFi side carries the EM energy
      // and clusters reliably. The TimeAlign output name
      // (EcalBarrelTimeAlignClusters) does not change, so every consumer
      // (E_calo weights[24], the tag bitmask, the unfolder) reads this as
      // before.
      "EcalBarrelScFiClusterFrame",
      "EcalEndcapNClusterFrame",
      "EcalEndcapPClusterFrame",
      "HcalBarrelClusterFrame",
      "HcalEndcapPInsertClusterFrame",
      "LFHCALClusterFrame",
      "EcalFarForwardZDCClusterFrame",
      "HcalFarForwardZDCClusterFrame",
      // Placeholders for the four slow-detector link inputs the shared
      // template declares. The calo instantiation deposits nothing, so
      // these stay empty; the inputs are optional and read back as nullptr.
      "", "", "", ""
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

  // Process-wide services: a thread-safe, frame-numbered cross-frame Si/B0
  // hit buffer, and the shared ONNX runtime used by the prefilter.
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

  // EventBuilder_factory: finds candidate physics events in a time frame and
  // assigns each a (t0, t0sigma). Runs at the Timeslice (time-frame) level.
  // Only fast detectors (TOF, MPGD, indices 0-5) take part in event finding
  // (kNumTriggerDet=6 in EventBuilder_factory.h). Si/B0 are wired as optional
  // inputs but are never read there. Cross-frame recovery of slow RecHits
  // happens downstream, in EventUnfolder. Raw hits use the current frame
  // only.
  //
  // This flat tag list is matched POSITIONALLY to the factory's declared
  // inputs: the 10-slot tracker variadic, then MCParticles, then the 9-slot
  // calo-cluster variadic (feeds the E_calo score, weights[24]). Adding an
  // input to EventBuilder_factory.h without extending this list shifts every
  // later binding, and any input past the shift silently becomes nullptr.
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
          .m_default_output_tags = {"EventCandidates", "EventBuilderFrameInfo"},
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // TrkTimeCoincidence (TimeCoincidence_factory<edm4eic::TrackerHit>): gates
  // slow (Si/B0) RecHits to candidate t0. Gates the current frame's
  // *TimeAlignRecHits against each candidate's t0 window (N sigma times the
  // hit's own getTimeError()) and produces a subset collection for
  // EventUnfolder. Cross-frame recovery happens downstream: EventUnfolder
  // merges in adjacent frames' slow hits from TimesliceBuffer_service.
  //
  // Input: EventCandidates plus 4 slow *TimeAlignRecHits (current frame,
  // Si/B0). Output: 4 *TimeCoincRecHits (subset), read by EventUnfolder's
  // slow RecHit inputs.
  const std::vector<std::string> coinc_inputs = {
      "EventCandidates",
      "SiBarrelVertexTimeAlignRecHits", "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits", "B0TrackerTimeAlignRecHits"
  };
  const std::vector<std::string> coinc_outputs = {
      "SiBarrelVertexTimeCoincRecHits", "SiBarrelTrackerTimeCoincRecHits",
      "SiEndcapTrackerTimeCoincRecHits", "B0TrackerTimeCoincRecHits"
  };
  app->Add(new JOmniFactoryGeneratorT<TrkTimeCoincidence>(
      JOmniFactoryGeneratorT<TrkTimeCoincidence>::TypedWiring{
          .m_tag                 = "coincidence",
          .m_default_input_tags  = coinc_inputs,
          .m_default_output_tags = coinc_outputs,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // TrkTimeCoincidence again: a second, separate instance of the same gate,
  // applied to the fast (TOF+MPGD) hits instead of slow (Si/B0). Keeps only
  // the fast hits that fall inside some candidate's window. Feeds
  // EventUnfolder's fast-detector inputs directly (m_trk_in indices 0-5,
  // *TrkCoincRecHits).
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
  app->Add(new JOmniFactoryGeneratorT<TrkTimeCoincidence>(
      JOmniFactoryGeneratorT<TrkTimeCoincidence>::TypedWiring{
          .m_tag                 = "trkcoincidence",
          .m_default_input_tags  = trkcoinc_inputs,
          .m_default_output_tags = trkcoinc_outputs,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // CalTimeCoincidence (TimeCoincidence_factory<edm4eic::Cluster>): the same
  // candidate-t0 gate, applied to calo clusters, plus an energy floor (off
  // by default). Feeds EventUnfolder's m_clu_in directly (*CoincClusters),
  // already gated to candidate t0 and to a minimum energy. 10 names, matching
  // EventBuilder_factory's own m_calo_collection_names. EcalLumiSpec has no
  // plugin yet and resolves to nullptr.
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
  app->Add(new JOmniFactoryGeneratorT<CalTimeCoincidence>(
      JOmniFactoryGeneratorT<CalTimeCoincidence>::TypedWiring{
          .m_tag                 = "calcoincidence",
          .m_default_input_tags  = calcoinc_inputs,
          .m_default_output_tags = calcoinc_outputs,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // EventPrefilter_factory: optional GNN background rejection between
  // EventBuilder and EventUnfolder. When eventbuilder:onnx_model is empty
  // (the default), every EventCandidate is forwarded unchanged and no ONNX
  // session is created. When a model path is set, only candidates with
  // probs[BKG=0] < eventbuilder:score_threshold survive.
  //
  // Hit/cluster features per candidate frame: [x/4000, y/4000, z/5000, t/50,
  // eDep-or-energy/10, det_id/6]. This list must stay in sync with
  // EventPrefilter_factory.h's own m_trk_hit_names/m_calo_cluster_names.
  //
  // All hit/cluster inputs are the coincidence-gated collections
  // (Trk/Time/CalCoincidence output), not the raw frame-wide *TimeAlign*
  // mirrors: the prefilter should see as little as possible, after
  // coincidence cleaning, not the whole frame.
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
          .m_default_output_tags = {"EventCandidatesFiltered", "PrefilterScores",
                                     "PrefilterLambdaK", "PrefilterSegLogits"},
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // Unfolder: materialises each surviving candidate as a PhysicsEvent.
  app->Add(new EventUnfolder());

  // Detector plugins (BTOF, MPGD, calorimeters, and so on) are not loaded
  // here. eventbuilder depends on their digitization output, including a
  // Timeslice-level, "Frame"-suffixed mirror of each detector's chain (see
  // e.g. src/detectors/BTOF/BTOF.cc) that TimeAlignment_factory reads. Each
  // detector's own InitPlugin already registers that mirror, alongside its
  // normal PhysicsEvent-level output. The standard eicrecon executable
  // already loads every default plugin, so calling a detector's InitPlugin
  // again here would only initialize it twice.
}
} // extern "C"
