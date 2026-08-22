// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>

#include <JANA/JApplication.h>
#include <JANA/JEventUnfolder.h>
#include <JANA/JException.h>

#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociation.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/eventbuilder/TimeAlignment_factory.h"
#include "factories/eventbuilder/TimeCoincidence_factory.h"
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

  // Acceptance-window controls. Fast dets (TOF/MPGD): half-width is computed
  // per-hit from edm4eic::TrackerHit::getTimeError(), which each detector's own
  // TrackerHitReconstruction_factory already populates from its digitization
  // config — no separate copy kept here. Slow dets (Si/B0): timeResolution_Silicon
  // is the MAPS shutter/integration-window half-width, a genuine detector
  // constant (not a per-hit measurement, so there's nothing on the hit to read);
  // mirrors TimeCoincidence_factory's kSiSigma_ns.
  Parameter<float> timeResolution_Silicon{this, "timeResolution_Silicon", 2000.0,
                                          "time resolution of Silicon detectors in ns"};
  Parameter<float> nsigma_window{this, "nsigma_window", 3.0,
                                 "N-sigma half-width for per-detector hit acceptance"};

  // Prefilter ONNX parameters — registered at eventbuilder: prefix so the user
  // sets -Peventbuilder:onnx_model=... (read by EventPrefilter_factory::Configure).
  Parameter<std::string> onnx_model{this, "onnx_model", "",
                                    "Path to frame-level ONNX prefilter model. Empty = pass-through."};
  Parameter<float>   onnx_score_threshold{this, "score_threshold", 0.5f,
                                          "Maximum BKG probability to pass (GNN class-0 score)."};
  Parameter<int32_t> onnx_max_hits{this, "max_hits", 0,
                                   "Fixed input length for padded models (e.g. 2500). 0 = use actual hit count."};
  Parameter<int32_t> onnx_k_neighbors{this, "k_neighbors", 16,
                                      "k for kNN graph construction (per DGCNN layer)."};
  Parameter<float>   onnx_time_weight{this, "time_weight", 1.0f,
                                      "Scale factor on the time dimension in layer-0 spacetime kNN."};

  // Cross-frame Si/B0 recovery window. Slow-detector hits of the adjacent
  // frames are fetched from the TimesliceBuffer_service (deposited there by
  // TrkTimeAlignment_factory) and gated per candidate like current-frame hits.
  // forward_frames adds no output latency: it only waits for the deposit of
  // frames that are already in flight (hence the max_inflight requirement).
  Parameter<int32_t> backward_frames{this, "backward_frames", 1,
      "Past frames whose Si/B0 hits are included when gating each candidate. 0 = off."};
  Parameter<int32_t> forward_frames{this, "forward_frames", 0,
      "Future frames whose Si/B0 hits are included when gating each candidate. 0 = off. "
      "Requires nthreads >= 2 and jana:max_inflight_timeslices > forward_frames."};

  std::shared_ptr<TimesliceBuffer_service> m_ts_buffer;
  bool m_checked_config = false;

  // Adjacent-frame slow-det hits, fetched once per parent (Unfold is sequential).
  uint64_t m_adjacent_parent = std::numeric_limits<uint64_t>::max();
  std::vector<TimesliceBuffer_service::FrameHits> m_adjacent;

  // -- tracker RecHits (10 detectors) ------------------------------------------
  // Fast detectors (TOF, MPGD): read directly from *TimeAlignRecHits.
  // Slow detectors (Si, B0): read from TimeCoincidence_factory output
  // (*TimeCoincRecHits, current frame); adjacent-frame hits are merged in
  // Unfold() from the TimesliceBuffer_service (backward_frames/forward_frames).
  std::vector<std::string> m_trk_in_names = {
      "TOFBarrelTimeAlignRecHits",          "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",         "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits", "ForwardMPGDEndcapTimeAlignRecHits",
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
  std::vector<std::string> m_raw_in_names = {
      "TOFBarrelRawHitDigi",          "TOFEndcapRawHitDigi",
      "MPGDBarrelRawHitDigi",         "OuterMPGDBarrelRawHitDigi",
      "BackwardMPGDEndcapRawHitDigi", "ForwardMPGDEndcapRawHitDigi",
      "SiBarrelVertexRawHitDigi",     "SiBarrelRawHitDigi",
      "SiEndcapTrackerRawHitDigi",    "B0TrackerRawHitDigi"};
  std::vector<std::string> m_raw_out_names = {
      "TOFBarrelRawHits",       "TOFEndcapRawHits",          "MPGDBarrelRawHits",
      "OuterMPGDBarrelRawHits", "BackwardMPGDEndcapRawHits", "ForwardMPGDEndcapRawHits",
      "SiBarrelVertexRawHits",  "SiBarrelRawHits",           "SiEndcapTrackerRawHits",
      "B0TrackerRawHits"};

  std::vector<std::string> m_link_in_names = {
      "TOFBarrelRawHitLinkDigi",          "TOFEndcapRawHitLinkDigi",
      "MPGDBarrelRawHitLinkDigi",         "OuterMPGDBarrelRawHitLinkDigi",
      "BackwardMPGDEndcapRawHitLinkDigi", "ForwardMPGDEndcapRawHitLinkDigi",
      "SiBarrelVertexRawHitLinkDigi",     "SiBarrelRawHitLinkDigi",
      "SiEndcapTrackerRawHitLinkDigi",    "B0TrackerRawHitLinkDigi"};
  std::vector<std::string> m_link_out_names = {
      "TOFBarrelRawHitLinks",       "TOFEndcapRawHitLinks",          "MPGDBarrelRawHitLinks",
      "OuterMPGDBarrelRawHitLinks", "BackwardMPGDEndcapRawHitLinks", "ForwardMPGDEndcapRawHitLinks",
      "SiBarrelVertexRawHitLinks",  "SiBarrelRawHitLinks",           "SiEndcapTrackerRawHitLinks",
      "B0TrackerRawHitLinks"};

  std::vector<std::string> m_asso_in_names = {
      "TOFBarrelRawHitAssociationDigi",          "TOFEndcapRawHitAssociationDigi",
      "MPGDBarrelRawHitAssociationDigi",         "OuterMPGDBarrelRawHitAssociationDigi",
      "BackwardMPGDEndcapRawHitAssociationDigi", "ForwardMPGDEndcapRawHitAssociationDigi",
      "SiBarrelVertexRawHitAssociationDigi",     "SiBarrelRawHitAssociationDigi",
      "SiEndcapTrackerRawHitAssociationDigi",    "B0TrackerRawHitAssociationDigi"};
  std::vector<std::string> m_asso_out_names = {
      "TOFBarrelRawHitAssociations",       "TOFEndcapRawHitAssociations",
      "MPGDBarrelRawHitAssociations",      "OuterMPGDBarrelRawHitAssociations",
      "BackwardMPGDEndcapRawHitAssociations", "ForwardMPGDEndcapRawHitAssociations",
      "SiBarrelVertexRawHitAssociations",  "SiBarrelRawHitAssociations",
      "SiEndcapTrackerRawHitAssociations", "B0TrackerRawHitAssociations"};

  // -- calo clusters + cluster->particle associations (4 detectors) ------------
  std::vector<std::string> m_clu_in_names = {
      "B0ECalTimeAlignClusters", "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters", "EcalEndcapPTimeAlignClusters"};
  std::vector<std::string> m_clu_out_names = {
      "B0ECalClusters", "EcalBarrelClusters", "EcalEndcapNClusters", "EcalEndcapPClusters"};
  std::vector<std::string> m_cluasso_in_names = {
      "B0ECalClusterAssociationDigi", "EcalBarrelClusterAssociationDigi",
      "EcalEndcapNClusterAssociationDigi", "EcalEndcapPClusterAssociationDigi"};
  std::vector<std::string> m_cluasso_out_names = {
      "B0ECalClusterAssociations", "EcalBarrelClusterAssociations",
      "EcalEndcapNClusterAssociations", "EcalEndcapPClusterAssociations"};

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

  PodioInput<edm4hep::EventHeader> m_event_header_in{
      this, {.name = "EventHeader", .level = JEventLevel::Timeslice, .is_optional = true}};
  PodioOutput<edm4hep::EventHeader> m_event_header_out{this, "EventHeader"};
  PodioOutput<edm4hep::EventHeader> m_event_header_ts_out{this, "EventHeader_TS"};

  PodioInput<edm4hep::MCParticle> m_mcparticles_in{
      this, {.name = "MCParticles", .level = JEventLevel::Timeslice}};
  PodioOutput<edm4hep::MCParticle> m_mcparticles_out{this, "MCParticles"};

  VariadicPodioInput<edm4eic::TrackerHit> m_trk_in{
      this, {.names = m_trk_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  VariadicPodioOutput<edm4eic::TrackerHit> m_trk_out{this, m_trk_out_names};

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

  // B0ECAL, BEMC, EEMC are all registered at Timeslice level; FEMC (EcalEndcapP)
  // has no plugin yet so its optional inputs will be nullptr (handled below).
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

  // Slow-detector (Si/B0) gate half-width -- see timeResolution_Silicon above.
  // (not const: Parameter::operator() is non-const)
  double slowDetWindow() { return double(nsigma_window()) * double(timeResolution_Silicon()); }

  // Fast (trigger) detectors are indices 0-5; slow (Si/B0) are 6-9.
  static constexpr size_t kNumFastDet = 6;

  // Validate the forward_frames configuration once. Waiting for a future
  // frame's deposit only terminates if that frame can be processed while this
  // one sits in the (sequential) unfold stage.
  void checkConfig() {
    if (m_checked_config)
      return;
    m_checked_config = true;
    const int32_t n_fwd = forward_frames();
    if (n_fwd <= 0)
      return;
    auto* app = GetApplication();
    const int nthreads = app->GetNThreads();
    const int inflight = app->GetParameterValue<int>("jana:max_inflight_timeslices");
    if (nthreads < 2 || inflight <= n_fwd)
      throw JException("eventbuilder:forward_frames=%d requires nthreads >= 2 and "
                       "jana:max_inflight_timeslices > forward_frames "
                       "(got nthreads=%d, max_inflight_timeslices=%d)",
                       n_fwd, nthreads, inflight);
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
      const int32_t n_fwd  = std::max(0, forward_frames());
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
    // references), windowed by the shared MAPS shutter half-width.
    const double slow_win = slowDetWindow();
    for (size_t det = 0; det < m_trk_in().size(); ++det) {
      auto& out = m_trk_out().at(det);
      const auto* in = m_trk_in().at(det);

      if (det < kNumFastDet) {
        out->setSubsetCollection(true);
        if (in == nullptr)
          continue;
        for (const auto& hit : *in) {
          const double win = double(nsigma_window()) * double(hit.getTimeError());
          if (std::abs(hit.getTime() - t0) <= win)
            out->push_back(hit);
        }
        continue;
      }

      std::vector<edm4eic::MutableTrackerHit> merged;
      if (in != nullptr)
        for (const auto& hit : *in)
          if (std::abs(hit.getTime() - t0) <= slow_win)
            merged.push_back(hit.clone());
      const size_t si = det - kNumFastDet;
      for (const auto& frame : m_adjacent)
        if (si < frame.size())
          for (const auto& hit : frame[si])
            if (std::abs(hit.getTime() - t0) <= slow_win)
              merged.push_back(hit.clone());
      std::sort(merged.begin(), merged.end(),
                [](const auto& a, const auto& b) { return a.getTime() < b.getTime(); });
      for (auto& hit : merged)
        out->push_back(hit);
    }

    // -- raw hits / links / hit associations: carry-through (subset) -----------
    copyAll(m_raw_in(), m_raw_out());
    copyAll(m_link_in(), m_link_out());
    copyAll(m_asso_in(), m_asso_out());

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

    // -- calo clusters + matched cluster->particle associations ----------------
    for (size_t cd = 0; cd < m_clu_in().size(); ++cd) {
      auto& cout = m_clu_out().at(cd);
      cout->setSubsetCollection(true);
      auto& aout = m_cluasso_out().at(cd);
      aout->setSubsetCollection(true);
      const auto* cin = m_clu_in().at(cd);
      if (cin == nullptr) continue; // optional: e.g. FEMC not yet registered
      std::unordered_set<uint32_t> included;
      for (const auto& clu : *cin) {
        cout->push_back(clu);
        included.insert(static_cast<uint32_t>(clu.getObjectID().index));
      }
      const auto* ain = (cd < m_cluasso_in().size()) ? m_cluasso_in().at(cd) : nullptr;
      if (ain == nullptr) continue;
      for (const auto& as : *ain)
        if (included.count(static_cast<uint32_t>(as.getRec().getObjectID().index)))
          aout->push_back(as);
    }

    // -- MC particles: clone all (truth bookkeeping) ---------------------------
    for (const auto& mcp : *m_mcparticles_in())
      m_mcparticles_out()->push_back(mcp.clone(true));

    // -- child EventHeader with the real t0/t0sigma --------------------------------
    // Deterministic numbering: with nthreads > 1 parents can reach the unfolder
    // out of order, so a running counter is not reproducible run-to-run. Encode
    // (frame, candidate) instead; candidates per 2 µs frame stay far below 1000.
    const uint64_t child_event_nr = parent_nr * 1000 + static_cast<uint64_t>(child_idx);
    child.SetEventNumber(child_event_nr);
    child.SetRunNumber(parent.GetRunNumber());

    edm4hep::MutableEventHeader hdr;
    hdr.setRunNumber(parent.GetRunNumber());
    hdr.setEventNumber(child_event_nr);
    hdr.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    hdr.setWeight(t0sigma);
    hdr.addToWeights(t0);
    hdr.addToWeights(t0sigma);
    hdr.addToWeights(phys);
    m_event_header_ts_out()->push_back(hdr);

    // Pass through the timeframe EventHeader as a reference, if present.
    if (m_event_header_in() != nullptr && m_event_header_in()->size() > 0) {
      m_event_header_out()->setSubsetCollection(true);
      m_event_header_out()->push_back(m_event_header_in()->at(0));
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
};


// ---------------------------------------------------------------------------
// Plugin wiring
// ---------------------------------------------------------------------------

void InitPlugin_BTOF(JApplication* app);
void InitPlugin_MPGD(JApplication* app);
void InitPlugin_BVTX(JApplication* app);
void InitPlugin_BTRK(JApplication* app);
void InitPlugin_ECTRK(JApplication* app);
void InitPlugin_ECTOF(JApplication* app);
void InitPlugin_B0TRK(JApplication* app);
// void InitPlugin_DIRC(JApplication* app);
// void InitPlugin_DRICH(JApplication* app);
void InitPlugin_FOFFMTRK(JApplication* app);
// void InitPlugin_PFRICH(JApplication* app);
// void InitPlugin_LOWQ2(JApplication* app);

void InitPlugin_B0ECAL(JApplication* app);
void InitPlugin_BEMC(JApplication* app);
void InitPlugin_EEMC(JApplication* app);
void InitPlugin_FEMC(JApplication* app);
void InitPlugin_EHCAL(JApplication* app);
// void InitPlugin_ECHAL(JApplication* app);
// void InitPlugin_BHCAL(JApplication* app);
// void InitPlugin_FHCAL(JApplication* app);
// void InitPlugin_LUMISPECCAL(JApplication* app);
// void InitPlugin_ZDC(JApplication* app);


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
      "TOFBarrelRecHitDigi",
      "TOFEndcapRecHitDigi",
      "MPGDBarrelRecHitDigi",
      "OuterMPGDBarrelRecHitDigi",
      "BackwardMPGDEndcapRecHitDigi",
      "ForwardMPGDEndcapRecHitDigi",
      "SiBarrelVertexRecHitDigi",
      "SiBarrelTrackerRecHitDigi",
      "SiEndcapTrackerRecHitDigi",
      "B0TrackerRecHitDigi"
  };
  // "TaggerTrackerRecHitDigi",
  // "DIRCBarRecHitDigi",
  // "DRICHRecHitDigi",
  // "ForwardOffMTrackerRecHitDigi",
  // "ForwardRomanPotRecHitDigi",
  // "LumiSpecRecHitDigi",
  // "RICHEndcapNRecHitDigi"

  std::vector<std::string> m_simcalocluster_collection_names_aligned = {
      "B0ECalTimeAlignClusters",
      "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters",
      "EcalEndcapPTimeAlignClusters"
  };
  // "EcalFarForwardZDCTimeAlignClusters",
  // "EcalLumiSpecTimeAlignClusters",
  // "HcalBarrelTimeAlignClusters",
  // "HcalEndcapNTimeAlignClusters",
  // "HcalEndcapPInsertTimeAlignClusters",
  // "HcalFarForwardZDCTimeAlignClusters",
  // "LFHCALTimeAlignClusters"

  std::vector<std::string> m_simcalocluster_collection_names = {
      "B0ECalClusterDigi",
      "EcalBarrelClusterDigi",
      "EcalEndcapNClusterDigi",
      "EcalEndcapPClusterDigi"
  };
  // "EcalFarForwardZDCClusterDigi",
  // "EcalLumiSpecClusterDigi",
  // "HcalBarrelClusterDigi",
  // "HcalEndcapNClusterDigi",
  // "HcalEndcapPInsertClusterDigi",
  // "HcalFarForwardZDCClusterDigi",
  // "LFHCALClusterDigi",
  // "EcalBarrelImagingClusterDigi",
  // "EcalBarrelScFiClusterDigi",
  // "EcalEndcapNImagingClusterDigi",
  // "EcalEndcapPImagingClusterDigi",
  // "EcalFarForwardZDCImagingClusterDigi",
  // "EcalLumiSpecImagingClusterDigi"

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
      "MCParticles"
  };
  app->Add(new JOmniFactoryGeneratorT<EventBuilder_factory>(
      JOmniFactoryGeneratorT<EventBuilder_factory>::TypedWiring{
          .m_tag                 = "eventBuilder",
          .m_default_input_tags  = m_eventbuilder_input_tags,
          .m_default_output_tags = {"EventCandidates"},
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // TimeCoincidence_factory: candidate-t0 gate for slow (Si/B0) RecHits.
  // Gates the current frame's *TimeAlignRecHits on each candidate's t0 window
  // (nsigma × 2000 ns) and produces a subset collection for EventUnfolder.
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
          .m_tag                 = "coinc",
          .m_default_input_tags  = coinc_inputs,
          .m_default_output_tags = coinc_outputs,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // EventPrefilter_factory: optional GNN background rejection between EventBuilder
  // and EventUnfolder. When eventbuilder:onnx_model is empty (the default), all
  // EventCandidates are forwarded unchanged and no ONNX session is created.
  // When a model path is set, only candidates with
  // probs[BKG=0] < eventbuilder:score_threshold survive.
  //
  // Hit features per candidate: [x/4000, y/4000, z/5000, t/50, eDep/10, det_id/6]
  // Input: first 6 aligned tracker collections (TOF + MPGD) for hit features.
  const std::vector<std::string> prefilter_inputs = {
      "EventCandidates",
      "TOFBarrelTimeAlignRecHits",          "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",         "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits", "ForwardMPGDEndcapTimeAlignRecHits"
  };
  app->Add(new JOmniFactoryGeneratorT<EventPrefilter_factory>(
      JOmniFactoryGeneratorT<EventPrefilter_factory>::TypedWiring{
          .m_tag                 = "prefilter",
          .m_default_input_tags  = prefilter_inputs,
          .m_default_output_tags = {"EventCandidatesFiltered"},
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // Unfolder: materialises each surviving candidate as a PhysicsEvent.
  app->Add(new EventUnfolder());

  InitPlugin_BTOF(app);
  InitPlugin_MPGD(app);
  InitPlugin_BVTX(app);
  InitPlugin_BTRK(app);
  InitPlugin_ECTRK(app);
  InitPlugin_ECTOF(app);
  InitPlugin_B0TRK(app);
  InitPlugin_FOFFMTRK(app);
  // InitPlugin_DIRC(app);
  // InitPlugin_DRICH(app);
  // InitPlugin_PFRICH(app);
  // InitPlugin_LOWQ2(app);

  InitPlugin_B0ECAL(app);
  InitPlugin_BEMC(app);
  InitPlugin_EEMC(app);
  InitPlugin_FEMC(app);
  InitPlugin_EHCAL(app);
  // InitPlugin_ECHAL(app);
  // InitPlugin_BHCAL(app);
  // InitPlugin_FHCAL(app);
  // InitPlugin_LUMISPECCAL(app);
  // InitPlugin_ZDC(app);

}
} // extern "C"
