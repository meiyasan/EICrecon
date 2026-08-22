// Copyright 2024, Jefferson Science Associates, LLC.
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// EventUnfolder
// -------------
// Minimal JANA unfolder. It performs NO analysis: it reads the candidate-event
// list produced by EventBuilderFactory ("EventCandidates", an
// edm4hep::EventHeaderCollection carrying t0/dt0) and materialises one
// PhysicsEvent per candidate.
//
// For candidate child_idx it:
//   * selects tracker RecHits per detector inside the resolution-matched window
//     [t0 - N*sigma_det, t0 + N*sigma_det] (TOF tight, Silicon wide),
//   * carries through raw hits / links / hit associations,
//   * carries through calo clusters and their (matched) particle associations,
//   * clones the MC particles,
//   * writes a child EventHeader holding the real t0 and dt0.
//
// All the event-finding / coincidence / topology logic lives in
// EventBuilderFactory; this class is just structural glue (only a JEventUnfolder
// can create finer-level child events in JANA).

#pragma once

#include <cmath>
#include <cstdint>
#include <unordered_set>

#include <JANA/JEventUnfolder.h>

#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociation.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>

struct EventUnfolder : public JEventUnfolder {

  // Acceptance-window controls. Same meaning/defaults as EventBuilderFactory:
  // the per-detector half-width is nsigma_window * timeResolution_<group>.
  Parameter<float> timeResolution_TOF{this, "timeResolution_TOF", 0.03,
                                      "time resolution of TOF detectors in ns"};
  Parameter<float> timeResolution_MPGD{this, "timeResolution_MPGD", 10.0,
                                       "time resolution of MPGD detectors in ns"};
  Parameter<float> timeResolution_Silicon{this, "timeResolution_Silicon", 2000.0,
                                          "time resolution of Silicon detectors in ns"};
  Parameter<float> nsigma_window{this, "nsigma_window", 3.0,
                                 "N-sigma half-width for per-detector hit acceptance"};

  size_t m_event_number_ts = 0;

  // -- tracker RecHits (10 detectors) ------------------------------------------
  std::vector<std::string> m_trk_in_names = {
      "TOFBarrelTimeAlignRecHits",          "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",         "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits", "ForwardMPGDEndcapTimeAlignRecHits",
      "SiBarrelVertexTimeAlignRecHits",     "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits",    "B0TrackerTimeAlignRecHits"};
  std::vector<std::string> m_trk_out_names = {
      "TOFBarrelRecHits",       "TOFEndcapRecHits",          "MPGDBarrelRecHits",
      "OuterMPGDBarrelRecHits", "BackwardMPGDEndcapRecHits", "ForwardMPGDEndcapRecHits",
      "SiBarrelVertexRecHits",  "SiBarrelTrackerRecHits",    "SiEndcapTrackerRecHits",
      "B0TrackerRecHits"};

  // -- raw hits / links / hit associations (10 detectors each) -----------------
  std::vector<std::string> m_raw_in_names = {
      "TOFBarrelRawHitDigi",       "TOFEndcapRawHitDigi",          "MPGDBarrelRawHitDigi",
      "OuterMPGDBarrelRawHitDigi", "BackwardMPGDEndcapRawHitDigi", "ForwardMPGDEndcapRawHitDigi",
      "SiBarrelVertexRawHitDigi",  "SiBarrelRawHitDigi",           "SiEndcapTrackerRawHitDigi",
      "B0TrackerRawHitDigi"};
  std::vector<std::string> m_raw_out_names = {
      "TOFBarrelRawHits",       "TOFEndcapRawHits",          "MPGDBarrelRawHits",
      "OuterMPGDBarrelRawHits", "BackwardMPGDEndcapRawHits", "ForwardMPGDEndcapRawHits",
      "SiBarrelVertexRawHits",  "SiBarrelRawHits",           "SiEndcapTrackerRawHits",
      "B0TrackerRawHits"};

  std::vector<std::string> m_link_in_names = {
      "TOFBarrelRawHitLinkDigi",       "TOFEndcapRawHitLinkDigi",
      "MPGDBarrelRawHitLinkDigi",      "OuterMPGDBarrelRawHitLinkDigi",
      "BackwardMPGDEndcapRawHitLinkDigi", "ForwardMPGDEndcapRawHitLinkDigi",
      "SiBarrelVertexRawHitLinkDigi",  "SiBarrelRawHitLinkDigi",
      "SiEndcapTrackerRawHitLinkDigi", "B0TrackerRawHitLinkDigi"};
  std::vector<std::string> m_link_out_names = {
      "TOFBarrelRawHitLinks",       "TOFEndcapRawHitLinks",          "MPGDBarrelRawHitLinks",
      "OuterMPGDBarrelRawHitLinks", "BackwardMPGDEndcapRawHitLinks", "ForwardMPGDEndcapRawHitLinks",
      "SiBarrelVertexRawHitLinks",  "SiBarrelRawHitLinks",           "SiEndcapTrackerRawHitLinks",
      "B0TrackerRawHitLinks"};

  std::vector<std::string> m_asso_in_names = {
      "TOFBarrelRawHitAssociationDigi",       "TOFEndcapRawHitAssociationDigi",
      "MPGDBarrelRawHitAssociationDigi",      "OuterMPGDBarrelRawHitAssociationDigi",
      "BackwardMPGDEndcapRawHitAssociationDigi", "ForwardMPGDEndcapRawHitAssociationDigi",
      "SiBarrelVertexRawHitAssociationDigi",  "SiBarrelRawHitAssociationDigi",
      "SiEndcapTrackerRawHitAssociationDigi", "B0TrackerRawHitAssociationDigi"};
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

  // -- candidate list from EventBuilderFactory ---------------------------------
  // All inputs live in the parent time-frame (Timeslice level); outputs are
  // written into the child PhysicsEvent. Hence .level = Timeslice on every input
  // so JANA fetches them from the parent rather than the (empty) child.
  PodioInput<edm4hep::EventHeader> m_candidates_in{
      this, {.name = "EventCandidates", .level = JEventLevel::Timeslice}};

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

  // TODO(calo): re-enable once the calo digi/clustering registrations are all
  // moved to the Timeslice level. Some calo cluster->particle association
  // factories are still registered at PhysicsEvent (a pre-existing, half-finished
  // migration), which makes fetching them from the time-frame parent throw a
  // level-mismatch. The EventBuilder itself uses only trackers, so calo
  // carry-through is disabled here to keep the tracker t0/dt0 path working.
  // VariadicPodioInput<edm4eic::Cluster> m_clu_in{
  //     this, {.names = m_clu_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  // VariadicPodioOutput<edm4eic::Cluster> m_clu_out{this, m_clu_out_names};
  //
  // VariadicPodioInput<edm4eic::MCRecoClusterParticleAssociation> m_cluasso_in{
  //     this, {.names = m_cluasso_in_names, .level = JEventLevel::Timeslice, .is_optional = true}};
  // VariadicPodioOutput<edm4eic::MCRecoClusterParticleAssociation> m_cluasso_out{
  //     this, m_cluasso_out_names};

  EventUnfolder() {
    SetTypeName(NAME_OF_THIS);
    SetParentLevel(JEventLevel::Timeslice);
    SetChildLevel(JEventLevel::PhysicsEvent);
  }

  // (not const: JANA Parameter::operator() is non-const)
  float detTimeRes(size_t detIdx) {
    if (detIdx < 2)
      return timeResolution_TOF();
    if (detIdx < 6)
      return timeResolution_MPGD();
    return timeResolution_Silicon();
  }

  Result Unfold(const JEvent& parent, JEvent& child, int child_idx) override {
    const auto& cands = *m_candidates_in();
    const int n_cand  = static_cast<int>(cands.size());
    if (n_cand == 0 || child_idx >= n_cand)
      return Result::KeepChildNextParent; // nothing (more) to emit in this frame

    auto cand          = cands.at(child_idx);
    const double t0    = cand.getWeights(0); // ns
    const double dt0   = cand.getWeights(1); // ns
    const double phys  = cand.getWeights(2);

    // -- tracker RecHits: per-detector resolution-matched window --------------
    for (size_t det = 0; det < m_trk_in().size(); ++det) {
      auto& out = m_trk_out().at(det);
      out->setSubsetCollection(true);
      const auto* in = m_trk_in().at(det);
      if (in == nullptr)
        continue;
      const double win = double(nsigma_window()) * double(detTimeRes(det));
      for (const auto& hit : *in)
        if (std::abs(hit.getTime() - t0) <= win)
          out->push_back(hit);
    }

    // -- raw hits / links / hit associations: carry-through (subset) ----------
    copyAll(m_raw_in(), m_raw_out());
    copyAll(m_link_in(), m_link_out());
    copyAll(m_asso_in(), m_asso_out());

    // -- calo clusters + matched cluster->particle associations ---------------
    // TODO(calo): temporarily disabled pending the calo digi/clustering
    // Timeslice-level migration (see member declarations above). The matching
    // logic is preserved here for restoration.
    // for (size_t cd = 0; cd < m_clu_in().size(); ++cd) {
    //   auto& cout = m_clu_out().at(cd);
    //   cout->setSubsetCollection(true);
    //   auto& aout = m_cluasso_out().at(cd);
    //   aout->setSubsetCollection(true);
    //   const auto* cin = m_clu_in().at(cd);
    //   if (cin == nullptr) continue;
    //   std::unordered_set<uint32_t> included;
    //   for (const auto& clu : *cin) {
    //     cout->push_back(clu);
    //     included.insert(static_cast<uint32_t>(clu.getObjectID().index));
    //   }
    //   const auto* ain = (cd < m_cluasso_in().size()) ? m_cluasso_in().at(cd) : nullptr;
    //   if (ain == nullptr) continue;
    //   for (const auto& as : *ain)
    //     if (included.count(static_cast<uint32_t>(as.getRec().getObjectID().index)))
    //       aout->push_back(as);
    // }

    // -- MC particles: clone all (truth bookkeeping) --------------------------
    for (const auto& mcp : *m_mcparticles_in())
      m_mcparticles_out()->push_back(mcp.clone(true));

    // -- child EventHeader with the real t0/dt0 -------------------------------
    child.SetEventNumber(m_event_number_ts);
    child.SetRunNumber(parent.GetRunNumber());

    edm4hep::MutableEventHeader hdr;
    hdr.setRunNumber(parent.GetRunNumber());
    hdr.setEventNumber(m_event_number_ts);
    hdr.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    hdr.setWeight(dt0);
    hdr.addToWeights(t0);
    hdr.addToWeights(dt0);
    hdr.addToWeights(phys);
    m_event_header_ts_out()->push_back(hdr);
    ++m_event_number_ts;

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
  template <typename InVec, typename OutVec> void copyAll(const InVec& in, OutVec& out) {
    for (size_t i = 0; i < out.size(); ++i) {
      out.at(i)->setSubsetCollection(true);
      if (i >= in.size() || in.at(i) == nullptr)
        continue;
      for (const auto& elem : *in.at(i))
        out.at(i)->push_back(elem);
    }
  }
};
