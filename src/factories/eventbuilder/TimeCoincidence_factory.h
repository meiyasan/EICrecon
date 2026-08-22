// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// TimeCoincidence_factory
// -----------------------
// Candidate-t0 gate for slow-detector (Si/B0) time-aligned hits.
//
// Sits between EventBuilder_factory and EventUnfolder. For each slow hit in
// the input *TimeAlignRecHits collection, it forwards the hit only if it falls
// within any EventCandidate's acceptance window:
//
//   |t_hit - t0_candidate| <= nsigma × kSiSigma_ns
//
// kSiSigma_ns (2000 ns) is the MAPS integration window — a detector constant
// equal to the frame width. It is not configurable because it is fixed by the
// detector design and must not drift out of sync with EventBuilder_factory.
//
// Cross-frame boundary recovery
// -----------------------------
// Handled DOWNSTREAM by the EventUnfolder (eventbuilder.cc), which merges the
// adjacent frames' slow-detector hits from the TimesliceBuffer_service. This
// factory is stateless and gates the current frame only — factories must not
// carry state across frames (JANA2 runs one instance per in-flight timeslice).
// Control: -Peventbuilder:backward_frames / -Peventbuilder:forward_frames.
//
// Output: *TimeCoincRecHits — PODIO subset collections of *TimeAlignRecHits.
// All hits in the output originate from the current frame's *TimeAlignRecHits
// collection, so subset references are valid.
//
// Input wiring convention (see eventbuilder.cc)
//   index 0:   EventCandidates (edm4hep::EventHeader from EventBuilder_factory)
//   index 1-4: slow *TimeAlignRecHits (SiBarrelVertex, SiBarrelTracker,
//              SiEndcapTracker, B0Tracker)
// Output wiring
//   4 *TimeCoincRecHits subset collections → EventUnfolder (det. indices 6-9)
//
// Parameter
//   nsigma  N: gate half-width = N × kSiSigma_ns

#pragma once

#include <cmath>
#include <vector>

#include <extensions/jana/JOmniFactory.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4eic/TrackerHitCollection.h>

struct TimeCoincidence_factory : public JOmniFactory<TimeCoincidence_factory> {

  // MAPS integration window = frame width = detector constant (ns).
  static constexpr float kSiSigma_ns = 2000.0f;

  // EventCandidates supplies the t0 values for the coincidence gate.
  PodioInput<edm4hep::EventHeader> m_candidates_in{this, "EventCandidates"};

  // Slow-detector *TimeAlignRecHits (current frame only).
  VariadicPodioInput<edm4eic::TrackerHit, true> m_slow_hits_in{this, {}};

  // One subset output per slow detector.
  VariadicPodioOutput<edm4eic::TrackerHit> m_hits_out{this, {}};

  Parameter<float> nsigma{
      this, "nsigma", 3.0f,
      "N-sigma half-width: gate = N × 2000 ns (MAPS integration window)."};

  void Configure() {}
  void ChangeRun(int32_t) {}

  void Process(int64_t, uint64_t) {
    const size_t n_det    = m_hits_out().size();
    const float  half_win = nsigma() * kSiSigma_ns;

    const auto* cands = m_candidates_in();
    std::vector<float> t0s;
    if (cands)
      for (const auto& h : *cands)
        t0s.push_back(static_cast<float>(h.getWeights(0)));

    auto passes = [&](float t) -> bool {
      for (const float t0 : t0s)
        if (std::fabs(t - t0) <= half_win) return true;
      return false;
    };

    for (size_t s = 0; s < n_det; ++s) {
      auto& out = m_hits_out().at(s);
      // Subset collection: every hit originates from the current frame's
      // *TimeAlignRecHits collection produced by TrkTimeAlignment_factory.
      out->setSubsetCollection(true);
      const auto* in = m_slow_hits_in().at(s);
      if (!in) continue;
      for (const auto& hit : *in)
        if (passes(hit.getTime()))
          out->push_back(hit);
    }
  }
};
