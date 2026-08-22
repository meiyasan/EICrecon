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
//   |t_hit - t0_candidate| <= nsigma × sigma_hit
//
// sigma_hit is edm4eic::TrackerHit::getTimeError(), populated per-hit by each
// detector's own digitization config (10 ns for Si Barrel/Endcap, 8 ns for
// B0), exactly like EventBuilder_factory already does for TOF/MPGD. This
// replaces a previous hardcoded 2000 ns "MAPS integration window" constant
// that was ~200-600x wider than the detectors' actual configured resolution
// (3 sigma of 2000 ns = 6 us against a 2 us frame width — far wider than the
// frame itself, so the gate was not discriminating anything).
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
// This output feeds EventUnfolder's own Si/B0 gate directly (eventbuilder.cc,
// m_trk_in indices 6-9): a hit dropped here can never be recovered downstream,
// so the gate must be at least as permissive as EventUnfolder's own. In
// particular, EventUnfolder widens the LATE side of its per-hit gate by
// late_nsigma_window_asym × sigma_hit to absorb the beta<1 arrival-delay
// tail (massive/curling particles arrive strictly late; see eventbuilder.cc's
// trk_late_nsigma_window_asym comment) — a plain symmetric cut here would
// silently discard those hits before EventUnfolder ever gets to keep them,
// defeating that widening for exactly the detectors it exists to help.
// late_nsigma_window_asym therefore mirrors EventUnfolder's parameter and
// default, widening ONLY the + side, same as there. It is a SIGMA MULTIPLE,
// not a flat ns constant: a flat ns add-on isn't comparable across detectors
// of very different resolution (2 ns is ~0.2 sigma of a 10 ns Si/B0 hit but
// tens of sigma of a TOF hit), so the extra LATE allowance scales with each
// hit's own getTimeError() instead.
//
// Input wiring convention (see eventbuilder.cc)
//   index 0:   EventCandidates (edm4hep::EventHeader from EventBuilder_factory)
//   index 1-4: slow *TimeAlignRecHits (SiBarrelVertex, SiBarrelTracker,
//              SiEndcapTracker, B0Tracker)
// Output wiring
//   4 *TimeCoincRecHits subset collections → EventUnfolder (det. indices 6-9)
//
// Parameters
//   nsigma                   N: gate half-width = N × sigma_hit (per-hit
//                             getTimeError())
//   late_nsigma_window_asym  extra LATE acceptance, in units of sigma_hit,
//                             added to the + side only (0 = legacy symmetric
//                             gate)

#pragma once

#include <cmath>
#include <vector>

#include <extensions/jana/JOmniFactory.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4eic/TrackerHitCollection.h>

struct TimeCoincidence_factory : public JOmniFactory<TimeCoincidence_factory> {

  // EventCandidates supplies the t0 values for the coincidence gate.
  PodioInput<edm4hep::EventHeader> m_candidates_in{this, "EventCandidates"};

  // Slow-detector *TimeAlignRecHits (current frame only).
  VariadicPodioInput<edm4eic::TrackerHit, true> m_slow_hits_in{this, {}};

  // One subset output per slow detector.
  VariadicPodioOutput<edm4eic::TrackerHit> m_hits_out{this, {}};

  Parameter<float> nsigma{
      this, "nsigma", 3.0f,
      "N-sigma half-width: gate = N × sigma_hit (per-hit getTimeError())."};
  Parameter<float> late_nsigma_window_asym{
      this, "late_nsigma_window_asym", 0.5f,
      "extra LATE acceptance, in units of that hit's own sigma "
      "(getTimeError()), added to the + side only (beta<1 arrival delay; "
      "mirrors EventUnfolder's eventbuilder:unfolder:trk_late_nsigma_window_asym). 0 = "
      "legacy symmetric gate. UNMEASURED PLACEHOLDER — retune via an "
      "efficiency/purity scan."};

  void Configure() {}
  void ChangeRun(int32_t) {}

  void Process(int64_t, uint64_t) {
    const size_t n_det = m_hits_out().size();
    const float  n     = nsigma();
    const float  late  = late_nsigma_window_asym();

    const auto* cands = m_candidates_in();
    std::vector<float> t0s;
    if (cands)
      for (const auto& h : *cands)
        t0s.push_back(static_cast<float>(h.getWeights(0)));

    auto passes = [&](float t, float sigma) -> bool {
      const float half_win = n * sigma;
      for (const float t0 : t0s) {
        const float dt = t - t0;
        if (dt >= -half_win && dt <= half_win + late * sigma) return true;
      }
      return false;
    };

    for (size_t s = 0; s < n_det; ++s) {
      auto& out = m_hits_out().at(s);
      // Subset collection: every hit originates from the current frame's
      // *TimeAlignRecHits collection produced by TrkTimeAlignment_factory.
      out->setSubsetCollection(true);
      const auto* in = m_slow_hits_in().at(s);
      if (!in) continue;
      for (const auto& hit : *in) {
        const float sigma = hit.getTimeError();
        if (sigma > 0.0f && passes(hit.getTime(), sigma))
          out->push_back(hit);
      }
    }
  }
};
