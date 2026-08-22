// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
//
// TrkCoincidence_factory
// -----------------------
// Candidate-t0 gate for FAST-detector (TOF+MPGD) time-aligned hits — the
// exact same mechanism as TimeCoincidence_factory (slow/Si+B0), extended to
// the fast hits themselves. Sits between EventBuilder_factory and
// EventPrefilter_factory. For each fast hit in the input *TimeAlignRecHits
// collection, it forwards the hit only if it falls within any
// EventCandidate's acceptance window:
//
//   |t_hit - t0_candidate| <= nsigma × sigma_hit
//
// Unlike the slow-hit gate, fast hits already participated in FORMING t0
// (EventBuilder_factory's coincidence scan is what computed it) — gating
// them again by the same criterion is not circular, it discards the rest of
// the frame's fast-hit population that never contributed to any candidate,
// which is most of it. This exists so downstream consumers (currently
// EventPrefilter_factory) see a minimized, coincidence-cleaned hit set
// instead of the whole frame — "as little as possible, after cleaning
// coincidence," not the frame-wide *TimeAlignRecHits mirror.
//
// Output: *RecHits (fast) — PODIO subset collections of *TimeAlignRecHits.
// Not currently wired into EventUnfolder (which still reads raw
// *TimeAlignRecHits for fast dets and applies its own trk_nsigma_window+
// trk_late_nsigma_window_asym gate directly) — only EventPrefilter_factory
// consumes this today. late_nsigma_window_asym is still offered here,
// mirroring TimeCoincidence_factory/EventUnfolder, so behavior stays
// consistent if that wiring ever changes. It is a SIGMA MULTIPLE, not a flat
// ns constant — see TimeCoincidence_factory.h's header comment for why a
// flat ns add-on isn't comparable across detectors of very different
// resolution.
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

struct TrkCoincidence_factory : public JOmniFactory<TrkCoincidence_factory> {

  // EventCandidates supplies the t0 values for the coincidence gate.
  PodioInput<edm4hep::EventHeader> m_candidates_in{this, "EventCandidates"};

  // Fast-detector *TimeAlignRecHits (current frame only).
  VariadicPodioInput<edm4eic::TrackerHit, true> m_fast_hits_in{this, {}};

  // One subset output per fast detector.
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
      const auto* in = m_fast_hits_in().at(s);
      if (!in) continue;
      for (const auto& hit : *in) {
        const float sigma = hit.getTimeError();
        if (sigma > 0.0f && passes(hit.getTime(), sigma))
          out->push_back(hit);
      }
    }
  }
};
