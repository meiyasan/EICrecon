// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// Time-alignment factory, shared by trackers and calorimeter clusters.
//
// For each input collection it clones every hit/cluster, subtracts a straight-
// line propagation term r/c (so all times refer to a common collision t0), and
// sorts the collection by time. The logic is identical for tracker hits
// (edm4eic::TrackerHit) and calo clusters (edm4eic::Cluster), so it is written
// once as a template and instantiated once per family.
//
// This factory is STATELESS by design. JANA2 creates one factory instance per
// in-flight timeslice (jana:max_inflight_timeslices, default nthreads) and
// assigns frames to instances arbitrarily, so factory members must never carry
// information across frames: with nthreads > 1 such state fragments across
// instances and sees frames out of order.
//
// Cross-frame boundary recovery (slow detectors only)
// ---------------------------------------------------
// Si/MAPS detectors have σ_Si ≈ 2 µs ≈ frame width, so hits near a frame
// boundary can be reconstructed into a neighbouring frame. For slow detectors
// (indices >= NFast) this factory deposits each frame's corrected hits into
// the process-wide TimesliceBuffer_service, keyed by absolute frame number.
// The EventUnfolder (see eventbuilder.cc) fetches the adjacent frames from
// that service while gating hits per candidate — window controlled by
// -Peventbuilder:backward_frames / -Peventbuilder:forward_frames there.
// The *TimeAlignRecHits output itself always contains the current frame only.
//
// Aliases (defined at the bottom of this file)
//   TrkTimeAlignment_factory = TimeAlignment_factory<edm4eic::TrackerHit, 6>
//     NFast=6: indices 0-5 are fast (TOF, MPGD); indices 6-9 are slow (Si/B0)
//     and are deposited into the TimesliceBuffer_service.
//   CalTimeAlignment_factory = TimeAlignment_factory<edm4eic::Cluster>
//     NFast defaults to UINT32_MAX: all clusters are "fast", nothing deposited.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <extensions/jana/JOmniFactory.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/ClusterCollection.h>

#include "services/eventbuilder/TimesliceBuffer_service.h"

template <typename HitT,
          uint32_t NFast = std::numeric_limits<uint32_t>::max()>
struct TimeAlignment_factory
    : public JOmniFactory<TimeAlignment_factory<HitT, NFast>> {

  using Base    = JOmniFactory<TimeAlignment_factory<HitT, NFast>>;
  using MutHitT = decltype(std::declval<const HitT &>().clone());

  // Only the tracker family deposits into the cross-frame buffer.
  static constexpr bool kDeposits = std::is_same_v<HitT, edm4eic::TrackerHit>;

  // Collection names supplied by the wiring (see eventbuilder.cc).
  typename Base::template VariadicPodioInput<HitT, true> m_in{this, {}};
  typename Base::template VariadicPodioOutput<HitT>      m_out{this, {}};

  std::shared_ptr<TimesliceBuffer_service> m_ts_buffer;

  void Configure() {
    if constexpr (kDeposits)
      m_ts_buffer = this->GetApplication()->template GetService<TimesliceBuffer_service>();
  }

  void ChangeRun(int32_t) {}

  // Apply the r/c propagation correction to a hit and return a mutable clone.
  MutHitT correct(const HitT& hit) const {
    MutHitT c = hit.clone();
    const double x = hit.getPosition()[0];
    const double y = hit.getPosition()[1];
    const double z = hit.getPosition()[2];
    c.setTime(hit.getTime() - std::sqrt(x*x + y*y + z*z) * 0.0034);
    return c;
  }

  void Process(int64_t, uint64_t frame_nr) {
    const size_t n_out  = m_out().size();
    const size_t nf     = static_cast<size_t>(NFast);
    const size_t n_slow = (nf < n_out) ? (n_out - nf) : 0u;

    TimesliceBuffer_service::FrameHits slow_frame;
    if constexpr (kDeposits)
      slow_frame.resize(n_slow);

    for (size_t i = 0; i < n_out; ++i) {
      const auto* in  = m_in().at(i);
      auto&       out = m_out().at(i);

      std::vector<MutHitT> current;
      if (in != nullptr) {
        current.reserve(in->size());
        for (const auto& hit : *in)
          current.push_back(correct(hit));
      }
      std::sort(current.begin(), current.end(),
                [](const auto& a, const auto& b) { return a.getTime() < b.getTime(); });

      if constexpr (kDeposits) {
        if (i >= nf) {
          // Slow detector: the output gets clones; the originals stay
          // untracked and go to the cross-frame deposit below.
          for (const auto& h : current)
            out->push_back(h.clone());
          slow_frame[i - nf] = std::move(current);
          continue;
        }
      }
      for (auto& h : current)
        out->push_back(h);
    }

    if constexpr (kDeposits) {
      if (n_slow > 0)
        m_ts_buffer->deposit(frame_nr, std::move(slow_frame));
    }
  }
};

// TrkTimeAlignment_factory: NFast=6 means indices 0-5 (TOF, MPGD) are fast and
// indices 6-9 (Si/B0) are slow — their corrected hits are deposited into the
// TimesliceBuffer_service for cross-frame recovery by the EventUnfolder.
using TrkTimeAlignment_factory = TimeAlignment_factory<edm4eic::TrackerHit, 6>;

// CalTimeAlignment_factory: all calo clusters treated as fast (no deposit).
using CalTimeAlignment_factory = TimeAlignment_factory<edm4eic::Cluster>;
