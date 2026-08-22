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
// Aliases (defined at the bottom of this file). The index order is fixed by
// the wiring lists in eventbuilder.cc (m_simtrackerhit_collection_names /
// m_simcalocluster_collection_names) — keep the three in sync.
//
//   TrkTimeAlignment_factory = TimeAlignment_factory<edm4eic::TrackerHit, 6>
//     Aligns raw HIT times — it runs BEFORE any tracking, on hits from the
//     tracking detectors. 10 subsystems wired today, NFast=6:
//       fast (trigger group, sigma read per hit):
//         0  TOFBarrel              AC-LGAD barrel time-of-flight  (~30 ps)
//         1  TOFEndcap              AC-LGAD forward TOF            (~25 ps)
//         2  MPGDBarrel             inner MPGD barrel              (~10 ns)
//         3  OuterMPGDBarrel        outer MPGD barrel              (~10 ns)
//         4  BackwardMPGDEndcap     electron-side MPGD disks       (~10 ns)
//         5  ForwardMPGDEndcap      hadron-side MPGD disks         (~10 ns)
//       slow (deposited into TimesliceBuffer_service for cross-frame recovery):
//         6  SiBarrelVertex         MAPS vertex barrel             (10 ns sim)
//         7  SiBarrelTracker        MAPS sagitta barrel            (10 ns sim)
//         8  SiEndcapTracker        MAPS endcap disks              (10 ns sim)
//         9  B0Tracker              far-forward B0 spectrometer    (8 ns sim)
//     Future candidates (commented in the wiring): TaggerTracker (LOW-Q2),
//     DIRC bars, DRICH — appending keeps existing indices stable; a new SLOW
//     detector must go at the end AND NFast stays 6 unless it is a genuine
//     trigger detector.
//
//   CalTimeAlignment_factory = TimeAlignment_factory<edm4eic::Cluster>
//     Same r/c correction for calorimeter clusters. 9 collections wired today;
//     NFast defaults to UINT32_MAX: all clusters are "fast", nothing deposited.
//         0  B0ECal                 far-forward EM calorimeter
//         1  EcalBarrel             barrel EM calorimeter
//         2  EcalEndcapN            electron-side EM calorimeter
//         3  EcalEndcapP            hadron-side EM calorimeter
//         4  HcalBarrel             barrel hadronic calorimeter (plugins/BHCAL.cc)
//         5  HcalEndcapPInsert      hadron-side insert HCal   (plugins/FHCAL.cc)
//         6  LFHCAL                 forward longitudinal HCal (plugins/FHCAL.cc)
//         7  EcalFarForwardZDC      ZDC EM section            (plugins/ZDC.cc)
//         8  HcalFarForwardZDC      ZDC hadronic section      (plugins/ZDC.cc)
//     Only LumiSpec remains among the anticipated calos (not yet wired).

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

  // Per-collection calibration constants [ns], subtracted AFTER the r/c
  // correction; order follows the wiring lists in eventbuilder.cc, missing
  // entries = 0. This is the home for the deterministic, particle-INdependent
  // offsets r/c cannot model — calorimeter shower-depth/light-collection
  // delays (+2..+5 ns measured for the EM calos), Cherenkov photon paths,
  // readout offsets. Measure with `make gatecheck` (truth-referenced peak
  // offsets), then set e.g.:
  //   -PCalTimeAlignment:offsets=3.6,3.6,1.9,3.3,...
  typename Base::template Parameter<std::vector<float>> m_offsets{this, "offsets", {},
      "per-collection time offsets [ns] subtracted after r/c (wiring order; "
      "empty/short = zeros)"};

  void Configure() {
    if constexpr (kDeposits)
      m_ts_buffer = this->GetApplication()->template GetService<TimesliceBuffer_service>();
  }

  void ChangeRun(int32_t) {}

  // Apply the r/c propagation correction to a hit and return a mutable clone.
  // 1/c = 0.0033356 ns/mm. (Was 0.0034 — a 2% error that biased barrel-TOF
  // times by ~40 ps at r=640 mm, larger than the TOF's own 30 ps sigma.)
  static constexpr double kInvC_ns_per_mm = 1.0 / 299.792458;
  MutHitT correct(const HitT& hit, double offset_ns) const {
    MutHitT c = hit.clone();
    const double x = hit.getPosition()[0];
    const double y = hit.getPosition()[1];
    const double z = hit.getPosition()[2];
    c.setTime(hit.getTime() - std::sqrt(x*x + y*y + z*z) * kInvC_ns_per_mm
              - offset_ns);
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

      const double off = (i < m_offsets().size()) ? double(m_offsets()[i]) : 0.0;
      std::vector<MutHitT> current;
      if (in != nullptr) {
        current.reserve(in->size());
        for (const auto& hit : *in)
          current.push_back(correct(hit, off));
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

// TrkTimeAlignment_factory: aligns raw hit times (pre-tracking). NFast=6 means
// indices 0-5 (TOF, MPGD) are fast and indices 6-9 (Si/B0) are slow — their
// corrected hits are deposited into the TimesliceBuffer_service for
// cross-frame recovery by the EventUnfolder.
using TrkTimeAlignment_factory = TimeAlignment_factory<edm4eic::TrackerHit, 6>;

// CalTimeAlignment_factory: all calo clusters treated as fast (no deposit).
using CalTimeAlignment_factory = TimeAlignment_factory<edm4eic::Cluster>;
