// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// Time-alignment factory, shared by trackers and calorimeter clusters.
//
// For each input collection, this clones every hit or cluster, subtracts a
// straight-line propagation term r/c (so all times refer to a common
// collision t0), and sorts the collection by time. The logic is identical
// for tracker hits (edm4eic::TrackerHit) and calo clusters
// (edm4eic::Cluster), so it is written once as a template and instantiated
// once per family.
//
// This factory is stateless by design. JANA2 creates one factory instance
// per in-flight timeslice (jana:max_inflight_timeslices, default nthreads)
// and assigns frames to instances arbitrarily. Factory members must never
// carry information across frames: with nthreads > 1, such state would
// fragment across instances and see frames out of order.
//
// Cross-frame boundary recovery (slow detectors only)
// ---------------------------------------------------
// Si/MAPS detectors have a resolution of about 2 microseconds, close to
// the frame width, so hits near a frame boundary can be reconstructed into
// a neighboring frame. For slow detectors (indices >= NFast), this factory
// deposits each frame's corrected hits into the process-wide
// TimesliceBuffer_service, keyed by absolute frame number. EventUnfolder
// (see eventbuilder.cc) fetches the adjacent frames from that service
// while gating hits per candidate; the window is controlled there by
// -Peventbuilder:unfolder:backward_frames and
// -Peventbuilder:unfolder:forward_frames. The *TimeAlignRecHits output
// itself always contains only the current frame.
//
// Parameter namespaces: gating and cross-frame recovery are the unfolder's
// business (eventbuilder:unfolder:*), while what MC truth gets written is
// its own concern (eventbuilder:truth:*). This factory reads
// eventbuilder:truth:simhits, which is legitimate under that split -- the
// switch describes the truth output, not the unfolding, and both this
// factory and the unfolder act on it.
//
// Aliases are defined at the bottom of this file. The index order is fixed
// by the wiring lists in eventbuilder.cc (m_simtrackerhit_collection_names
// and m_simcalocluster_collection_names). Keep all three in sync.
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
//     Future candidates, commented out in the wiring: TaggerTracker (low
//     Q2), DIRC bars, DRICH. Appending keeps existing indices stable; a new
//     slow detector must go at the end, and NFast stays 6 unless the new
//     detector is a genuine trigger detector.
//
//   CalTimeAlignment_factory = TimeAlignment_factory<edm4eic::Cluster>
//     Same r/c correction for calorimeter clusters. 9 collections wired
//     today. NFast defaults to UINT32_MAX: every cluster is "fast", so
//     nothing is deposited.
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

#include <unordered_map>

#include <extensions/jana/JOmniFactory.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <podio/LinkCollection.h>

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

  using SlowLink = podio::Link<edm4eic::RawTrackerHit, edm4hep::SimTrackerHit>;

  // Truth links for the slow detectors, one per slow index (NFast..NFast+3).
  // Read here, and not in the EventUnfolder, because this is the last point
  // at which a frame's own link -> SimTrackerHit -> MCParticle chain is
  // still alive: the unfolder sees adjacent frames only through
  // TimesliceBuffer_service, long after their collections are gone.
  //
  // Four fixed (non-variadic) inputs rather than one variadic input: this
  // JOmniFactory shim splits the wiring's input-name list evenly across
  // however many variadic inputs a factory declares, so a second variadic
  // input would silently take half of m_in's RecHit names. Non-variadic
  // inputs each consume exactly one trailing name instead, which leaves
  // m_in's share intact. They are appended to both wiring lists in
  // eventbuilder.cc; the calo wiring passes empty names, and being optional
  // they then simply read back as nullptr.
  static constexpr size_t kNumSlowLink = 4;
  typename Base::template PodioInput<SlowLink, true> m_slow_link0{this};
  typename Base::template PodioInput<SlowLink, true> m_slow_link1{this};
  typename Base::template PodioInput<SlowLink, true> m_slow_link2{this};
  typename Base::template PodioInput<SlowLink, true> m_slow_link3{this};

  const typename PodioTypeMap<SlowLink>::collection_t* slowLinks(size_t si) {
    switch (si) {
    case 0:  return m_slow_link0();
    case 1:  return m_slow_link1();
    case 2:  return m_slow_link2();
    case 3:  return m_slow_link3();
    default: return nullptr;
    }
  }

  std::shared_ptr<TimesliceBuffer_service> m_ts_buffer;

  // (collectionID, index) of a raw hit, matching EventUnfolder::rawHitKey.
  template <typename ObjT>
  static uint64_t rawHitKey(const ObjT& h) {
    const auto id = h.getObjectID();
    return (static_cast<uint64_t>(static_cast<uint32_t>(id.collectionID)) << 32) |
           static_cast<uint32_t>(id.index);
  }

  // Per-collection calibration constants, in ns, subtracted after the r/c
  // correction. Order follows the wiring lists in eventbuilder.cc; missing
  // entries default to 0. This holds the deterministic, particle-independent
  // offsets that r/c cannot model: calorimeter shower-depth and
  // light-collection delays, Cherenkov photon paths, readout offsets.
  // Measure with `make gatecheck` (truth-referenced peak offsets), then set,
  // for example:
  //   -PCalTimeAlignment:offsets=3.6,3.6,1.9,3.3,...
  typename Base::template Parameter<std::vector<float>> m_offsets{this, "offsets", {},
      "per-collection time offsets [ns] subtracted after r/c (wiring order; "
      "empty/short = zeros)"};

  // Buffer the full truth chain (sim hit + MCParticle) alongside each slow
  // hit, so cross-frame-recovered hits can resolve
  // link -> SimTrackerHit -> MCParticle in the child event. Costs memory in
  // every retained frame, so it follows the unfolder's own switch rather than
  // being always on: the weight-encoded label does not need it.
  bool m_buffer_truth = false;

  void Configure() {
    if constexpr (kDeposits) {
      m_ts_buffer = this->GetApplication()->template GetService<TimesliceBuffer_service>();
      // Read the unfolder's parameter rather than declaring a second one, so
      // there is a single switch. Reading an unset parameter throws, so check.
      auto* pm = this->GetApplication()->GetJParameterManager();
      if (pm->Exists("eventbuilder:truth:simhits"))
        m_buffer_truth =
            this->GetApplication()->template GetParameterValue<bool>("eventbuilder:truth:simhits");
    }
  }

  void ChangeRun(int32_t) {}

  // Applies the r/c propagation correction to a hit and returns a mutable
  // clone. 1/c = 0.0033356 ns/mm. This constant must stay precise: a 2%
  // error here would bias barrel-TOF times by about 40 ps at r=640 mm,
  // larger than TOF's own 30 ps resolution.
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

  void Process(int32_t, uint64_t frame_nr) {
    const size_t n_out  = m_out().size();
    const size_t nf     = static_cast<size_t>(NFast);
    const size_t n_slow = (nf < n_out) ? (n_out - nf) : 0u;

    TimesliceBuffer_service::FrameHits slow_frame;
    if constexpr (kDeposits)
      slow_frame.resize(n_slow);

    for (size_t i = 0; i < n_out; ++i) {
      // Bounds-guarded, like m_offsets() just below: a variadic input list can
      // be SHORTER than the output list when a wired collection has no
      // producer (EcalLumiSpec has no plugin yet), and an unguarded .at(i)
      // then throws std::out_of_range out of the factory. The nullptr branch
      // below already handles "no input for this slot", so treat a missing
      // entry the same way instead of aborting the job.
      const auto* in  = (i < m_in().size()) ? m_in().at(i) : nullptr;
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
          // untracked and go to the cross-frame deposit below, together
          // with the truth resolved here (see SlowHit's comment for why
          // the unfolder cannot do this itself).
          for (const auto& h : current)
            out->push_back(h.clone());

          const size_t si = i - nf;

          // raw-hit key -> (generatorStatus, MCParticle time) for this frame.
          std::unordered_map<uint64_t, std::pair<float, float>> truth;
          // Optional full chain, deduplicated within this frame: one sim hit
          // backs several raw hits and one particle backs many sim hits, so
          // clone each object once and share the handle.
          struct Chain {
            edm4hep::MutableSimTrackerHit sim;
            edm4hep::MutableMCParticle    particle;
            uint64_t sim_key = 0;
            uint32_t mcp_key = 0;
            bool     has_particle = false;
          };
          std::unordered_map<uint64_t, Chain>                        chain_by_raw;
          std::unordered_map<uint64_t, edm4hep::MutableSimTrackerHit> sim_clones;
          std::unordered_map<uint32_t, edm4hep::MutableMCParticle>    mcp_clones;

          if (const auto* links = slowLinks(si); links != nullptr) {
            for (const auto& l : *links) {
              const auto sim = l.getTo();
              if (!l.getFrom().isAvailable() || !sim.isAvailable() ||
                  !sim.getParticle().isAvailable())
                continue;
              const auto p = sim.getParticle();
              const uint64_t rk = rawHitKey(l.getFrom());
              // Secondaries inherit their first labelled ancestor's
              // (status, time) -- see effectiveProvenance().
              const auto [est, etime] = effectiveProvenance(p);
              truth.emplace(rk, std::make_pair(static_cast<float>(est),
                                               static_cast<float>(etime)));
              if (!m_buffer_truth)
                continue;
              const uint64_t skey = rawHitKey(sim);
              const uint32_t pkey = static_cast<uint32_t>(p.getObjectID().index);
              auto sit = sim_clones.find(skey);
              if (sit == sim_clones.end())
                // relations dropped: they point into this frame, which the
                // consumer never sees. The consumer re-links sim -> particle
                // from the buffered clones instead.
                sit = sim_clones.emplace(skey, sim.clone(false)).first;
              auto pit = mcp_clones.find(pkey);
              if (pit == mcp_clones.end())
                pit = mcp_clones.emplace(pkey, p.clone(false)).first;
              chain_by_raw.emplace(rk, Chain{sit->second, pit->second, skey, pkey, true});
            }
          }

          std::vector<TimesliceBuffer_service::SlowHit> slow;
          slow.reserve(current.size());
          for (auto& h : current) {
            TimesliceBuffer_service::SlowHit sh;
            sh.frame = frame_nr;   // part of the dedup key downstream; see SlowHit
            sh.hit = std::move(h);
            const auto rh = sh.hit.getRawHit();
            if (rh.isAvailable()) {
              // Untracked clone: the producing frame's RawTrackerHit
              // collection does not outlive this deposit.
              sh.raw     = rh.clone(false);
              sh.has_raw = true;
              const uint64_t rk = rawHitKey(rh);
              auto it = truth.find(rk);
              if (it != truth.end()) {
                sh.status    = it->second.first;
                sh.ptime     = it->second.second;
                sh.has_truth = true;
              }
              if (m_buffer_truth) {
                auto cit = chain_by_raw.find(rk);
                if (cit != chain_by_raw.end()) {
                  sh.sim          = cit->second.sim;
                  sh.sim_key      = cit->second.sim_key;
                  sh.has_sim      = true;
                  sh.particle     = cit->second.particle;
                  sh.mcp_key      = cit->second.mcp_key;
                  sh.has_particle = cit->second.has_particle;
                }
              }
            }
            slow.push_back(std::move(sh));
          }
          slow_frame[si] = std::move(slow);
          continue;
        }
      }
      for (auto& h : current)
        out->push_back(h);
    }

    if constexpr (kDeposits) {
      // Deposit unconditionally, even when this frame has no slow hits at all.
      // The buffer answers "is frame N coming?" from a contiguous watermark
      // over deposited frame numbers, and a frame that never deposits stalls
      // that watermark forever -- every later fetch would fall back to the
      // timeout, which is the nondeterminism the watermark exists to remove.
      // An empty deposit is a few bytes and says "frame N happened, and it
      // had nothing for you".
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
