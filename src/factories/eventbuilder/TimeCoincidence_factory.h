// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// Candidate-t0 coincidence gate, shared by tracker hits and calo clusters.
//
// Sits between EventBuilder_factory and EventUnfolder. For each input hit
// or cluster, it forwards the element only if it falls within any
// EventCandidate's acceptance window:
//
//   tracker: |t_hit - t0_candidate| <= nsigma x sigma_hit (symmetric)
//   calo:    |t_hit - t0_candidate| <= nsigma x sigma_hit, with the + side
//            widened further by late_nsigma_window_asym x sigma_hit
//
// sigma_hit is the element's own getTimeError(): per-hit or per-cluster,
// never a group-level constant.
//
// late_nsigma_window_asym is calo only (see kHasEnergyGate below).
// Shower-development and light-collection delay push calo hit times
// systematically late relative to a straight-line r/c correction, an
// asymmetry confirmed by `make rescheck`'s time-difference panels
// (external resolution-check tooling). The same panels, run on tracker
// families, show no such asymmetry, so the tracker gate stays plain and
// symmetric, with no late-side widening. late_nsigma_window_asym is a
// sigma multiple, not a flat ns constant, since a flat ns add-on would not
// be comparable across calo families of very different resolution.
//
// The logic is identical for tracker hits (edm4eic::TrackerHit) and calo
// clusters (edm4eic::Cluster), apart from the late-side widening and the
// energy floor below, so it is written once as a template and instantiated
// per family below. Both calo-only extras are compiled out entirely for
// the tracker instantiation, through `if constexpr`, not just skipped at
// runtime.
//
// This factory is stateless by design. JANA2 creates one factory instance
// per in-flight timeslice and assigns frames to instances arbitrarily, so
// factory members must never carry information across frames.
//
// Output: PODIO subset collections of the input (already time-aligned)
// collections. All elements originate from the current frame, so subset
// references are valid.
//
// Cross-frame boundary recovery, for slow tracker detectors only, is
// handled downstream by EventUnfolder (eventbuilder.cc), which merges the
// adjacent frames' hits from TimesliceBuffer_service. This factory gates
// the current frame only.
//
// Energy floor (calo only)
// -------------------------
// A cluster must pass both conditions:
//   1. |t_cluster - t0_candidate| <= the time gate above
//   2. cluster.energy >= min_energy
// min_energy is a per-cluster floor, not a background-rate measurement: it
// is set to the EIC Yellow Report's minimum detectable-photon energy
// (arXiv:2103.05419): 50 MeV near-barrel (|eta|<1), 30 MeV
// forward/backward (1<=|eta|<4). Below that floor, a cluster is not a
// physically real single deposit, regardless of any background
// measurement, so this floor does not need the retuning a summed-quantity
// threshold would (compare EventBuilder_factory's
// eventbuilder:trigger:min_cal_energy, which floors e_calo, a sum over
// clusters, and does need a measured background scan).
//
// Barrel-ness is read off the cluster's own detector index (is_barrel_det
// below), not a cluster-level eta cut, since clusters carry no eta
// directly and the detector split is already an exact proxy for the YR's
// angular regions here. debug_gating (calo only) logs the per-detector,
// per-frame cutflow (usable sigma, then energy floor, then time gate), so
// a retune can see exactly which cut removes clusters.
//
// Aliases are defined at the bottom of this file. Index order is fixed by
// the wiring lists in eventbuilder.cc; keep them in sync.
//
//   TrkTimeCoincidence = TimeCoincidence_factory<edm4eic::TrackerHit>
//     Wired twice, under different tags, same type, different i/o:
//       tag `trkcoincidence`: fast hits (TOF+MPGD). Feeds EventUnfolder's
//         m_trk_in directly (indices 0-5, *TrkCoincRecHits).
//       tag `coincidence`: slow hits (Si/B0). Feeds EventUnfolder's
//         m_trk_in (indices 6-9, *TimeCoincRecHits), but not as a direct
//         passthrough like the fast-hit case above: EventUnfolder re-gates
//         these a second time inline (its own trk_nsigma_window, through
//         slowHitWindow()), so two nsigma windows apply in series. A hit
//         dropped here can never be recovered downstream, so this gate
//         must stay at least as permissive as EventUnfolder's own per-hit
//         gate.
//
//   CalTimeCoincidence = TimeCoincidence_factory<edm4eic::Cluster>
//     Tag `calcoincidence`. Feeds EventUnfolder's m_clu_in directly
//     (*CoincClusters), already gated to candidate t0 and to energy, so
//     EventUnfolder's e_calo sum needs no further inline gate.
//
// Parameters (registered per instance, under that instance's own tag)
//   nsigma                   N: gate half-width = N x sigma (per-element
//                             getTimeError())
//   late_nsigma_window_asym  calo only: extra LATE acceptance, in units of
//                             sigma, added to the + side only. No effect on
//                             TrkTimeCoincidence, whose gate stays plain
//                             and symmetric.
//   min_energy_barrel        calo only: min cluster energy (GeV), for
//                             EcalBarrel/HcalBarrel clusters. No effect on
//                             TrkTimeCoincidence.
//   min_energy_forward       calo only: min cluster energy (GeV), for every
//                             other calo collection. No effect on
//                             TrkTimeCoincidence.
//   debug_gating              calo only: 1 = log per-detector, per-frame
//                             cutflow. No effect on TrkTimeCoincidence.

#pragma once

#include <cmath>
#include <type_traits>
#include <vector>

#include <JANA/JLogger.h>
#include <extensions/jana/JOmniFactory.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/ClusterCollection.h>

template <typename HitT>
struct TimeCoincidence_factory : public JOmniFactory<TimeCoincidence_factory<HitT>> {

  using Base = JOmniFactory<TimeCoincidence_factory<HitT>>;

  // Only clusters carry an energy field or need the energy-floor cutflow:
  // edm4eic::TrackerHit has no getEnergy(). Compiled out, not just skipped,
  // for the tracker instantiation.
  static constexpr bool kHasEnergyGate = std::is_same_v<HitT, edm4eic::Cluster>;

  // EventCandidates supplies the t0 values for the coincidence gate.
  typename Base::template PodioInput<edm4hep::EventHeader> m_candidates_in{this, "EventCandidates"};

  // Time-aligned input (current frame only).
  typename Base::template VariadicPodioInput<HitT, true, true> m_hits_in{this, {}};

  // One subset output per detector.
  typename Base::template VariadicPodioOutput<HitT> m_hits_out{this, {}};

  typename Base::template Parameter<float> nsigma{
      this, "nsigma", 3.0f,
      "N-sigma half-width: gate = N x sigma (per-element getTimeError())."};
  // Calo only (see kHasEnergyGate): shower-development and light-collection
  // delay push calo hit times systematically late, an asymmetry measured
  // in res_check.py's time-difference panels. Declared unconditionally, so
  // both instantiations share one definition; TrkTimeCoincidence never
  // reads it, since its own residuals measured symmetric (see the file
  // header).
  typename Base::template Parameter<float> late_nsigma_window_asym{
      this, "late_nsigma_window_asym", 0.5f,
      "CALO ONLY: extra LATE acceptance, in units of that element's own sigma "
      "(getTimeError()), added to the + side only. No effect on the tracker "
      "instantiation (measured symmetric, no late widening applied there). "
      "UNMEASURED PLACEHOLDER — retune via an efficiency/purity scan."};

  typename Base::template Parameter<float> min_energy_barrel{
      this, "min_energy_barrel", 0.050f,
      "CALO ONLY: minimum cluster energy (GeV) to pass for EcalBarrel/HcalBarrel "
      "clusters -- YR near-barrel (|eta|<1) minimum detectable-photon energy. "
      "No effect on the tracker instantiation."};
  typename Base::template Parameter<float> min_energy_forward{
      this, "min_energy_forward", 0.030f,
      "CALO ONLY: minimum cluster energy (GeV) to pass for every other calo "
      "collection -- YR forward/backward (1<=|eta|<4) minimum detectable-photon "
      "energy. No effect on the tracker instantiation."};
  typename Base::template Parameter<bool> debug_gating{
      this, "debug_gating", false,
      "CALO ONLY: 1 = log per-detector, per-frame cluster survival at each cut "
      "stage in sequence (usable sigma -> energy floor -> time gate). No effect "
      "on the tracker instantiation."};

  void Configure() {}
  void ChangeRun(int32_t) {}

  // Detector index to YR region: true = near-barrel, false =
  // forward/backward. Indices match eventbuilder.cc's calcoinc_inputs
  // wiring order (B0ECal, EcalBarrel, EcalEndcapN, EcalEndcapP, HcalBarrel,
  // HcalEndcapPInsert, LFHCAL, EcalFarForwardZDC, HcalFarForwardZDC,
  // EcalLumiSpec). Only indices 1 (EcalBarrel) and 4 (HcalBarrel) are barrel.
  static bool is_barrel_det(size_t s) { return s == 1 || s == 4; }

  void Process(int32_t, uint64_t) {
    const size_t n_det = m_hits_out().size();
    const float  n     = nsigma();
    // Tracker gate is plain symmetric (see the file header); late-side
    // widening is calo-only.
    const float  late  = kHasEnergyGate ? late_nsigma_window_asym() : 0.0f;

    const auto* cands = m_candidates_in();
    std::vector<float> t0s;
    if (cands)
      for (const auto& h : *cands)
        t0s.push_back(static_cast<float>(h.getWeights(0)));

    auto passes_time = [&](float t, float sigma) -> bool {
      const float half_win = n * sigma;
      for (const float t0 : t0s) {
        const float dt = t - t0;
        if (dt >= -half_win && dt <= half_win + late * sigma) return true;
      }
      return false;
    };

    for (size_t s = 0; s < n_det; ++s) {
      auto& out = m_hits_out().at(s);
      // Subset collection: every element originates from the current frame's
      // input collection (produced by TrkTimeAlignment_factory/
      // CalTimeAlignment_factory).
      out->setSubsetCollection(true);
      const auto* in = m_hits_in().at(s);
      if (!in) continue;

      if constexpr (kHasEnergyGate) {
        const float e_min = is_barrel_det(s) ? min_energy_barrel() : min_energy_forward();
        const bool dbg = debug_gating();
        // Cut stages are counted in sequence: each stage's denominator is
        // the previous stage's survivors, not the full input. A true
        // cutflow, not independent counts.
        int n_in = 0, n_sigma_ok = 0, n_energy_ok = 0, n_time_ok = 0;
        for (const auto& hit : *in) {
          ++n_in;
          const float sigma = hit.getTimeError();
          if (!(sigma > 0.0f)) continue;
          ++n_sigma_ok;
          if (!(hit.getEnergy() >= e_min)) continue;
          ++n_energy_ok;
          if (!passes_time(hit.getTime(), sigma)) {
            if (dbg) {
              float best_dt = 1e30f;
              for (const float t0 : t0s) {
                const float d = hit.getTime() - t0;
                if (std::fabs(d) < std::fabs(best_dt)) best_dt = d;
              }
              jout << "[eb:calcoinc:miss] det=" << s << " t_clu=" << hit.getTime()
                   << " sigma=" << sigma << " nearest_t0_dt=" << best_dt
                   << " half_win=" << (n * sigma) << jendl;
            }
            continue;
          }
          ++n_time_ok;
          out->push_back(hit);
        }
        if (dbg)
          jout << "[eb:calcoinc:gating] det=" << s
               << " in=" << n_in << " sigma_ok=" << n_sigma_ok
               << " energy_ok=" << n_energy_ok << " time_ok=" << n_time_ok << jendl;
      } else {
        for (const auto& hit : *in) {
          const float sigma = hit.getTimeError();
          if (sigma > 0.0f && passes_time(hit.getTime(), sigma))
            out->push_back(hit);
        }
      }
    }
  }
};

using TrkTimeCoincidence = TimeCoincidence_factory<edm4eic::TrackerHit>;
using CalTimeCoincidence = TimeCoincidence_factory<edm4eic::Cluster>;
