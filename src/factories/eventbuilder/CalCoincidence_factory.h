// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
//
// CalCoincidence_factory
// -----------------------
// Candidate-t0 gate for calo clusters, same mechanism as TimeCoincidence_
// factory/TrkCoincidence_factory, PLUS a minimum-energy floor — the two
// conditions a cluster must pass, both required:
//
//   1. |t_cluster - t0_candidate| <= nsigma × sigma_cluster  (time)
//   2. cluster.energy >= min_energy                          (energy)
//
// Today's E_calo (EventBuilder_factory::emitCandidate, weights[24]) applies
// condition 1 only, inline, collapsed immediately into a summed scalar — no
// energy floor exists anywhere on calo clusters today (checked directly:
// EventBuilder_factory.h's e_calo accumulation sums every cluster in the
// window regardless of size). This factory is the same gate promoted to a
// real, reusable, first-class collection instead of a scalar side effect —
// so downstream consumers (currently EventPrefilter_factory) see the
// coincidence-cleaned clusters themselves, not the frame-wide *TimeAlign
// Clusters mirror.
//
// min_energy is a PER-CLUSTER floor, not a background-rate measurement: it's
// set to the EIC Yellow Report's minimum detectable-photon energy
// (arXiv:2103.05419) — 50 MeV near-barrel (|eta|<1), 30 MeV forward/backward
// (1<=|eta|<4), same split data_stats.py already uses for the STAGE 2 calo
// truth denominator. Below that floor a cluster isn't a physically real
// single deposit regardless of any background measurement, so this doesn't
// need retuning the way a summed-quantity threshold would (see
// EventBuilder_factory's eventbuilder:trigger:min_cal_energy, which floors
// e_calo — a SUM over clusters — and does need a measured background scan;
// applying a single-photon floor to a multi-cluster sum is the wrong unit).
// Barrel-ness is read off the cluster's own detector index (see is_barrel_det
// below) rather than a cluster-level eta cut: clusters don't carry eta
// directly, and the detector split is already an exact proxy for the YR's
// angular regions here.
//
// Output: *Clusters — PODIO subset collections of *TimeAlignClusters.
// Not currently wired into EventUnfolder (which still reads raw
// *TimeAlignClusters and applies its own cal_nsigma_window+
// cal_late_nsigma_window_asym gate directly) — only EventPrefilter_factory
// consumes this today. late_nsigma_window_asym is still offered here,
// mirroring TimeCoincidence_factory/EventUnfolder, so behavior stays
// consistent if that wiring ever changes. It is a SIGMA MULTIPLE, not a flat
// ns constant — see TimeCoincidence_factory.h's header comment for why a
// flat ns add-on isn't comparable across detectors/clusters of very
// different resolution.
//
// Parameters
//   nsigma                   N: gate half-width = N × sigma_cluster
//                             (getTimeError())
//   late_nsigma_window_asym  extra LATE acceptance, in units of sigma_cluster,
//                             added to the + side of the time gate only
//                             (0 = legacy symmetric gate)
//   min_energy_barrel         minimum cluster energy (GeV) to pass, EcalBarrel/
//                             HcalBarrel clusters (YR near-barrel |eta|<1
//                             minimum detectable photon energy).
//   min_energy_forward        minimum cluster energy (GeV) to pass, every other
//                             calo collection (YR forward/backward endcap
//                             1<=|eta|<4 minimum detectable photon energy).

#pragma once

#include <cmath>
#include <vector>

#include <JANA/JLogger.h>
#include <extensions/jana/JOmniFactory.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4eic/ClusterCollection.h>

struct CalCoincidence_factory : public JOmniFactory<CalCoincidence_factory> {

  // EventCandidates supplies the t0 values for the coincidence gate.
  PodioInput<edm4hep::EventHeader> m_candidates_in{this, "EventCandidates"};

  // Calo *TimeAlignClusters (current frame only).
  VariadicPodioInput<edm4eic::Cluster, true> m_calo_in{this, {}};

  // One subset output per calo detector.
  VariadicPodioOutput<edm4eic::Cluster> m_calo_out{this, {}};

  Parameter<float> nsigma{
      this, "nsigma", 3.0f,
      "N-sigma half-width: gate = N × sigma_cluster (per-cluster getTimeError())."};
  Parameter<float> late_nsigma_window_asym{
      this, "late_nsigma_window_asym", 0.5f,
      "extra LATE acceptance, in units of that cluster's own sigma "
      "(getTimeError()), added to the + side of the time gate only (mirrors "
      "EventUnfolder's eventbuilder:unfolder:cal_late_nsigma_window_asym). 0 = legacy "
      "symmetric gate. UNMEASURED PLACEHOLDER — retune via an efficiency/"
      "purity scan."};
  Parameter<float> min_energy_barrel{
      this, "min_energy_barrel", 0.050f,
      "minimum cluster energy (GeV) to pass for EcalBarrel/HcalBarrel clusters "
      "-- YR near-barrel (|eta|<1) minimum detectable-photon energy."};
  Parameter<float> min_energy_forward{
      this, "min_energy_forward", 0.030f,
      "minimum cluster energy (GeV) to pass for every other calo collection "
      "-- YR forward/backward endcap (1<=|eta|<4) minimum detectable-photon "
      "energy."};
  Parameter<bool> debug_gating{
      this, "debug_gating", false,
      "1 = log per-detector, per-frame cluster survival at each cut stage "
      "in sequence (usable sigma -> energy floor -> time gate) -- lets you "
      "see exactly which cut is removing clusters, instead of just in/out "
      "totals. 0 (default) = silent."};

  void Configure() {}
  void ChangeRun(int32_t) {}

  // Detector index -> YR "near-barrel" (true) vs "forward/backward" (false).
  // Indices match eventbuilder.cc's calcoinc_inputs wiring order (B0ECal,
  // EcalBarrel, EcalEndcapN, EcalEndcapP, HcalBarrel, HcalEndcapPInsert,
  // LFHCAL, EcalFarForwardZDC, HcalFarForwardZDC, EcalLumiSpec) -- only
  // indices 1 (EcalBarrel) and 4 (HcalBarrel) are barrel.
  static bool is_barrel_det(size_t s) { return s == 1 || s == 4; }

  void Process(int64_t, uint64_t) {
    const size_t n_det = m_calo_out().size();
    const float  n     = nsigma();
    const float  late  = late_nsigma_window_asym();
    const float  e_min_bar = min_energy_barrel();
    const float  e_min_fwd = min_energy_forward();

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

    const bool dbg = debug_gating();
    for (size_t s = 0; s < n_det; ++s) {
      auto& out = m_calo_out().at(s);
      // Subset collection: every cluster originates from the current
      // frame's *TimeAlignClusters collection produced by
      // CalTimeAlignment_factory.
      out->setSubsetCollection(true);
      const auto* in = m_calo_in().at(s);
      if (!in) continue;
      const float e_min = is_barrel_det(s) ? e_min_bar : e_min_fwd;
      // Cut stages counted IN SEQUENCE (each stage's denominator is the
      // previous stage's survivors), not independently -- so this is a true
      // cutflow, not three separate pass rates against the same n_in.
      int n_in = 0, n_sigma_ok = 0, n_energy_ok = 0, n_time_ok = 0;
      for (const auto& clu : *in) {
        ++n_in;
        const float sigma = clu.getTimeError();
        if (!(sigma > 0.0f)) continue;
        ++n_sigma_ok;
        if (!(clu.getEnergy() >= e_min)) continue;
        ++n_energy_ok;
        if (!passes_time(clu.getTime(), sigma)) {
          if (dbg) {
            float best_dt = 1e30f;
            for (const float t0 : t0s) {
              const float d = clu.getTime() - t0;
              if (std::fabs(d) < std::fabs(best_dt)) best_dt = d;
            }
            jout << "[eb:calcoinc:miss] det=" << s << " t_clu=" << clu.getTime()
                 << " sigma=" << sigma << " nearest_t0_dt=" << best_dt
                 << " half_win=" << (n * sigma) << jendl;
          }
          continue;
        }
        ++n_time_ok;
        out->push_back(clu);
      }
      if (dbg)
        jout << "[eb:calcoinc:gating] det=" << s
             << " in=" << n_in << " sigma_ok=" << n_sigma_ok
             << " energy_ok=" << n_energy_ok << " time_ok=" << n_time_ok << jendl;
    }
  }
};
