// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// EventBuilder_factory
// -------------------
// Streaming "event building" stage. Runs once per time-frame (JANA Timeslice
// level) over the time-aligned tracker hits and *locates* candidate physics
// events. For each candidate it produces a precise event time t0 and its
// half-width t0sigma.
//
// t0  = inverse-variance weighted mean of the coincident trigger-detector
//       hit times,  t0 = Sum(t_i / sigma_i^2) / Sum(1 / sigma_i^2)
// t0sigma = nsigma_window * sigma_t0,  sigma_t0 = 1 / sqrt(Sum 1/sigma_i^2)
//
// The SAME nsigma_window also derives the coincidence windows that FIND a
// candidate in the first place (dt, dt_tof -- see Process()): each is
// nsigma_window * sqrt(2) * the worst per-hit resolution actually present in
// that pool this frame, not a fixed constant.
//
// Per-hit sigma_i is edm4eic::TrackerHit::getTimeError(), which each detector's
// own TrackerHitReconstruction_factory (or equivalent) already populates from
// its digitization config's timeResolution. This factory does NOT keep its own
// copy of that number: reading it off the hit means it can never drift out of
// sync with what the detector plugin actually configured.
//
// TOF (~30 ps) dominates the weighting, so t0 is effectively set by the fastest
// detector. Silicon (~2 us) is carried by the detector but is NOT used to find
// events (it cannot time-resolve them); only TOF + MPGD drive the coincidence.
//
// This factory does assembly/finding only -- NO reconstruction and it does NOT
// create child events. The downstream EventUnfolder turns each surviving
// candidate into a JANA PhysicsEvent.
//
// Output: an edm4hep::EventHeaderCollection "EventCandidates", one entry per
// candidate, ordered in time. Encoding (lossless):
//   timeStamp  = llround(t0_ns * 1000)   (t0 in ps, integer)
//   weight     = t0sigma_ns
//   weights[0] = t0_ns   weights[1] = t0sigma_ns   weights[2] = phys/fake flag
//                (1 = matched the tagged primary collision, 2 = matched
//                nothing (accidental), 3 = matched a pileup collision --
//                a genuine OTHER real ep collision, cross-mixed in under
//                PILEUP_FREQ, just not the one tagged primary; see the
//                pileup_times block in Process(). weights[30] = the number
//                of pileup collisions injected into this candidate's frame
//                (same value on every candidate of one frame) -- together
//                with a count of flag==3 candidates this gives pileup
//                efficiency/purity, parallel to the flag==1 pair above.
//                weights[35] = n_mc_coincident, the PER-CANDIDATE count of real
//                collisions (primary + pileup) actually within +-dt of
//                THIS candidate's own t0 -- unlike weights[2]/[5], which
//                only ever record ONE collision per candidate (nearest
//                wins, matching is exclusive), n_mc_coincident>=2 flags genuine
//                combined pile-up: multiple real collisions landing in the
//                SAME candidate window, which flag/mc_t cannot represent.
//                weights[36] = coincident_streams, a bitmask of WHICH pileup
//                class(es) (status/1000, see vendor/sro/README.md's "Signal
//                pileup" table) are behind that count -- decode against that
//                table to name them (e.g. bit 14 = dvcs).

#pragma once

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <string>
#include <vector>

#include "TMath.h"

#include <extensions/jana/JOmniFactory.h>

#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/ClusterCollection.h>

struct EventBuilder_factory : public JOmniFactory<EventBuilder_factory> {

  // Single confidence-level knob for the whole trigger side (the unfold-stage
  // per-hit content gate is a SEPARATE parameter, EventUnfolder's own
  // "eventbuilder:unfolder:nsigma_window" -- deliberately independent, since
  // it gates PhysicsEvent content against an already-fixed t0, a different
  // statistical object from anything below; see marco.md Q11/S6.2).
  // Two jobs, both scaling with the SAME nsigma_window:
  //   1. t0sigma = nsigma_window * sigma_t0 (reported candidate precision).
  //   2. The finding-stage coincidence windows dt/dt_tof, derived per frame
  //      in Process() from this and each pool's own measured resolution --
  //      see the derivation there. Raising this to widen the TOF-only
  //      window's confidence widens the combined window proportionally too
  //      (one knob, not two); the combined path's topology:threshold floor
  //      still protects it regardless of window width.
  float m_nsigma_window = 3.0f;
  ParameterRef<float> nsigma_window{this, "nsigma_window", m_nsigma_window,
                                    "N-sigma confidence level for t0sigma AND for deriving "
                                    "the finding-stage coincidence windows (dt = nsigma_window * "
                                    "sqrt(2) * max per-hit resolution in the pool)"};

  // Topology trigger: minimum coincident hits in one staggered (theta,phi)
  // bin to fire. Raise to cut accidental (fake) candidates from dense machine background — Poisson
  // probability of N random hits sharing one bin falls steeply with N.
  //
  // Registered by hand below (Configure()) as "eventbuilder:topology:*"
  // instead of via ParameterRef -- ParameterRef always prefixes with this
  // factory's own tag ("eventbuilder:trigger:*"), which would nest a second
  // ":"-separated word under it (e.g. "eventbuilder:trigger:topology:threshold",
  // 4 segments). "topology" is conceptually its own component of the
  // canonical eventbuilder:<component>:<parameter> namespace, on equal
  // footing with trigger/coincidence/unfolder/prefilter, even though it is
  // implemented inside EventBuilder_factory.
  // Tracklet windows + tag thresholds (see Configure() for the parameter
  // names and semantics; score-side only, never the trigger decision).
  //
  // Why these are parameters and not derived from geometry (yet): each has a
  // different natural source, and only some are geometric —
  //   dr_min/dz_min    lever-arm floors that exclude same-layer pairs; the
  //                    true values are the actual layer spacings (MPGD barrel
  //                    r ~ 550/745 mm, TOF ~ 640; endcap planes dz ~ 125 mm)
  //                    -> DERIVABLE from the DD4hep geometry service (future).
  //   dphi_*           curvature windows: dphi ~ 0.3*B*dr / (2*pT_min) plus
  //                    multiple scattering -> derivable from the field map +
  //                    layer radii GIVEN a physics choice of pT_min; the
  //                    physics choice itself cannot come from geometry.
  //   z0_max           vertex-region half-length: a machine/optics parameter
  //                    (bunch length / IR design), not detector geometry.
  //   r_tol            projective tolerance: MS-dominated, depends on
  //                    material budget and pT spectrum -> tuned, not derived.
  //   zone_nhits/b0/zdc detector-response thresholds (shower sizes) from the
  //                    trigger-study reference; re-tune per campaign.
  // Until a geometry-service derivation exists, they are runtime parameters
  // with offline-tuned defaults — never silent hardcodes.
  double m_trk_dphi_barrel = 0.05;
  double m_trk_dphi_endcap = 0.10;
  double m_trk_dr_min      = 50.0;
  double m_trk_dz_min      = 50.0;
  double m_trk_z0_max      = 100.0;
  double m_trk_r_tol       = 20.0;
  int m_tag_zone_nhits  = 10;
  int m_tag_b0_min_hits = 4;
  int m_tag_zdc_min_hits = 50;

  // STAGE 1: TRIGGER — offline TIME && (TRK || CAL) thresholds on the
  // already-stored n_tracklets/e_calo (weights[25]/[24]). Score-only, same
  // status as tags:*/tracklet:* above: NOT applied to phys_flag or
  // candidate emission. Exist so the threshold a run was actually
  // configured with travels with the data (weights[31]/[32]) instead of
  // being hardcoded a second time, separately, in offline analysis
  // (data_stats.py's TRIGGER = TIME && (TRK || CAL) evaluation).
  int    m_trig_min_tracklets  = 1;
  double m_trig_min_cal_energy = 0.0;  // GeV; 0 = off, UNMEASURED PLACEHOLDER -- thresholds
                                        // e_calo, a SUM over clusters; needs a background scan
                                        // (labeled fake candidates' e_calo distribution), same
                                        // discipline as everything else in this plugin. Distinct
                                        // from cal_min_energy_barrel/forward below, which floor
                                        // ONE cluster at a time before it ever enters the sum.

  // Per-cluster calo gate feeding e_calo (see Process()'s calo_te loop): SAME
  // formula as CalCoincidence_factory (nsigma x cluster.getTimeError() + late
  // asym on the + side, per-cluster min_energy split by YR region), applied
  // inline here instead of consuming that factory's output because
  // CalCoincidence_factory depends on EventCandidates, which THIS factory is
  // still in the middle of producing -- see the eventbuilder README for why
  // that's a hard ordering constraint, not a style choice. Replaces the old
  // e_calo window, which reused the trigger-finding dt (built from TOF/MPGD
  // hit resolution only -- never calo-aware, see max_sigma_all above).
  double m_cal_nsigma            = 3.0;   // mirrors CalCoincidence_factory::nsigma
  double m_cal_late_nsigma_asym  = 0.5;   // mirrors CalCoincidence_factory::late_nsigma_window_asym; UNMEASURED PLACEHOLDER
  double m_cal_min_energy_barrel = 0.050; // GeV; YR near-barrel (|eta|<1) minimum detectable-photon energy
  double m_cal_min_energy_forward = 0.030; // GeV; YR forward/backward (1<=|eta|<4) minimum detectable-photon energy
  // For-later knob, default OFF: 1 = actually apply TIME && (TRK || CAL) —
  // a candidate failing BOTH thresholds above is dropped before it ever
  // becomes an EventCandidates row, so EventUnfolder never unfolds it and
  // ACTS/CaloIsland never run on it (consistent by construction — nothing
  // downstream needs separate updating). DESTRUCTIVE for any candidate that
  // would have matched a real collision (flags[i]==1, see Process()): a
  // vetoed real candidate is genuinely lost, not just hidden, and is NOT
  // recoverable without reprocessing with veto=0 — "may be accepted later"
  // means re-running with this off, not un-vetoing after the fact.
  bool m_trig_veto = false;

  int m_threshold  = 2;
  int m_theta_bins = 12;
  int m_phi_bins   = 8;

  // Per-cell background model (equal solid angle != equal expected
  // background): per-cell mu_i with the statistics ladder cell -> band ->
  // global, gated by mu_min_hits (min frame hits for an estimate to be
  // trusted). See Process() step 3.
  // t0 estimator ("mean" | "leading"):
  //   mean    — inverse-variance weighted mean of ALL window hits (default;
  //             the one-sided beta<1 late tail biases it late, ~+1 ns
  //             measured vs truth, but its breakdown tail is the smallest)
  //   leading — coincidence-confirmed leading edge: anchor on the earliest
  //             precise hit pair and average only hits time-compatible with
  //             it (the beta~1 population) — bias-free median, larger
  //             breakdown tail (A/B 2026-07). Planned: "cell" — mean over
  //             the winning topology cell's own hits.
  std::string m_t0_estimator = "mean";

  // In-acceptance definition for the expected-track count stamped into
  // weights[27]: stable charged signal MC with pT > pt_min and |eta| <
  // eta_abs. Defaults follow the EIC Yellow Report tracking matrix
  // (|eta| <= 3.5; minimum pT "100 MeV pi" — arXiv:2103.05419, detector
  // requirements table). Truth-side bookkeeping only, never a trigger cut.
  double m_acc_pt_min  = 0.1; // GeV/c
  double m_acc_eta_abs = 3.5;

  // eventbuilder:statistics:null — which null model the STORED significance
  // (weights[3]) and tail probability (weights[23]) are evaluated against,
  // named as <law>_<granularity>. The TRIGGER decision (count thresholds,
  // always the cell->band->global mu ladder) is identical in every mode.
  //   "poisson_cell"   (default) exact Poisson tail on the per-cell LOO mu
  //                    ladder — legacy; sharpest conditioning, but assumes
  //                    independent hit arrivals, which the measured
  //                    overdispersion (Fano > 1: conversion/shower hit
  //                    families) contradicts — its p-values overstate;
  //   "nbinom_band"    negative-binomial tail moment-fitted to the frame's
  //                    sideband windowed counts (per-band mean + variance) —
  //                    the count law of Poisson-arriving hit FAMILIES;
  //                    honest dispersion, band-level conditioning;
  //   "nbinom_cell"    hybrid: NB tail at the per-cell LOO mean with the
  //                    band-measured dispersion (variance = Fano_band x
  //                    mu_cell) — cell conditioning AND honest dispersion,
  //                    no sampling ceiling; degenerates to poisson_cell
  //                    wherever Fano <= 1;
  //   "empirical_band" rank of cmax against the frame's own sideband
  //                    windowed-count distribution (disjoint dt-wide
  //                    windows, empty ones included; add-one p =
  //                    (#{count >= c} + 1)/(N + 1)). Assumption-free and
  //                    self-calibrating; p floor per frame =
  //                    1/(nphi*nwindows + 1) ~ 6e-5 — pool frames for more.
  // Bare aliases: "poisson" = poisson_cell, "nbinom" = nbinom_band,
  // "empirical" = empirical_band.
  // Regardless of mode, weights[28] = empirical S and weights[29] = band
  // Fano factor are always stored, so one run carries the full A/B.
  std::string m_stat_null = "poisson_cell";

  // The trigger's background ladder granularity is fixed (cell -> band ->
  // global, self-degrading via mu_min_hits); the old topology:mu_mode knob
  // was folded into statistics:null's <law>_<granularity> naming.
  static constexpr int m_mu_mode = 2;
  int m_mu_min_hits = 50;
  // Leave-one-out scoring null (see Process(): the stored S/p_tail are
  // evaluated against mu recomputed WITHOUT the tested window). 1 = on
  // (validated: removes the self-masking circularity, S-AUC 0.770 -> 0.888
  // at 1x1); 0 = legacy frame-wide null. Runtime-switchable so the choice
  // can be revisited without a rebuild.
  int m_mu_loo = 1;

  // Optional soft-hit rejection before event finding: synchrotron/beam-gas
  // hits are predominantly low-eDep. 0 (default) = off.
  bool m_debug_calo = false;
  ParameterRef<bool> debug_calo{
      this, "debug_calo", m_debug_calo,
      "1 = dump each frame's time-aligned calo clusters (>0.1 GeV) to the log — "
      "diagnostic for cluster-time pathologies; score/trigger behavior unchanged"};

  float m_min_edep = 0.0f;
  ParameterRef<float> min_edep{
      this, "min_edep", m_min_edep,
      "minimum hit eDep (GeV) for a hit to participate in event finding; 0 = off"};

  // Adaptive-threshold trigger mode: instead of the fixed topology:threshold
  // count, derive the count threshold per frame from the frame's own background density. The
  // expected background count per (bin x window) is
  //   mu = (hits.size() x dt / frame_span) / ncells
  //   (dt = the derived combined-group coincidence window, see Process();
  //    ncells = theta_bins x phi_bins, runtime-configurable)
  // and the trigger fires when a bin's count is statistically incompatible
  // with Poisson(mu): c >= smallest k with P(X >= k | Poisson(mu)) < p_fake.
  // p_fake is the target accidental probability per (bin x window position) —
  // a fixed p_fake holds the fake rate constant across background conditions
  // (gold-coating vs vacuum vs vanilla) instead of a fixed hit count being
  // simultaneously too loose for dense frames and too tight for quiet ones.
  // The threshold is floored at topology:threshold so quiet frames never
  // trigger on a bare hit pair. The candidate's significance
  //   S = -10 * log10(p_tail),  p_tail = P(X >= cmax | Poisson(mu))
  // is a dB-scale rarity score built directly from the exact Poisson tail
  // (see poissonSurvival) — deliberately NOT a Gaussian sigma/z-score, since
  // our operating mu is routinely <1, far below the mu>~5 a Pearson-residual
  // z-score reading would need. It is computed in BOTH modes and stored in
  // the candidate (weights[3]) so the operating point can also be chosen
  // OFFLINE by cutting on S.
  // p_fake alone selects the mode: a value in (0,1) enables the adaptive
  // per-frame threshold targeting that accidental probability; 0 (the
  // default) or >=1 means fixed topology:threshold. One knob, no separate
  // enable flag.
  int m_adaptive_threshold = 0;
  ParameterRef<int> adaptive_threshold{
      this, "adaptive_threshold", m_adaptive_threshold,
      "DEPRECATED, ignored — adaptive mode is enabled by p_fake in (0,1); "
      "kept registered only so old command lines don't fail parameter strictness"};
  float m_p_fake = 0.0f;
  ParameterRef<float> p_fake{
      this, "p_fake", m_p_fake,
      "in (0,1): adaptive per-frame Poisson threshold targeting this accidental probability "
      "per (bin x window), smaller = tighter; 0 (default) or >=1: fixed topology:threshold"};

  // The single source of truth for the trigger mode (see p_fake above).
  bool adaptiveMode() {
    const double p = double(p_fake());
    return p > 0.0 && p < 1.0;
  }

  // Resolution-group refinement of the adaptive-threshold mode. The combined
  // window dt is sized for MPGD (sigma ~10 ns, see the derivation in
  // Process()) and wastes the TOF's ~300x better timing: a TOF pair inside
  // its own tight dt_tof is enormously more significant than the same pair
  // judged against dt. In adaptive mode the trigger therefore fires on EITHER
  // group: the combined group (all fast hits within dt — the existing,
  // unsubscripted test), OR the TOF group (TOF hits only within dt_tof), each
  // tested against its own Poisson background (same p_fake). Silicon has no
  // trigger group by construction — it cannot time-resolve events. Both group
  // significances are stored per candidate regardless of mode. Neither dt nor
  // dt_tof is a separate parameter: both derive per frame from the single
  // nsigma_window confidence level above (see Process()).
  // Floor on the TOF-group count threshold — the TOF-group counterpart of
  // topology:threshold, and registered next to it in Configure() as
  // "eventbuilder:topology:threshold_tof" (same manual-registration rationale
  // as the other topology:* parameters above).
  int m_threshold_tof = 2;

  // Trigger detectors = TOF (0,1) + MPGD (2..5). Silicon (6..9) is carried by
  // the unfolder but excluded from event finding.
  static constexpr size_t kNumTriggerDet = 6;

  std::vector<std::string> m_trk_collection_names = {
      "TOFBarrelTimeAlignRecHits",          "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",         "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits", "ForwardMPGDEndcapTimeAlignRecHits",
      "SiBarrelVertexTimeAlignRecHits",     "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits",    "B0TrackerTimeAlignRecHits"};

  VariadicPodioInput<edm4eic::TrackerHit, true> m_trk_in{this, m_trk_collection_names};
  PodioInput<edm4hep::MCParticle, true> m_mc_in{this, "MCParticles"};

  // Calorimeter clusters (time-aligned, frame level) — consumed ONLY to
  // compute the stored calorimeter-coincidence energy (weights[24]): the
  // summed cluster energy within the candidate's own combined window dt of
  // t0. Never part of the trigger decision, so efficiency is untouched.
  // Motivation (measured offline, ddis gold-coating, 100 frames): a real
  // collision deposits GeV-scale synchronous calo energy, a background
  // splash deposits ~MeV — windowed E_sum alone separates real from fake at
  // rank AUC 0.914 where every count-based statistic saturates at ~0.77.
  // See the eventbuilder README, "Outlook / TODO".
  // 10 entries, deliberately: JOmniFactory distributes a flat wiring-tag list
  // EVENLY across variadic inputs, and m_trk_in above has 10 — an uneven pair
  // throws "can't be distributed evenly" at wiring time. EcalLumiSpec is the
  // anticipated-but-not-yet-wired calo (see TimeAlignment_factory.h); as an
  // optional input it resolves to nullptr until someone wires it, at which
  // point it starts contributing automatically.
  std::vector<std::string> m_calo_collection_names = {
      "B0ECalTimeAlignClusters", "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters", "EcalEndcapPTimeAlignClusters",
      "HcalBarrelTimeAlignClusters", "HcalEndcapPInsertTimeAlignClusters",
      "LFHCALTimeAlignClusters", "EcalFarForwardZDCTimeAlignClusters",
      "HcalFarForwardZDCTimeAlignClusters", "EcalLumiSpecTimeAlignClusters"};
  VariadicPodioInput<edm4eic::Cluster, true> m_clu_in{this, m_calo_collection_names};

  PodioOutput<edm4hep::EventHeader> m_candidates_out{this, "EventCandidates"};

  void Configure() {
    auto* pm = GetApplication()->GetJParameterManager();
    pm->SetDefaultParameter("eventbuilder:topology:threshold", m_threshold,
        "minimum time-coincident hits in one (theta,phi) bin to fire the trigger; "
        "0 or 1 = 'no topology tightening' and clamp to the physical minimum (2), "
        "since a coincidence and a weighted t0 need at least a pair of hits");
    pm->SetDefaultParameter("eventbuilder:topology:theta_bins", m_theta_bins,
        "topology grid: number of polar-angle bins over [0, pi]");
    pm->SetDefaultParameter("eventbuilder:topology:phi_bins", m_phi_bins,
        "topology grid: number of azimuth bins over [0, 2pi)");
    pm->SetDefaultParameter("eventbuilder:topology:threshold_tof", m_threshold_tof,
        "minimum time-coincident TOF-group hits (within dt_tof, see Process()) in one "
        "(theta,phi) bin; same semantics and clamp-to-2 as topology:threshold");
    pm->SetDefaultParameter("eventbuilder:trigger:t0_estimator", m_t0_estimator,
        "t0 estimator: 'mean' = weighted mean of all window hits (default; "
        "late-tail biased ~+1 ns, smallest breakdown tail); 'leading' = "
        "coincidence-confirmed leading-edge anchored mean (bias-free median)");
    pm->SetDefaultParameter("eventbuilder:findable:pt_min", m_acc_pt_min,
        "in-acceptance expected-track count (weights[27]): minimum MC pT (GeV/c). "
        "Default 0.1 = EIC Yellow Report tracking requirement ('100 MeV pi'). "
        "Truth bookkeeping only — never a trigger cut");
    pm->SetDefaultParameter("eventbuilder:findable:eta_abs", m_acc_eta_abs,
        "in-acceptance expected-track count (weights[27]): maximum MC |eta|. "
        "Default 3.5 = EIC Yellow Report tracking coverage. Truth bookkeeping only");
    pm->SetDefaultParameter("eventbuilder:statistics:null", m_stat_null,
        "null model for the STORED significance weights[3]/[23] (trigger decision "
        "unchanged; granularity is part of the name, <law>_<granularity>): "
        "'poisson_cell' = exact Poisson tail, per-cell LOO mu ladder (legacy; assumes "
        "independent arrivals); 'nbinom_band' = negative-binomial tail moment-fitted "
        "to the frame's per-band sidebands (honest dispersion); 'nbinom_cell' = NB at "
        "the per-cell LOO mean with band-measured Fano dispersion (hybrid); "
        "'empirical_band' = rank against the per-band sideband windowed-count "
        "distribution (assumption-free, self-calibrating). Bare 'poisson'/'nbinom'/"
        "'empirical' alias to poisson_cell/nbinom_band/empirical_band. "
        "weights[28] (empirical S) and weights[29] (band Fano) are stored in every mode");
    pm->SetDefaultParameter("eventbuilder:topology:mu_min_hits", m_mu_min_hits,
        "minimum frame hits in a cell (or band) for its own mu estimate to be used; "
        "below this the ladder falls back to the band (then global) estimate");
    pm->SetDefaultParameter("eventbuilder:topology:mu_loo", m_mu_loo,
        "1 = leave-one-out scoring null (stored S/p_tail evaluated against mu computed "
        "without the tested window's own hits, removing self-masking circularity); "
        "0 = legacy frame-wide null. Trigger thresholds are frame-wide either way");
    // Tracklet-pointing windows (weights[25], score-only). Sized to absorb
    // solenoid curvature + multiple scattering, NOT detector resolution;
    // defaults tuned offline on sim quantities (barrel AUC 0.864, +endcap
    // 0.911) — retune on reco quantities pending.
    pm->SetDefaultParameter("eventbuilder:tracklet:dphi_barrel", m_trk_dphi_barrel,
        "tracklet: max |dphi| (rad) for a barrel fast-hit pair");
    pm->SetDefaultParameter("eventbuilder:tracklet:dphi_endcap", m_trk_dphi_endcap,
        "tracklet: max |dphi| (rad) for a same-side endcap fast-hit pair");
    pm->SetDefaultParameter("eventbuilder:tracklet:dr_min", m_trk_dr_min,
        "tracklet: min radial lever arm (mm) between barrel layers");
    pm->SetDefaultParameter("eventbuilder:tracklet:dz_min", m_trk_dz_min,
        "tracklet: min |dz| lever arm (mm) between endcap planes");
    pm->SetDefaultParameter("eventbuilder:tracklet:z0_max", m_trk_z0_max,
        "tracklet: max |z0| (mm) of the barrel pair's beamline extrapolation");
    pm->SetDefaultParameter("eventbuilder:tracklet:r_tol", m_trk_r_tol,
        "tracklet: projective r tolerance (mm) for endcap pairs, |r2 - r1*z2/z1|");
    // Trigger-tag bitmask thresholds (weights[26], inclusion tags, score-only).
    pm->SetDefaultParameter("eventbuilder:tags:zone_nhits", m_tag_zone_nhits,
        "tag bits 0/2/4: min cluster nhits for an ECal zone region tag");
    pm->SetDefaultParameter("eventbuilder:tags:b0_min_hits", m_tag_b0_min_hits,
        "tag bit 6: min in-window B0 hits for the far-forward tag");
    pm->SetDefaultParameter("eventbuilder:tags:zdc_min_hits", m_tag_zdc_min_hits,
        "tag bit 7: min summed in-window ZDC-ECal cluster hits");
    // STAGE 1 trigger TRK/CAL thresholds (weights[31]/[32], score-only —
    // see the member declarations above for why these exist unapplied).
    pm->SetDefaultParameter("eventbuilder:trigger:min_tracklets", m_trig_min_tracklets,
        "TRK term of offline TIME && (TRK || CAL): min n_tracklets (weights[25]) "
        "to count as tracker content. Not applied to phys_flag/emission.");
    pm->SetDefaultParameter("eventbuilder:trigger:min_cal_energy", m_trig_min_cal_energy,
        "CAL term of offline TIME && (TRK || CAL): min E_calo (GeV, weights[24]) "
        "to count as calo content. 0 = off (UNMEASURED PLACEHOLDER). Not applied "
        "to phys_flag/emission.");
    pm->SetDefaultParameter("eventbuilder:trigger:cal_nsigma", m_cal_nsigma,
        "N-sigma half-width for the per-cluster gate feeding e_calo (weights[24]): "
        "gate = N x cluster.getTimeError(). Mirrors CalCoincidence_factory::nsigma.");
    pm->SetDefaultParameter("eventbuilder:trigger:cal_late_nsigma_asym", m_cal_late_nsigma_asym,
        "extra LATE acceptance for the e_calo per-cluster gate, in units of "
        "that cluster's own sigma, added to the + side only. Mirrors "
        "CalCoincidence_factory::late_nsigma_window_asym. UNMEASURED PLACEHOLDER "
        "-- retune via an efficiency/purity scan.");
    pm->SetDefaultParameter("eventbuilder:trigger:cal_min_energy_barrel", m_cal_min_energy_barrel,
        "per-cluster energy floor (GeV) for e_calo, EcalBarrel/HcalBarrel clusters "
        "-- YR near-barrel (|eta|<1) minimum detectable-photon energy.");
    pm->SetDefaultParameter("eventbuilder:trigger:cal_min_energy_forward", m_cal_min_energy_forward,
        "per-cluster energy floor (GeV) for e_calo, every other calo collection "
        "-- YR forward/backward (1<=|eta|<4) minimum detectable-photon energy.");
    pm->SetDefaultParameter("eventbuilder:trigger:veto", m_trig_veto,
        "1 = ACTUALLY apply TIME && (TRK || CAL): a candidate failing both "
        "min_tracklets and min_cal_energy is dropped before emission (ACTS/"
        "CaloIsland never run on it). DESTRUCTIVE for real candidates -- not "
        "reversible without reprocessing at veto=0. 0 (default) = score-only, "
        "matches weights[33]/[34] without affecting output.");
    // A coincidence requires >= 2 hits (fewer cannot form a cluster, and t0
    // over 0 hits is undefined). Treat threshold 0/1 as the sentinel "no
    // topology cut / loosest": clamp to the minimum pair. This keeps behaviour
    // and the self-describing stored value (weights[12]) consistent.
    if (m_threshold < 2)
      m_threshold = 2;
    if (m_threshold_tof < 2)
      m_threshold_tof = 2;
  }

  void ChangeRun(int32_t /*run_nr*/) {}

  // One time-aligned trigger hit, flattened across detectors.
  struct FastHit {
    double time;
    float sigma;
    bool is_tof; // detector index 0/1 (TOF barrel/endcap) vs MPGD (2-5)
    Int_t theta1, phi1, theta2, phi2;
    // Cylindrical position + barrel/endcap split, consumed only by the
    // tracklet-pointing score at emit (weights[25]) — never by the trigger.
    // Barrel layers = TOFBarrel(0), MPGDBarrel(2), OuterMPGDBarrel(3);
    // endcap disks = TOFEndcap(1), Backward(4)/Forward(5) MPGDEndcap.
    float r_cyl, z_pos, phi_raw;
    bool is_barrel;
  };

  // Incremental per-cell windowed-count tracker (one instance per (grid,
  // cell) -- see tof1/tof2 in Process(), today's only instantiation).
  // Maintains the best count of hits achievable in any `window`-wide span
  // among the hits currently active in THIS cell (i.e. within the outer
  // scan's [lo, hi)), updated on Insert()/EvictFront() as hits enter/leave,
  // instead of rescanned from scratch on every query the way the O(n)
  // tofSubwindowMax used to be. Reading `best` is O(1); a caller's combined
  // max (e.g. cmax_tof in Process()) is then the max of `best` over all
  // cells, O(ncells) per query -- the same order the outer scan's own
  // per-iteration argmax already costs. Generic in the hit pool it tracks
  // (today: the TOF-only pool), so a future second detector group with its
  // own window width can reuse this unchanged by instantiating a second set
  // of trackers over its own pool.
  //
  // Correctness argument:
  //  - Insert(idx) treats idx as a candidate WINDOW END (backward-looking):
  //    advance back_lo (the earliest active hit still within `window` of
  //    idx) forward as needed, then the window ending at idx has count
  //    active.size()-back_lo. Recording the max of this over every insertion
  //    enumerates the same set of candidate windows the original forward-
  //    anchored scan did (try every start, grow forward) -- just walked in
  //    the opposite direction -- so both agree on the true max.
  //  - EvictFront() removes the globally-earliest active hit in this cell,
  //    mirroring the outer scan's lo advancing. If it wasn't part of the
  //    window that currently achieves `best` (best_pos > 0 before the
  //    shift), that window is still intact among what remains, so `best` is
  //    kept (indices shifted down by one to stay valid positions into the
  //    now-shorter deque). If it WAS part of that window (best_pos == 0),
  //    `best` may no longer be achievable, and Recompute() reruns the same
  //    two-pointer scan tofSubwindowMax always used, now bounded by this
  //    cell's own (typically small, since background spreads across cells)
  //    active size rather than the whole window's hit count.
  //
  //  Amortized cost over one frame: O(n) for Insert (each hit entered once)
  //  plus the sum of every Recompute()'s cost. In the common case
  //  (background roughly spread across cells) that sum stays small, since
  //  each cell's own active size is a fraction of n. In the pathological
  //  case where essentially all tracked hits land in one cell, Recompute()
  //  can fire on most evictions with a shrinking-but-still-large active size
  //  each time, degrading toward the same O(n) per query the incremental
  //  fix (tof_idx-only version) already gives -- i.e. this is a real
  //  improvement in the typical case, not an unconditional asymptotic
  //  guarantee independent of how background is distributed across cells.
  struct CellTracker {
    std::deque<size_t> active;   // global hit indices in this cell, time-ordered
    size_t back_lo   = 0;        // position in `active`: earliest hit still
                                  // within `window` of the most recent insertion
    int best         = 0;        // best count-in-window among all insertions so far
    size_t best_pos  = 0;        // position in `active` where that window starts

    void Insert(size_t idx, const std::vector<FastHit>& hits, float window) {
      active.push_back(idx);
      while (hits[idx].time - hits[active[back_lo]].time > double(window))
        ++back_lo;
      const int count = int(active.size() - back_lo);
      if (count > best) {
        best     = count;
        best_pos = back_lo;
      }
    }

    void EvictFront(const std::vector<FastHit>& hits, float window) {
      active.pop_front();
      if (back_lo > 0)
        --back_lo;
      if (best_pos == 0)
        Recompute(hits, window);
      else
        --best_pos;
    }

    void Recompute(const std::vector<FastHit>& hits, float window) {
      best = 0;
      back_lo = 0;
      best_pos = 0;
      if (active.empty())
        return;
      size_t i = 0, j = 0;
      while (i < active.size()) {
        while (j < active.size() &&
               hits[active[j]].time - hits[active[i]].time <= double(window))
          ++j;
        const int count = int(j - i);
        if (count > best) {
          best     = count;
          best_pos = i;
        }
        ++i;
        if (j < i)
          j = i;
      }
    }
  };

  void Process(int64_t run_number, uint64_t event_number) {
    // Topology grid dimensions (runtime; two staggered grids of ntheta x nphi).
    const int ntheta = std::max(1, m_theta_bins);
    const int nphi   = std::max(1, m_phi_bins);
    const int ncells = ntheta * nphi;

    // -- 1. Flatten the fast (TOF+MPGD) hits, with their (theta,phi) bins ------
    std::vector<FastHit> hits;
    for (size_t det = 0; det < kNumTriggerDet && det < m_trk_in().size(); ++det) {
      const auto* coll = m_trk_in().at(det);
      if (coll == nullptr)
        continue;
      for (const auto& hit : *coll) {
        const float sigma = hit.getTimeError();
        if (!(sigma > 0.0f))
          continue; // unusable (unset/zero) resolution -- would blow up the 1/sigma^2 weight
        if (min_edep() > 0.0f && hit.getEdep() < min_edep())
          continue; // soft background hit -- excluded from event finding only
        const double hx = hit.getPosition()[0], hy = hit.getPosition()[1],
                     hz = hit.getPosition()[2];
        FastHit fh{hit.getTime(), sigma, det < 2, 999, 999, 999, 999,
                   float(std::sqrt(hx * hx + hy * hy)), float(hz),
                   float(std::atan2(hy, hx)),
                   det == 0 || det == 2 || det == 3};
        thetaPhiBinCalc(hit, ntheta, nphi, fh.theta1, fh.phi1, fh.theta2, fh.phi2);
        hits.push_back(fh);
      }
    }
    std::sort(hits.begin(), hits.end(),
              [](const FastHit& a, const FastHit& b) { return a.time < b.time; });

    // -- 2. Per-frame resolution-derived coincidence windows --------------------
    // dt/dt_tof are NOT fixed constants: each is the widest pairwise time gap
    // two hits of THIS FRAME's own measured resolution could genuinely differ
    // by at nsigma_window confidence. For a pair with resolutions sigma_a,
    // sigma_b the matched (quadrature) width is sqrt(sigma_a^2+sigma_b^2),
    // which is maximised -- over every pair the pool could produce -- when
    // both hits sit at the pool's own worst resolution; mixing in any
    // better-resolution hit only narrows it. So max(sigma) then sqrt(2)* is
    // exact, not a conservative stand-in: it IS the worst-case pairwise
    // width, not an approximation of it. max_sigma_all is dominated by MPGD
    // (~10 ns) whenever any MPGD hit is present in the frame, which
    // background alone guarantees in practice; max_sigma_tof is the TOF-only
    // subset's own (~30 ps) resolution. Measured directly off this frame's
    // hits (getTimeError()), not a hardcoded reference constant, so this adds
    // no parameter beyond nsigma_window itself.
    float max_sigma_all = 0.0f, max_sigma_tof = 0.0f;
    size_t n_tof = 0;
    for (const auto& h : hits) {
      max_sigma_all = std::max(max_sigma_all, h.sigma);
      if (h.is_tof) {
        ++n_tof;
        max_sigma_tof = std::max(max_sigma_tof, h.sigma);
      }
    }
    const float dt     = float(m_nsigma_window) * std::sqrt(2.0f) * max_sigma_all;
    const float dt_tof = float(m_nsigma_window) * std::sqrt(2.0f) * max_sigma_tof;

    // -- MC collision times (for phys/fake labelling) ---------------------------
    // Signal-stream particles keep their raw generator status; the merger
    // offsets background streams by +1000*stream (2001, 3001, ...). Status 1
    // (generator-level final state) exists in EVERY generator — unlike the
    // pythia-only beam-remnant status 61 used before, which silently labelled
    // ALL rapgap/lAger/DVCS candidates fake. Every status-1 particle of one
    // collision carries that collision's frame time, so dedupe within the
    // (frame-derived) coincidence window to get one entry per injected signal
    // collision.
    std::vector<double> mc_times;
    std::vector<int> mc_nexp; // per-collision in-acceptance stable charged
                              // (the weights[27] expected-track count)
    if (m_mc_in() != nullptr) {
      const double dedup = double(dt);
      for (const auto& mcp : *m_mc_in()) {
        if (mcp.getGeneratorStatus() != 1)
          continue;
        const double t = mcp.getTime();
        int ci = -1;
        for (size_t j = 0; j < mc_times.size(); ++j)
          if (std::abs(mc_times[j] - t) < dedup) {
            ci = int(j);
            break;
          }
        if (ci < 0) {
          mc_times.push_back(t);
          mc_nexp.push_back(0);
          ci = int(mc_times.size()) - 1;
        }
        if (std::abs(mcp.getCharge()) > 0.01) {
          const auto p    = mcp.getMomentum();
          const double pt = std::hypot(double(p.x), double(p.y));
          if (pt > m_acc_pt_min &&
              std::abs(std::asinh(double(p.z) / pt)) < m_acc_eta_abs)
            ++mc_nexp[ci];
        }
      }
    }

    // -- Pileup collision times (cross-class SIGNAL_FREQ mixing) ---------------
    // Same idea as mc_times above, but for the OTHER real collisions
    // SignalBackgroundMerger cross-mixed in under PILEUP_FREQ (status offsets
    // 7000+, one +1000 block per source class — see README's "Signal pileup"
    // status table; status%1000==1 is that stream's final-state tag, same
    // convention the +1000*stream comment above already documents for
    // background). These are genuine additional ep collisions, not detector
    // noise, so a candidate matching one of these is a different thing from a
    // candidate matching nothing at all (see flag=3 below) — lumping them
    // into "fake" would hide whether the eventbuilder is actually capable of
    // recovering pileup collisions as their own events, not just correctly
    // rejecting them as "not the tagged one." Deduped by (stream, time), not
    // time alone: two different source classes' pileup collisions landing
    // close together in time must NOT merge into one entry.
    std::vector<double> pileup_times;
    std::vector<int> pileup_nexp;
    std::vector<int> pileup_stream;
    if (m_mc_in() != nullptr) {
      const double dedup = double(dt);
      for (const auto& mcp : *m_mc_in()) {
        const int status = mcp.getGeneratorStatus();
        if (status < 7000 || status % 1000 != 1)
          continue; // not a pileup-signal final-state particle
        const int stream = status / 1000;
        const double t = mcp.getTime();
        int ci = -1;
        for (size_t j = 0; j < pileup_times.size(); ++j)
          if (pileup_stream[j] == stream && std::abs(pileup_times[j] - t) < dedup) {
            ci = int(j);
            break;
          }
        if (ci < 0) {
          pileup_times.push_back(t);
          pileup_nexp.push_back(0);
          pileup_stream.push_back(stream);
          ci = int(pileup_times.size()) - 1;
        }
        if (std::abs(mcp.getCharge()) > 0.01) {
          const auto p    = mcp.getMomentum();
          const double pt = std::hypot(double(p.x), double(p.y));
          if (pt > m_acc_pt_min &&
              std::abs(std::asinh(double(p.z) / pt)) < m_acc_eta_abs)
            ++pileup_nexp[ci];
        }
      }
    }

    // Calorimeter clusters of this frame, sorted by time — consumed only at
    // emit for the windowed energy sum (weights[24]) and the trigger-tag
    // bitmask (weights[26]). `coll` indexes m_calo_collection_names (1 =
    // barrel ECal, 2 = backward ECal, 3 = forward ECal, 7 = ZDC ECal — the
    // zones the bitmask reads); nhits = cluster hit count, the regional
    // hit-multiplicity the bitmask thresholds on.
    struct CaloCl {
      double t, e;
      float sigma; // clu.getTimeError() -- the calo-specific per-cluster
                   // resolution e_calo's gate is built from below, NOT dt
                   // (which is TOF/MPGD hit resolution -- see max_sigma_all).
      int coll, nhits;
      bool operator<(const CaloCl& o) const { return t < o.t; }
    };
    std::vector<CaloCl> calo_te;
    for (size_t cd = 0; cd < m_clu_in().size(); ++cd) {
      const auto* cin = m_clu_in().at(cd);
      if (cin == nullptr)
        continue; // optional input: detector variant without this calo
      for (const auto& clu : *cin)
        calo_te.push_back({double(clu.getTime()), double(clu.getEnergy()),
                           clu.getTimeError(), int(cd), int(clu.getNhits())});
    }
    std::sort(calo_te.begin(), calo_te.end()); // kept for the debug_calo dump; the
                                                // per-cluster gate below no longer
                                                // needs sorted input (was: binary search)
    if (debug_calo()) {
      // Opt-in diagnostic (eventbuilder:trigger:debug_calo=1): dump the
      // frame's time-aligned clusters so cluster-time pathologies (e.g. a
      // calo whose cluster time is unset and lands at -r/c after alignment)
      // are visible without persisting Timeslice-level collections.
      jout << "EB debug_calo frame=" << event_number << " nclusters=" << calo_te.size() << jendl;
      for (const auto& c : calo_te)
        if (c.e > 0.1)
          jout << "  coll=" << c.coll << " t=" << c.t << " E=" << c.e
               << " nhits=" << c.nhits << jendl;
    }

    // B0 hit times (sorted) — trk input index 9, not part of event finding;
    // consumed only by the bitmask's far-forward tag (bit 6).
    std::vector<double> b0_times;
    if (m_trk_in().size() > 9 && m_trk_in().at(9) != nullptr)
      for (const auto& hit : *m_trk_in().at(9))
        b0_times.push_back(double(hit.getTime()));
    std::sort(b0_times.begin(), b0_times.end());

    // -- 3. O(N) sliding-window coincidence + (theta,phi) topology trigger -----
    // Per-frame background density estimate (used by the adaptive threshold and by the
    // stored per-candidate significance in both modes). Signal contaminates
    // mu by <~3% (one cluster among ~dT/dt window positions) — ignored.
    double mu = 0.0, mu_tof = 0.0;
    double dT = double(dt);
    if (hits.size() >= 2) {
      dT     = std::max(double(hits.back().time - hits.front().time), double(dt));
      mu     = (double(hits.size()) * double(dt) / dT) / double(ncells);
      mu_tof = (double(n_tof) * double(dt_tof) / dT) / double(ncells);
    }

    // INHOMOGENEOUS Poisson background model. The equal-solid-angle grid is
    // kept, but equal solid angle does NOT imply equal expected background
    // (the measured background is beam-axis-peaked, ~19x across theta rows):
    // the Poisson MEAN is allowed to vary with detector position. Each cell
    // of each staggered grid carries its own local null
    //   H_{0,i}: X_i ~ Poisson(mu_i),   mu_i = N_i * dt / dT
    // (N_i = this frame's hits in cell i, normalized to the coincidence
    // window) and, in adaptive mode, its own threshold k_i.
    //
    // ESTIMATOR INDEPENDENCE: these frame-wide mu_i drive the trigger
    // thresholds only, where the tested window's self-contamination is
    // O(w/N) ~ 2% and irrelevant to the firing decision. The STORED
    // significance of each candidate is evaluated against a LEAVE-ONE-OUT
    // estimate instead — the same ladder recomputed with the candidate's own
    // window excluded (counts minus the window's contribution, time span
    // dT - dt) — so no candidate is ever scored against a null it
    // contributed to. A cross-frame running background (via
    // TimesliceBuffer_service) is the natural upgrade if per-frame
    // statistics ever become limiting; see the README Outlook.
    //
    // Statistics ladder (topology:mu_min_hits):
    //   N_i           >= mu_min_hits  -> per-cell mu_i
    //   N_band(theta) >= mu_min_hits  -> cos-theta-band mu/nphi
    //   otherwise                     -> global mu
    // At 1x1 every path degenerates to the global mu by construction.
    std::vector<double> mu1(ncells, mu), mu2(ncells, mu);
    std::vector<int> ncell1(ncells, 0), ncell2(ncells, 0);
    std::vector<int> nband1(ntheta, 0), nband2(ntheta, 0);
    const double mu_norm = (dT > 0) ? double(dt) / dT : 0.0;
    if (m_mu_mode > 0 && ncells > 1 && hits.size() >= 2) {
      for (const auto& h : hits) {
        ++ncell1[h.theta1 * nphi + h.phi1];
        ++ncell2[h.theta2 * nphi + h.phi2];
        ++nband1[h.theta1];
        ++nband2[h.theta2];
      }
      for (int i = 0; i < ncells; ++i) {
        const int th = i / nphi;
        auto pick = [&](int ni, int nband) -> double {
          if (m_mu_mode >= 2 && ni >= m_mu_min_hits)
            return double(ni) * mu_norm;                       // per-cell
          if (nband >= m_mu_min_hits)
            return double(nband) * mu_norm / double(nphi);     // cos-theta band
          return mu;                                           // global fallback
        };
        mu1[i] = pick(ncell1[i], nband1[th]);
        mu2[i] = pick(ncell2[i], nband2[th]);
      }
    }

    // Sideband/dispersion statistics (eventbuilder:statistics), hoisted here
    // (from what used to be "step 3c", after the scan) so the adaptive
    // k-thresholds just below can use the SAME per-band dispersion that
    // weights[28]/[29] are stored from later, instead of a second estimate.
    // Windowed-count distribution of the frame itself: disjoint dt-wide
    // windows tiling [t_first, t_last], counted per grid-1 cell, histogrammed
    // per cos-theta band WITH empty windows (they carry the bulk of the
    // probability mass at mu << 1). Correlated hit families (conversions,
    // showers) land in the sideband tail exactly as they land in candidate
    // windows, so this prices them correctly where a plain Poisson tail
    // cannot. Self-contamination: a candidate's own hits occupy O(1) of the
    // nphi x nwin windows of its band — no LOO needed at this sample size.
    // Always computed (not gated on adaptive mode): weights[28]/[29] are
    // stored every frame regardless, so this is not new cost added by
    // adaptive k-determination, only a reuse of cost already paid.
    // dt can be 0 when the frame has fewer than 2 fast hits (max_sigma_all
    // stays at its zero-initialized default, see step 2 above), which would
    // make dT/dt a 0/0 division. Guard it explicitly rather than relying on
    // std::max(1, ...) to mask int(NaN), which is undefined behavior.
    const int nwin = (dt > 0.0) ? std::max(1, int(std::ceil(dT / double(dt)))) : 1;
    std::vector<std::vector<long long>> band_hist(ntheta, std::vector<long long>(1, 0));
    std::vector<double> band_mean(ntheta, 0.0), band_var(ntheta, 0.0), band_fano(ntheta, 0.0);
    const long long band_N = (long long)(nphi) * nwin;
    if (!hits.empty()) {
      const double t_lo = double(hits.front().time);
      std::vector<uint32_t> wc(size_t(ncells) * size_t(nwin), 0);
      for (const auto& h : hits) {
        int wj = int((double(h.time) - t_lo) / double(dt));
        wj     = std::min(std::max(wj, 0), nwin - 1);
        ++wc[size_t(h.theta1 * nphi + h.phi1) * size_t(nwin) + size_t(wj)];
      }
      for (int th = 0; th < ntheta; ++th) {
        auto& H = band_hist[th];
        double s1 = 0.0, s2 = 0.0;
        for (int ph = 0; ph < nphi; ++ph)
          for (int j = 0; j < nwin; ++j) {
            const uint32_t c = wc[size_t(th * nphi + ph) * size_t(nwin) + size_t(j)];
            if (H.size() <= size_t(c))
              H.resize(size_t(c) + 1, 0);
            ++H[c];
            s1 += c;
            s2 += double(c) * double(c);
          }
        band_mean[th] = s1 / double(band_N);
        band_var[th]  = s2 / double(band_N) - band_mean[th] * band_mean[th];
        band_fano[th] = (band_mean[th] > 0.0) ? band_var[th] / band_mean[th] : 0.0;
      }
    }
    // Empirical add-one upper-tail p of a count c in a band.
    auto empiricalTail = [&](int c, int band) -> double {
      if (band < 0 || band >= ntheta || c <= 0)
        return 1.0;
      const auto& H = band_hist[band];
      long long ge = 0;
      for (size_t k = std::min(size_t(c), H.size()); k < H.size(); ++k)
        ge += H[k];
      return double(ge + 1) / double(band_N + 1);
    };

    // Per-cell applied thresholds k_i: in adaptive mode each cell's own
    // P(X_i >= k_i) < p_fake (floored at topology:threshold); in fixed mode
    // every cell uses the configured floor. k1/k2 (combined-group, the ones
    // that actually gate firing — see the scan below) use a NEGATIVE-BINOMIAL
    // tail: mean = the cell/band/global mu ladder (mu1[i]/mu2[i], unchanged),
    // variance = band_fano[band] * mean — i.e. the mean keeps whatever
    // resolution the ladder already earned, but the WIDTH is corrected by the
    // separately measured local overdispersion (Fano, real and measured, see
    // weights[29]). A plain Poisson k systematically UNDERSTATES the true
    // accidental rate wherever background is bursty (Fano > 1: conversion/
    // shower hit families, not independent arrivals) — nbinomMinK degenerates
    // exactly to poissonMinK wherever the local Fano is <= 1, so this is a
    // strict widening, never a tighter threshold than Poisson would give.
    // k/k_tof (the GLOBAL scalars) stay on plain Poisson: k is cosmetic only
    // (never gates firing — see weights[9]'s comment, the applied per-cell
    // k1[i]/k2[i] do that job); k_tof's pool is comparatively sparse and its
    // own dispersion hasn't been separately measured (v1 scope decision).
    const bool adaptive = adaptiveMode();
    int k = m_threshold;  // representative/global value (weights[9] stores the winning cell's k_i)
    int k_tof = std::numeric_limits<int>::max(); // TOF group fires only in adaptive mode
    std::vector<int> k1(ncells, m_threshold), k2(ncells, m_threshold);
    if (adaptive) {
      k     = std::max(m_threshold, poissonMinK(mu, double(p_fake())));
      k_tof = std::max(m_threshold_tof, poissonMinK(mu_tof, double(p_fake())));
      for (int i = 0; i < ncells; ++i) {
        const int th1 = i / nphi;
        k1[i] = std::max(m_threshold,
                         nbinomMinK(mu1[i], band_fano[th1] * mu1[i], double(p_fake())));
        // mu2/k2 index the SAME (theta,phi) cell id space as mu1/k1 (see
        // Process() step 3's mu1/mu2 loop, both sized ncells and indexed by
        // i = theta*nphi+phi) even though they're accumulated from grid-2
        // (staggered) hit assignments — band_fano is grid-1-only (no
        // separate grid-2 sideband measurement exists), so grid-2 cells
        // reuse grid-1's band Fano as a proxy. Same simplification already
        // implicit in the stored "nbinom_cell" significance at emission.
        k2[i] = std::max(m_threshold,
                         nbinomMinK(mu2[i], band_fano[th1] * mu2[i], double(p_fake())));
      }
    }

    // Candidates found by the scan, buffered so MC matching can be done
    // globally (nearest candidate per collision) once every candidate of the
    // frame is known — see the matching pass after the scan loop.
    struct Pending {
      double t0, t0sigma, S, cluster_size, S_tof;
      int cmax, cmax_tof, win_theta, win_phi, win_grid;
      double e_calo;
      int n_tracklets, tagbits;
      double mu_win;  // winning cell's LEAVE-ONE-OUT scoring null (weights[7]) —
                      // the mu that S/p_tail were computed against
      int k_win;      // winning cell's applied trigger threshold (weights[9], frame-wide)
    };
    std::vector<Pending> pending;

    std::vector<Int_t> occupancy1(ncells, 0);
    std::vector<Int_t> occupancy2(ncells, 0);
    // Per-cell incremental windowed-count trackers (see CellTracker above),
    // instantiated here for the TOF-only pool -- one per cell, per grid.
    // Kept in lockstep with occupancy1/occupancy2 below: Insert() when a TOF
    // hit enters via hi, EvictFront() when one leaves via lo.
    std::vector<CellTracker> tof1(ncells), tof2(ncells);
    const auto cell = [nphi](Int_t theta, Int_t phi) { return theta * nphi + phi; };
    size_t lo = 0, hi = 0;
    while (lo < hits.size()) {
      while (hi < hits.size() && (hits[hi].time - hits[lo].time) <= dt) {
        ++occupancy1[cell(hits[hi].theta1, hits[hi].phi1)];
        ++occupancy2[cell(hits[hi].theta2, hits[hi].phi2)];
        if (hits[hi].is_tof) {
          tof1[cell(hits[hi].theta1, hits[hi].phi1)].Insert(hi, hits, dt_tof);
          tof2[cell(hits[hi].theta2, hits[hi].phi2)].Insert(hi, hits, dt_tof);
        }
        ++hi;
      }

      // Track the ARGMAX cell too (grid, theta bin, phi bin of the winning
      // cell) — stored per candidate so the offline diagnostics can build
      // the empirical occupancy map of where triggers actually fire, the
      // ground truth for how non-uniform the per-cell background share is
      // (the flat mu assumes 1/ncells everywhere; see weights[20..22]).
      // Winning cell = the MOST SIGNIFICANT one under its own local null
      // H_{0,i} (min P(X >= occ_i | mu_i)), not the max raw count: with
      // non-uniform mu_i a hot beam-axis cell needs more hits than a quiet
      // one to be equally surprising. The trigger fires when any cell
      // exceeds ITS OWN k_i. With mu_mode=0 (or 1x1) mu_i == mu everywhere
      // and this reduces exactly to the legacy max-count behavior.
      int cmax = 0, win_cell = 0, win_grid = 1;
      bool fired = false;
      double best_logp = 1.0; // log10 p is <= 0; 1.0 = "nothing yet"
      for (int i = 0; i < ncells; ++i) {
        const int c1 = occupancy1[i], c2 = occupancy2[i];
        if (c1 >= k1[i]) {
          const double lp = poissonLogSurvival10(c1, mu1[i]);
          if (!fired || lp < best_logp) { best_logp = lp; cmax = c1; win_cell = i; win_grid = 1; }
          fired = true;
        }
        if (c2 >= k2[i]) {
          const double lp = poissonLogSurvival10(c2, mu2[i]);
          if (!fired || lp < best_logp) { best_logp = lp; cmax = c2; win_cell = i; win_grid = 2; }
          fired = true;
        }
      }
      if (!fired) {
        // No cell over threshold: keep the max-count cell as the diagnostic
        // "winning" record (same as legacy) without firing.
        for (int i = 0; i < ncells; ++i) {
          if (int(occupancy1[i]) > cmax) { cmax = occupancy1[i]; win_cell = i; win_grid = 1; }
          if (int(occupancy2[i]) > cmax) { cmax = occupancy2[i]; win_cell = i; win_grid = 2; }
        }
      }
      const int win_theta = win_cell / nphi, win_phi = win_cell % nphi;
      const int k_win     = (win_grid == 1) ? k1[win_cell] : k2[win_cell];
      // cmax_tof: max over every cell's already-maintained `best` (O(1)
      // each), not a fresh O(n_tof) rescan -- see CellTracker.
      int cmax_tof = 0;
      for (int i = 0; i < ncells; ++i) {
        cmax_tof = std::max(cmax_tof, tof1[i].best);
        cmax_tof = std::max(cmax_tof, tof2[i].best);
      }
      const bool trigger = fired || (cmax_tof >= k_tof);

      if (trigger) {
        double sumw = 0.0, sumwt = 0.0;
        if (m_t0_estimator == "leading") {
          // Leading edge: sigma_min identifies the most precise detector class
          // in the window; the anchor is the EARLIEST hit of that class (the
          // beta~1 population arrives first — the late tail is one-sided).
          // Average only hits whose time is compatible with the anchor within
          // their own resolution; stragglers (slow hadrons, curlers) are
          // excluded instead of dragging t0 late at full 1/sigma^2 weight.
          double sigma_min = std::numeric_limits<double>::infinity();
          for (size_t j = lo; j < hi; ++j)
            sigma_min = std::min(sigma_min, double(hits[j].sigma));
          // COINCIDENCE-CONFIRMED anchor: the earliest precise hit that has
          // at least one OTHER precise hit within 3*(sigma_i+sigma_j). A
          // lone early accidental cannot hijack the anchor (the naive
          // min-anchor killed the +1 ns bias but raised the >10 ns
          // breakdown fraction from 32% to 51% in the 100-frame A/B).
          double t_lead = std::numeric_limits<double>::infinity();
          for (size_t j = lo; j < hi; ++j) {
            const double sj = double(hits[j].sigma);
            if (sj > 3.0 * sigma_min || double(hits[j].time) >= t_lead)
              continue;
            for (size_t k = lo; k < hi; ++k) {
              if (k == j)
                continue;
              const double sk = double(hits[k].sigma);
              if (sk > 3.0 * sigma_min)
                continue;
              if (std::abs(double(hits[k].time) - double(hits[j].time))
                  <= 3.0 * (sj + sk)) {
                t_lead = double(hits[j].time);
                break;
              }
            }
          }
          if (std::isfinite(t_lead)) {
            for (size_t j = lo; j < hi; ++j) {
              const double sig = double(hits[j].sigma);
              if (double(hits[j].time) - t_lead > 3.0 * (sig + sigma_min))
                continue;
              const double w = 1.0 / (sig * sig);
              sumw  += w;
              sumwt += w * hits[j].time;
            }
          }  // no confirmed anchor -> sums stay 0 -> legacy-mean fallback below
        }
        if (sumw <= 0.0) {  // t0_mode 0, or degenerate leading-edge selection
          for (size_t j = lo; j < hi; ++j) {
            const double w = 1.0 / (double(hits[j].sigma) * double(hits[j].sigma));
            sumw  += w;
            sumwt += w * hits[j].time;
          }
        }
        const double t0      = sumwt / sumw;
        const double sigmat0 = std::sqrt(1.0 / sumw);
        const double t0sigma     = double(nsigma_window()) * sigmat0;

        // S/S_tof: -10*log10(p_tail), a dB-scale significance built directly
        // from the EXACT one-sided Poisson tail, valid at any mu. Replaces
        // the earlier Pearson residual (c-mu)/sqrt(mu), which only reads as
        // a Gaussian sigma asymptotically at large mu — our operating mu is
        // routinely <1, well outside that regime. Deliberately NOT a
        // z-score/sigma count: it is a monotone ranking/rarity score on a
        // log scale, always well-defined and never reliant on the Gaussian
        // approximation. Computed via poissonLogSurvival10 (log-space upper
        // tail), NOT -log10(poissonSurvival(...)): the 1-cdf form saturates
        // at double precision (p ~ 2e-16 -> S caps at ~157, and exact zeros
        // clamp to the 1e-300 floor -> S piles up at 3000), which destroyed
        // the ranking of every strong candidate.
        // Local null, LEAVE-ONE-OUT: the stored significance is evaluated
        // against the winning cell's mu_i recomputed WITHOUT this window's
        // own hits (counts minus the window's contribution, off-window time
        // span dT - dt) — the tested window never contributes to the null it
        // is tested against. The frame-wide mu_win/k_win (trigger side)
        // remain what fired; only the scoring null is decontaminated.
        const int w_hits = int(hi - lo);
        const double dT_loo = std::max(dT - double(dt), double(dt));
        const double norm_loo = double(dt) / dT_loo;
        double mu_loo;
        if (m_mu_loo == 0) {
          // legacy: score against the frame-wide winning-cell null
          mu_loo = (win_grid == 1) ? mu1[win_cell] : mu2[win_cell];
        } else {
          const int N_loo = std::max(int(hits.size()) - w_hits, 0);
          const double mu_glob_loo = double(N_loo) * norm_loo / double(ncells);
          if (m_mu_mode > 0 && ncells > 1) {
            // Window contribution to the winning cell is exactly cmax; its
            // contribution to the band is bounded by w_hits (subtracting
            // w_hits is conservative for the band estimate).
            const auto& nc = (win_grid == 1) ? ncell1 : ncell2;
            const auto& nb = (win_grid == 1) ? nband1 : nband2;
            const int ni_loo = std::max(nc[win_cell] - cmax, 0);
            const int nb_loo = std::max(nb[win_theta] - w_hits, 0);
            if (m_mu_mode >= 2 && nc[win_cell] >= m_mu_min_hits)
              mu_loo = double(ni_loo) * norm_loo;
            else if (nb[win_theta] >= m_mu_min_hits)
              mu_loo = double(nb_loo) * norm_loo / double(nphi);
            else
              mu_loo = mu_glob_loo;
          } else {
            mu_loo = mu_glob_loo;
          }
          if (!(mu_loo > 0.0))
            mu_loo = mu_glob_loo; // empty off-window cell: fall back to global LOO
        }
        const double S = -10.0 * poissonLogSurvival10(cmax, mu_loo);
        const double S_tof = -10.0 * poissonLogSurvival10(cmax_tof, mu_tof);
        // Calorimeter-coincidence energy (weights[24]) and the per-zone
        // cluster tags feeding the bitmask: one binary-search range over the
        // time-sorted cluster list. Zone tag = any in-window cluster of that
        // zone's ECal with nhits >= m_tag_zone_nhits (a >= 10-hit region IS a
        // shower, so the cluster is the regional hit count already
        // aggregated); ZDC tag thresholds the summed in-window ZDC-ECal
        // cluster hit count instead (its showers fragment across clusters).
        double e_calo = 0.0;
        bool zone_bwd = false, zone_bar = false, zone_fwd = false;
        int zdc_nhits = 0;
        // Per-cluster gate: N-sigma x THIS cluster's own getTimeError() (+
        // late-side asym allowance), same formula as CalCoincidence_factory,
        // applied inline (see m_cal_nsigma's declaration comment for why).
        // Energy floor is the YR minimum detectable-photon energy, split
        // barrel (coll==1 EcalBarrel, coll==4 HcalBarrel) vs forward/
        // backward -- same indices the zone tags below already use.
        for (const auto& c : calo_te) {
          if (!(c.sigma > 0.0f))
            continue; // unusable (unset/zero) resolution
          const double e_min = (c.coll == 1 || c.coll == 4) ? m_cal_min_energy_barrel
                                                              : m_cal_min_energy_forward;
          if (c.e < e_min)
            continue;
          const double half_win = m_cal_nsigma * double(c.sigma);
          const double dtc = c.t - t0;
          if (dtc < -half_win || dtc > half_win + m_cal_late_nsigma_asym * double(c.sigma))
            continue;
          e_calo += c.e;
          if (c.nhits >= m_tag_zone_nhits) {
            zone_bar |= (c.coll == 1);
            zone_bwd |= (c.coll == 2);
            zone_fwd |= (c.coll == 3);
          }
          if (c.coll == 7)
            zdc_nhits += c.nhits;
        }
        // Vertex-pointing tracklets (weights[25]) from the gated fast hits of
        // THIS window [lo,hi) — score only, never in the trigger decision.
        // Barrel pairs: straight-line z0 extrapolation to the beamline.
        // Endcap pairs (same side): projective consistency r2 ~ r1*z2/z1
        // (the short two-plane lever arm makes z0 itself ill-conditioned).
        // Window sizes absorb solenoid curvature (~0.05-0.1 rad at 0.5 GeV),
        // NOT detector resolution; values validated offline (AUC 0.911).
        int trk_barrel = 0, trk_bwd = 0, trk_fwd = 0;
        for (size_t a = lo; a < hi; ++a) {
          for (size_t b = a + 1; b < hi; ++b) {
            const FastHit &h1 = hits[a], &h2 = hits[b];
            double dphi = std::abs(h1.phi_raw - h2.phi_raw);
            if (dphi > M_PI)
              dphi = 2.0 * M_PI - dphi;
            if (h1.is_barrel && h2.is_barrel) {
              const double dr = double(h2.r_cyl) - double(h1.r_cyl);
              if (std::abs(dr) < m_trk_dr_min || dphi >= m_trk_dphi_barrel)
                continue;
              const FastHit& in = (dr > 0) ? h1 : h2;
              const FastHit& out = (dr > 0) ? h2 : h1;
              const double z0 = double(in.z_pos) -
                                double(in.r_cyl) * (double(out.z_pos) - double(in.z_pos)) /
                                    std::abs(dr);
              if (std::abs(z0) < m_trk_z0_max)
                ++trk_barrel;
            } else if (!h1.is_barrel && !h2.is_barrel &&
                       (h1.z_pos > 0) == (h2.z_pos > 0)) {
              if (std::abs(double(h2.z_pos) - double(h1.z_pos)) < m_trk_dz_min ||
                  dphi >= m_trk_dphi_endcap)
                continue;
              const FastHit& near_ = (std::abs(h1.z_pos) < std::abs(h2.z_pos)) ? h1 : h2;
              const FastHit& far_ = (std::abs(h1.z_pos) < std::abs(h2.z_pos)) ? h2 : h1;
              const double r_exp =
                  double(near_.r_cyl) * std::abs(double(far_.z_pos) / double(near_.z_pos));
              if (std::abs(double(far_.r_cyl) - r_exp) < m_trk_r_tol)
                ((far_.z_pos < 0) ? trk_bwd : trk_fwd)++;
            }
          }
        }
        const int n_tracklets = trk_barrel + trk_bwd + trk_fwd;
        // In-window B0 hit count (bit 6's ingredient).
        const auto b0_lo = std::lower_bound(b0_times.begin(), b0_times.end(), t0 - double(dt));
        const auto b0_hi = std::upper_bound(b0_times.begin(), b0_times.end(), t0 + double(dt));
        const int n_b0 = int(b0_hi - b0_lo);
        // Trigger-tag bitmask (weights[26]) — INCLUSION tags, each an
        // independent positive excess; the builder cuts on none of them.
        //   bit 0: backward-ECal region tag —  in-window backward-ECal cluster, nhits >= 10
        //   bit 1: bit 0 with track match   —  ... AND >= 1 backward endcap tracklet
        //   bit 2: barrel-ECal region tag
        //   bit 3: bit 2 with track match (barrel tracklet)
        //   bit 4: forward-ECal region tag
        //   bit 5: bit 4 with track match (forward endcap tracklet)
        //   bit 6: far-forward B0 tag       —  >= 4 in-window B0 hits
        //   bit 7: ZDC tag                  —  in-window ZDC-ECal cluster hits sum >= 50
        // Combined tags (evaluated OFFLINE from these bits, never stored):
        //   cTag1 = (#set{0,2,4} >= 2) && bit6      cTag2 = (#set{1,3,5} >= 2) && bit6
        //   cTag3 = (#set{0,2,4} >= 2) && bit7      cTag4 = (#set{1,3,5} >= 2) && bit7
        //   cTag5 = (#set{0,2,4} == 3)              cTag6 = (#set{1,3,5} >= 2)
        int tagbits = 0;
        if (zone_bwd)                  tagbits |= 1 << 0;
        if (zone_bwd && trk_bwd > 0)   tagbits |= 1 << 1;
        if (zone_bar)                  tagbits |= 1 << 2;
        if (zone_bar && trk_barrel > 0) tagbits |= 1 << 3;
        if (zone_fwd)                  tagbits |= 1 << 4;
        if (zone_fwd && trk_fwd > 0)   tagbits |= 1 << 5;
        if (n_b0 >= m_tag_b0_min_hits)        tagbits |= 1 << 6;
        if (zdc_nhits >= m_tag_zdc_min_hits)  tagbits |= 1 << 7;

        // Defer emission: MC matching is done GLOBALLY after the scan (each
        // collision goes to its NEAREST candidate), so no flag is decided
        // here. Greedy in-scan matching mislabeled real events whenever an
        // earlier background window also sat within +-dt of the collision
        // (observed: a S=0.7 background window 25 ns early stole the match
        // from the collision's own S=2.9 window 17 ns away).
        pending.push_back({t0, t0sigma, S, double(hi - lo), S_tof,
                           cmax, cmax_tof, win_theta, win_phi, win_grid,
                           e_calo, n_tracklets, tagbits, mu_loo, k_win});

        // consume the cluster: clear counts and jump past it
        for (size_t j = lo; j < hi; ++j) {
          --occupancy1[cell(hits[j].theta1, hits[j].phi1)];
          --occupancy2[cell(hits[j].theta2, hits[j].phi2)];
          if (hits[j].is_tof) {
            tof1[cell(hits[j].theta1, hits[j].phi1)].EvictFront(hits, dt_tof);
            tof2[cell(hits[j].theta2, hits[j].phi2)].EvictFront(hits, dt_tof);
          }
        }
        lo = hi;
      } else {
        // slide: drop the earliest hit and advance
        --occupancy1[cell(hits[lo].theta1, hits[lo].phi1)];
        --occupancy2[cell(hits[lo].theta2, hits[lo].phi2)];
        if (hits[lo].is_tof) {
          tof1[cell(hits[lo].theta1, hits[lo].phi1)].EvictFront(hits, dt_tof);
          tof2[cell(hits[lo].theta2, hits[lo].phi2)].EvictFront(hits, dt_tof);
        }
        ++lo;
        if (hi < lo)
          hi = lo;
      }
    }

    // Sideband/dispersion statistics (band_hist, band_mean, band_var,
    // band_fano, empiricalTail) were computed earlier in this function, right
    // before "Per-cell applied thresholds" — moved there in step 3 so the
    // adaptive-mode k-thresholds below can use the SAME per-band dispersion
    // that weights[28]/[29] are stored from, instead of a second estimate.

    // -- 4. Global MC matching (nearest candidate), then emission --------------
    // Each injected collision is assigned to the candidate with the SMALLEST
    // |t0 - t_mc| among those within +-dt, skipping candidates already
    // matched to another collision. Candidates left unmatched are fakes.
    std::vector<int> flags(pending.size(), 2);
    std::vector<double> mcts(pending.size(), kNoMCTime);
    std::vector<int> nexps(pending.size(), 0); // weights[27] per candidate
    for (size_t m = 0; m < mc_times.size(); ++m) {
      const double mct = mc_times[m];
      int best = -1;
      double best_d = double(dt);
      for (size_t i = 0; i < pending.size(); ++i) {
        if (flags[i] == 1)
          continue; // one collision per candidate; next-nearest takes over
        const double dd = std::abs(pending[i].t0 - mct);
        if (dd < best_d) {
          best_d = dd;
          best   = int(i);
        }
      }
      if (best >= 0) {
        flags[best] = 1;
        mcts[best]  = mct;
        nexps[best] = mc_nexp[m];
      }
    }
    // Second pass: match still-unmatched candidates (flags[i]==2) against the
    // pileup collisions the same way — nearest candidate, exclusive
    // assignment, one pileup collision per candidate. flag=3 means "found a
    // real collision, just not the tagged primary one" — distinct from
    // flag=2 "matched nothing" (pure accidental). Candidates already claimed
    // by the primary (flags[i]==1) are never reconsidered here.
    for (size_t m = 0; m < pileup_times.size(); ++m) {
      const double mct = pileup_times[m];
      int best = -1;
      double best_d = double(dt);
      for (size_t i = 0; i < pending.size(); ++i) {
        if (flags[i] != 2)
          continue;
        const double dd = std::abs(pending[i].t0 - mct);
        if (dd < best_d) {
          best_d = dd;
          best   = int(i);
        }
      }
      if (best >= 0) {
        flags[best] = 3;
        mcts[best]  = mct;
        nexps[best] = pileup_nexp[m];
      }
    }
    const int n_pileup_inj = int(pileup_times.size());
    // Per-candidate true-overlap count (weights[35]): how many REAL collisions
    // (primary mc_times + cross-class pileup_times) fall within +-dt of THIS
    // candidate's own t0 -- independent of the EXCLUSIVE one-to-one matching
    // above. flags[i]/mcts[i] only ever record ONE collision per candidate
    // (nearest wins, "next-nearest takes over" for a DIFFERENT candidate), so
    // two collisions genuinely landing in the same window still get spread
    // across flags/mcts as if they were separate events -- this count is what
    // actually answers "did more than one real collision coincide here,"
    // i.e. whether the eventbuilder is seeing genuine combined pile-up (one
    // candidate, several collisions) as opposed to merely resolving several
    // nearby-but-separable collisions into their own candidates. ==1 is the
    // normal single-collision case (matches what flags/mcts already show);
    // ==0 is a pure-noise candidate (flags[i]==2); >=2 is the case flags/mcts
    // cannot represent today.
    std::vector<int> n_mc_coincident(pending.size(), 0);
    // WHICH pileup class(es) coincided (weights[36]): bit `pileup_stream`
    // set whenever that stream's collision falls in THIS candidate's window.
    // Streams are status/1000 (7..29 today, see vendor/sro/README.md's
    // "Signal pileup" table for the stream->class name mapping -- that
    // mapping is a convention of the merger tooling, not known to eicrecon,
    // so decoding a bit back to e.g. "dvcs" has to happen offline against
    // that table). The primary/tagged class isn't a "stream" (native status
    // 1, no shift) and isn't given a bit here -- phys_flag/mc_t already say
    // whether the primary coincided; this field is only about the OTHER,
    // cross-mixed classes riding along in the same window.
    std::vector<uint32_t> coincident_streams(pending.size(), 0);
    for (size_t i = 0; i < pending.size(); ++i) {
      const double t0i = pending[i].t0;
      int cnt = 0;
      uint32_t smask = 0;
      for (double mct : mc_times)
        if (std::abs(t0i - mct) <= double(dt))
          ++cnt;
      for (size_t m = 0; m < pileup_times.size(); ++m)
        if (std::abs(t0i - pileup_times[m]) <= double(dt)) {
          ++cnt;
          if (pileup_stream[m] >= 0 && pileup_stream[m] < 32)
            smask |= (1u << pileup_stream[m]);
        }
      n_mc_coincident[i]   = cnt;
      coincident_streams[i] = smask;
    }
    for (size_t i = 0; i < pending.size(); ++i) {
      const Pending& c = pending[i];
      const int band     = std::min(std::max(c.win_theta, 0), ntheta - 1);
      const double p_emp = empiricalTail(c.cmax, band);
      const double s_emp = -10.0 * std::log10(p_emp);
      const double fano  = band_fano[band];
      // Mode-selected stored significance/tail; trigger decision already made.
      // <law>_<granularity> grammar; bare names alias to their canonical
      // granularity (poisson->cell, nbinom->band, empirical->band).
      double s_stored = c.S;
      // c.S was computed via poissonLogSurvival10 (log-space, non-saturating)
      // against the exact same (cmax, mu_win) pair used here. Deriving
      // p_stored from it (p = 10^(-S/10)) keeps the two fields consistent by
      // construction; calling poissonSurvival(c.cmax, c.mu_win) directly
      // would saturate around 2e-16 (S capped near 157) and disagree with
      // c.S for any candidate significant beyond that floor.
      double p_stored = std::pow(10.0, -s_stored / 10.0);
      if (m_stat_null == "empirical" || m_stat_null == "empirical_band") {
        s_stored = s_emp;
        p_stored = p_emp;
      } else if (m_stat_null == "nbinom" || m_stat_null == "nbinom_band") {
        p_stored = nbinomSurvival(c.cmax, band_mean[band], band_var[band]);
        s_stored = -10.0 * std::log10(std::max(p_stored, 1e-300));
      } else if (m_stat_null == "nbinom_cell") {
        // Hybrid: per-cell LOO mean, band-measured dispersion. Falls back to
        // the Poisson tail automatically wherever Fano_band <= 1.
        p_stored = nbinomSurvival(c.cmax, c.mu_win, fano * c.mu_win);
        s_stored = -10.0 * std::log10(std::max(p_stored, 1e-300));
      }
      // Early veto (eventbuilder:trigger:veto, default off — see the member
      // declaration above): TIME already kept this window (it's in
      // `pending`), but if neither TRK nor CAL clears its threshold, drop
      // it here, before emitCandidate() ever creates the EventCandidates
      // row. flags[i]/mcts[i]/nexps[i] were already assigned by the MC-
      // matching pass above and are simply discarded for this i — a
      // flags[i]==1 candidate vetoed here is a real collision genuinely
      // lost, not recoverable offline (see the member comment).
      if (m_trig_veto && c.n_tracklets < m_trig_min_tracklets &&
          c.e_calo < m_trig_min_cal_energy) {
        continue;
      }
      emitCandidate(run_number, c.t0, c.t0sigma, flags[i], s_stored, c.cluster_size,
                    mcts[i], c.S_tof, c.mu_win, mu_tof, c.k_win, k_tof, dt, dt_tof,
                    c.cmax, c.cmax_tof, c.win_theta, c.win_phi, c.win_grid,
                    c.e_calo, c.n_tracklets, c.tagbits, nexps[i],
                    p_stored, s_emp, fano, n_pileup_inj, n_mc_coincident[i],
                    coincident_streams[i]);
    }
  }

  // ---------------------------------------------------------------------------
  // Sentinel for "no matched MC collision" in the candidate encoding (ns).
  static constexpr double kNoMCTime = -1.0e9;

  // Tracklet-pointing windows (weights[25]) and trigger-tag thresholds
  // (weights[26]) are runtime parameters — see Configure():
  //   eventbuilder:tracklet:{dphi_barrel,dphi_endcap,dr_min,dz_min,z0_max,r_tol}
  //   eventbuilder:tags:{zone_nhits,b0_min_hits,zdc_min_hits}
  // Defaults were tuned offline (barrel AUC 0.864, +endcap 0.911 on sim
  // quantities; retune on reco pending). All score-side, never the trigger.

  // (MC matching is done globally in Process() step 4 — nearest candidate
  // per collision — replacing the old greedy per-candidate matchPhysics,
  // which let an earlier background window steal a collision's label.)

  // Smallest k with P(X >= k | Poisson(mu)) < p — the adaptive-mode count threshold.
  static int poissonMinK(double mu, double p) {
    if (!(mu > 0.0))
      return 2;
    double pmf = std::exp(-mu), cdf = pmf;
    for (int k = 0; k < 100; ++k) {
      if (1.0 - cdf < p)
        return k + 1;
      pmf *= mu / double(k + 1);
      cdf += pmf;
    }
    return 100;
  }

  // Exact upper-tail probability P(X >= c | Poisson(mu)) — the HONEST
  // rarity of a c-hit clump under background, valid at any mu. Stored raw as
  // weights[23] so offline cuts have the untransformed tail probability.
  // The trigger's own firing decision already uses this tail via
  // poissonMinK. Returns 1.0 for c<=0, and clamps to a tiny floor so
  // log10(p) downstream never hits -inf. NOTE: as a probability this
  // saturates at double precision (1-cdf cancellation floors p at ~2e-16);
  // the stored significance S therefore uses poissonLogSurvival10 below,
  // which stays exact far past that.
  // Negative-binomial upper tail P(X >= c) moment-fitted to (mean m,
  // variance v) — the count law of Poisson-arriving hit FAMILIES
  // (overdispersed, Fano = v/m > 1). Degenerates to the Poisson tail when
  // the measured dispersion is consistent with independence (v <= m).
  static double nbinomSurvival(int c, double m, double v) {
    if (c <= 0)
      return 1.0;
    if (!(m > 0.0))
      return 0.0;
    if (v <= m * (1.0 + 1e-9))
      return poissonSurvival(c, m);
    const double q = 1.0 - m / v;     // NB 'success' parameter, in (0,1)
    const double r = m * m / (v - m); // dispersion (smaller = clumpier)
    double pmf = std::pow(1.0 - q, r), cdf = pmf; // P(X=0)
    for (int k = 0; k + 1 < c && k < 100000; ++k) {
      pmf *= q * (double(k) + r) / double(k + 1);
      cdf += pmf; // after the loop: cdf = P(X <= c-1)
    }
    return std::max(1.0 - cdf, 1e-300);
  }

  // Smallest k with P(X >= k | NegBinom(mean, var)) < p — the FAR-honest
  // analog of poissonMinK once background is measurably overdispersed
  // (Fano = var/mean > 1). Degenerates to poissonMinK exactly (same code
  // path, via nbinomSurvival's own v<=m guard) wherever the local dispersion
  // is consistent with independence, so this is a strict widening: never a
  // LOWER k than plain Poisson would give at the same (mean, p).
  static int nbinomMinK(double mean, double var, double p) {
    if (!(mean > 0.0))
      return 2;
    if (var <= mean * (1.0 + 1e-9))
      return poissonMinK(mean, p);
    const double q = 1.0 - mean / var;
    const double r = mean * mean / (var - mean);
    double pmf = std::pow(1.0 - q, r), cdf = pmf;
    for (int k = 0; k < 100; ++k) {
      if (1.0 - cdf < p)
        return k + 1;
      pmf *= q * (double(k) + r) / double(k + 1);
      cdf += pmf;
    }
    return 100;
  }

  static double poissonSurvival(int c, double mu) {
    if (c <= 0)
      return 1.0;
    if (!(mu > 0.0))
      return (c > 0) ? 0.0 : 1.0; // empty-background frame: any hit is "impossible" under H0
    double pmf = std::exp(-mu), cdf = pmf; // P(X=0)
    for (int k = 1; k < c; ++k) {
      pmf *= mu / double(k);
      cdf += pmf;
    }
    return std::max(1.0 - cdf, 1e-300);
  }

  // log10 P(X >= c | Poisson(mu)), computed entirely in log space — the
  // basis of the stored significance S = -10*log10(p_tail). Unlike
  // poissonSurvival's 1-cdf (which cancels catastrophically once the tail
  // drops below ~2e-16, capping S at ~157 and piling saturated candidates
  // into one bin), this sums the upper-tail series directly from k=c:
  //   P(X>=c) = pmf(c) * [1 + mu/(c+1) + mu^2/((c+1)(c+2)) + ...]
  // with log pmf(c) = -mu + c*ln(mu) - lgamma(c+1). The bracket converges
  // geometrically (ratio <= mu/(c+1) < 1 whenever the tail is small), and
  // the result is exact to double precision at ANY c — a c=100 clump at
  // mu=0.15 gets its true S of ~2400 instead of a meaningless 3000 clamp.
  // c<=0 -> 0.0 (p=1); mu<=0 with c>0 -> -300 (the same conventional floor
  // as poissonSurvival's clamp: empty background, maximally significant).
  static double poissonLogSurvival10(int c, double mu) {
    if (c <= 0)
      return 0.0;
    if (!(mu > 0.0))
      return -300.0;
    // For mu >= c the tail is O(1) — the direct 1-cdf has no cancellation
    // risk there, and the series below would converge slowly.
    if (mu >= double(c))
      return std::log10(poissonSurvival(c, mu));
    const double logpmf = -mu + double(c) * std::log(mu) - std::lgamma(double(c) + 1.0);
    double term = 1.0, series = 1.0;
    for (int k = 1; k < 1000; ++k) {
      term *= mu / double(c + k);
      series += term;
      if (term < 1e-17 * series)
        break;
    }
    return (logpmf + std::log(series)) / std::log(10.0);
  }

  // Candidate encoding (EventCandidates; the unfolder copies ALL entries
  // verbatim into the per-child EventBuilderInfo and appends the frame's
  // candidate count — keep eb_diagnostics.py in sync):
  //   timeStamp   = llround(t0_ns * 1000)  (t0 in ps, integer)
  //   weights[0]  = t0 (ns)          weights[1] = t0sigma (ns)
  //   weights[2]  = phys/fake flag (MC-only: 1 = matched the tagged primary,
  //                 2 = matched nothing (accidental), 3 = matched a pileup
  //                 collision -- a genuine other real collision, not the
  //                 tagged one; see pileup_times above)
  //   weights[3]  = S, combined-group (dt) significance = -10*log10(p_tail);
  //                 a dB-scale rarity score, NOT a Gaussian sigma count (see
  //                 poissonSurvival)
  //   weights[4]  = cluster size (fast hits in the consumed window)
  //   weights[5]  = matched MC collision time (ns; MC-only, kNoMCTime if none)
  //   weights[6]  = S_tof, TOF-group (dt_tof) significance
  //   -- self-describing frame/config block (files carry their own settings) --
  //   weights[7]  = mu       (measured combined-group background, hits/cell/window)
  //   weights[8]  = mu_tof   (measured TOF-group background)
  //   weights[9]  = k        (combined-group count threshold actually applied)
  //   weights[10] = k_tof    (TOF-group count threshold; huge when disabled)
  //   weights[11] = dt (ns) — the combined-group coincidence window THIS
  //                 CANDIDATE was found with; derived per frame from
  //                 nsigma_window and the frame's own measured resolution
  //                 (see Process()), NOT a fixed config constant, so this can
  //                 vary candidate to candidate if the frame's hit population
  //                 does. Read per-candidate offline, never assumed constant
  //                 across a file.
  //   weights[12] = topology:threshold
  //   weights[13] = topology:theta_bins   weights[14] = topology:phi_bins
  //   weights[15] = adaptive mode flag (derived: 1 iff p_fake in (0,1))
  //   weights[16] = p_fake
  //   weights[17] = dt_tof (ns) — same per-candidate caveat as weights[11]
  //   -- audit/instrumentation block (raw trigger counts + winning cell) --
  //   weights[18] = cmax      (combined-group coincidence count that fired —
  //                            the raw numerator of S, per the statistical
  //                            audit: needed offline for the (mu, c) plane)
  //   weights[19] = cmax_tof  (TOF-group coincidence count)
  //   weights[20] = winning cell theta bin   weights[21] = winning cell phi bin
  //   weights[22] = winning grid (1 = plain, 2 = half-cell staggered)
  //   weights[23] = p_tail = P(X >= cmax | Poisson(mu)) — the EXACT upper-tail
  //                 background probability of the combined-group clump; the
  //                 untransformed rarity that weights[3]'s S is now directly
  //                 derived from (S = -10*log10(p_tail)).
  //   weights[24] = E_calo (GeV) — calorimeter-coincidence energy: summed
  //                 time-aligned cluster energy with |t_cluster - t0| <= dt.
  //                 SCORE ONLY, never in the trigger decision (efficiency is
  //                 untouched); measured real-vs-fake rank AUC 0.981 at reco
  //                 cluster level where count statistics saturate at ~0.77.
  //   weights[25] = n_tracklets — vertex-pointing fast-hit pairs in the fired
  //                 window (barrel z0 + endcap projective, see Process()).
  //                 SCORE ONLY. Offline AUC 0.911 (barrel+endcap).
  //   weights[26] = trigger-tag bitmask (bits 0-7, see Process()): per-zone
  //                 ECal region tags (+track-matched variants), B0 and ZDC
  //                 far-forward tags. INCLUSION tags for downstream stream
  //                 assignment (union-of-tags); the builder cuts on none.
  //   weights[27] = n_expected — in-acceptance stable charged MC of the
  //                 MATCHED collision (pT > findable:pt_min, |eta| <
  //                 findable:eta_abs; YR defaults 0.1 GeV / 3.5). The
  //                 tracking-efficiency denominator, per candidate. 0 for
  //                 fakes (no matched collision) and MC-less running.
  //   weights[28] = S_empirical — significance of cmax against THIS frame's
  //                 per-band sideband windowed-count distribution
  //                 (eventbuilder:statistics; stored in every null mode for
  //                 offline A/B against the weights[3] of the active mode).
  //   weights[29] = Fano factor (variance/mean of sideband windowed counts)
  //                 of the winning cell's cos-theta band; 1 = Poisson-like,
  //                 > 1 = overdispersed (hit families).
  //   weights[30] = n_pileup_inj — number of pileup collisions (cross-class
  //                 SIGNAL_FREQ mixing, status 7000+) injected into this
  //                 candidate's frame; same value on every candidate of one
  //                 frame. Pairs with a count of flag==3 candidates (see
  //                 weights[2]) to give pileup efficiency/purity, the same
  //                 shape as the flag==1/N_inj pair already used for the
  //                 tagged primary.
  //   weights[31] = eventbuilder:trigger:min_tracklets — the TRK threshold a
  //                 run was configured with (min n_tracklets, weights[25],
  //                 to count as tracker content). Self-describing so offline
  //                 TIME && (TRK || CAL) evaluation reads the actual run
  //                 config instead of a second, possibly stale, hardcoded
  //                 copy. Score-only UNLESS eventbuilder:trigger:veto=1, in
  //                 which case candidates failing both this and weights[32]
  //                 are dropped before emission (never reach this row at
  //                 all, in fact — see weights[33]/[34] and Process()).
  //   weights[32] = eventbuilder:trigger:min_cal_energy — the CAL threshold
  //                 (GeV, min E_calo, weights[24]) a run was configured
  //                 with. Same veto caveat as weights[31]. 0 = off
  //                 (UNMEASURED PLACEHOLDER default).
  //   weights[33] = TRK pass/fail (1.0/0.0): n_tracklets (weights[25]) >=
  //                 weights[31], precomputed so a future veto/consumer
  //                 doesn't re-derive the comparison from raw fields.
  //   weights[34] = CAL pass/fail (1.0/0.0): E_calo (weights[24]) >=
  //                 weights[32]. Note weights[33]/[34] are computed and
  //                 stored REGARDLESS of veto — they describe every emitted
  //                 candidate's own TRK/CAL status, veto applied or not;
  //                 they are never both false on a veto=1 run's output,
  //                 since such candidates were dropped before this point.
  //   weights[35] = n_mc_coincident — count of REAL collisions (primary mc_times +
  //                 cross-class pileup_times) within +-dt of THIS candidate's
  //                 OWN t0, independent of the exclusive one-collision-per-
  //                 candidate matching that sets weights[2]/[5]/[27] above.
  //                 1 = the normal single-collision case; 0 = pure-noise
  //                 candidate (phys_flag==2); >=2 = genuine combined pile-up
  //                 (multiple real collisions landing in the SAME candidate
  //                 window) — a case phys_flag/mc_t cannot represent, since
  //                 that matching always assigns at most one collision per
  //                 candidate and pushes any other in-window collision onto
  //                 a different candidate if one exists nearby, or drops it
  //                 as unrecovered otherwise.
  //   weights[36] = coincident_streams — bitmask (bit N = 1<<N) of pileup
  //                 stream ids (status/1000) whose collision falls within
  //                 +-dt of THIS candidate's t0, i.e. WHICH cross-mixed
  //                 class(es) contributed to weights[35] beyond the primary.
  //                 The primary/tagged class has no stream id (native status
  //                 1, unshifted) and is never given a bit here -- phys_flag/
  //                 mc_t already say whether it coincided. Decoding a bit to
  //                 a class name (e.g. bit 14 = dvcs) requires the stream->
  //                 class table in vendor/sro/README.md's "Signal pileup"
  //                 section -- that mapping is a convention of the merger
  //                 tooling, not known to eicrecon itself. 0 = no cross-
  //                 mixed class coincided (n_mc_coincident counted only the
  //                 primary, or nothing at all).
  //                 weights[last] (appended by the unfolder) = candidates in
  //                 this frame.
  //  NOTE weights[3] and weights[23] carry the significance/tail of the
  //  ACTIVE eventbuilder:statistics:null mode (default poisson = legacy).
  void emitCandidate(int64_t run_number, double t0, double t0sigma, int phys_flag,
                     double signif, double cluster_size, double mc_t, double signif_tof,
                     double mu, double mu_tof, int k, int k_tof, float dt, float dt_tof,
                     int cmax, int cmax_tof, int win_theta, int win_phi, int win_grid,
                     double e_calo, int n_tracklets, int tagbits, int n_expected,
                     double p_stored, double s_emp, double fano, int n_pileup_inj,
                     int n_mc_coincident, uint32_t coincident_streams) {
    auto cand = m_candidates_out()->create();
    cand.setRunNumber(static_cast<uint32_t>(run_number));
    cand.setEventNumber(m_candidates_out()->size() - 1);
    cand.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    cand.setWeight(t0sigma);                                                 // ns (legacy)
    cand.addToWeights(t0);                                               // [0]
    cand.addToWeights(t0sigma);                                              // [1]
    cand.addToWeights(static_cast<double>(phys_flag));                   // [2]
    cand.addToWeights(signif);                                           // [3]
    cand.addToWeights(cluster_size);                                     // [4]
    cand.addToWeights(mc_t);                                             // [5]
    cand.addToWeights(signif_tof);                                       // [6]
    cand.addToWeights(mu);                                               // [7]
    cand.addToWeights(mu_tof);                                           // [8]
    cand.addToWeights(double(k));                                        // [9]
    cand.addToWeights(double(std::min(k_tof, 999)));                     // [10]
    cand.addToWeights(double(dt));                                       // [11]
    cand.addToWeights(double(m_threshold));                              // [12]
    cand.addToWeights(double(m_theta_bins));                             // [13]
    cand.addToWeights(double(m_phi_bins));                               // [14]
    cand.addToWeights(adaptiveMode() ? 1.0 : 0.0);                       // [15]
    cand.addToWeights(double(p_fake()));                                 // [16]
    cand.addToWeights(double(dt_tof));                                   // [17]
    cand.addToWeights(double(cmax));                                     // [18]
    cand.addToWeights(double(cmax_tof));                                 // [19]
    cand.addToWeights(double(win_theta));                                // [20]
    cand.addToWeights(double(win_phi));                                  // [21]
    cand.addToWeights(double(win_grid));                                 // [22]
    cand.addToWeights(p_stored);                                         // [23]
    cand.addToWeights(e_calo);                                           // [24]
    cand.addToWeights(double(n_tracklets));                              // [25]
    cand.addToWeights(double(tagbits));                                  // [26]
    cand.addToWeights(double(n_expected));                               // [27]
    cand.addToWeights(s_emp);                                            // [28]
    cand.addToWeights(fano);                                             // [29]
    cand.addToWeights(double(n_pileup_inj));                             // [30]
    cand.addToWeights(double(m_trig_min_tracklets));                     // [31]
    cand.addToWeights(m_trig_min_cal_energy);                            // [32]
    // Precomputed pass/fail, not just the raw threshold above: a future
    // online early-veto reads THIS, not weights[24]/[25] vs [31]/[32] —
    // one already-evaluated boolean instead of re-deriving the comparison
    // in a second code path. Still score-only today (unused by phys_flag).
    cand.addToWeights(n_tracklets >= m_trig_min_tracklets ? 1.0 : 0.0);   // [33]
    cand.addToWeights(e_calo >= m_trig_min_cal_energy ? 1.0 : 0.0);       // [34]
    cand.addToWeights(double(n_mc_coincident));                               // [35]
    cand.addToWeights(double(coincident_streams));                            // [36]
  }

  inline void thetaPhiBinCalc(edm4eic::TrackerHit hit, int ntheta, int nphi,
                              Int_t& thetaID1, Int_t& phiID1,
                              Int_t& thetaID2, Int_t& phiID2) {
    Double_t hitX = hit.getPosition()[0];
    Double_t hitY = hit.getPosition()[1];
    Double_t hitZ = hit.getPosition()[2];
    Double_t hitR = TMath::Sqrt(hitX * hitX + hitY * hitY + hitZ * hitZ);
    Double_t hitTheta = (hitR > 0.) ? TMath::ACos(hitZ / hitR) : 0.;
    if (hitTheta > TMath::Pi())
      hitTheta = TMath::Pi();
    if (hitTheta < 0.)
      hitTheta = 0.;

    Double_t hitPhi = TMath::ATan2(hitY, hitX);
    if (hitPhi < 0)
      hitPhi += 2 * TMath::Pi();

    // Grid 1: plain binning. Grid 2: shifted by half a cell in both angles
    // (the "stagger") so a cluster on a grid-1 edge is central in grid 2.
    //
    // theta is binned in COS(theta), not raw theta: equal-width cos(theta)
    // bins subtend equal solid angle (dOmega = dcos(theta) dphi), whereas
    // equal-width theta bins do not — the poles subtend far less solid angle
    // than the equator, so with raw-theta bins an isotropic background piles
    // up in the equatorial rows and the flat mu = N/ncells share (which
    // assumes every cell gets 1/ncells) systematically under-estimates the
    // local rate there, inflating those cells' S. Binning in cos(theta)
    // restores the equal-share assumption for an isotropic background. u runs
    // 0..ntheta over cos(theta) in [+1,-1] (theta 0..pi), so u increases with
    // theta, preserving the bin ordering.
    const Double_t u    = (1.0 - TMath::Cos(hitTheta)) * 0.5 * double(ntheta); // [0,ntheta]
    const Double_t dphi = 2. * TMath::Pi() / double(nphi);
    thetaID1 = int(u);
    thetaID2 = int(u + 0.5);
    phiID1   = hitPhi / dphi;
    phiID2   = (hitPhi + dphi / 2.) / dphi;

    // clamp/wrap into the grid (staggered bins can overflow by one; theta
    // clamps at the poles, phi wraps around)
    thetaID1 = std::min(std::max(thetaID1, 0), ntheta - 1);
    thetaID2 = std::min(std::max(thetaID2, 0), ntheta - 1);
    phiID1   = std::min(std::max(phiID1, 0), nphi - 1);
    phiID2   = std::min(std::max(phiID2 % nphi, 0), nphi - 1);
  }
};
