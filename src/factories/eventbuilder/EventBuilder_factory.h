// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// EventBuilder_factory
// -------------------
// Streaming "event building" stage. Runs once per time frame (JANA Timeslice
// level) over the time-aligned tracker hits and finds candidate physics
// events. For each candidate it computes a precise event time t0 and its
// half-width t0sigma.
//
//   t0      = inverse-variance weighted mean of the coincident
//             trigger-detector hit times: t0 = Sum(t_i/sigma_i^2) / Sum(1/sigma_i^2)
//   t0sigma = nsigma_window * sigma_t0, where sigma_t0 = 1/sqrt(Sum 1/sigma_i^2)
//
// nsigma_window also sets the coincidence windows (dt, dt_tof; see Process())
// that find a candidate in the first place. Each window is nsigma_window *
// sqrt(2) times the worst per-hit resolution present in that frame's pool,
// not a fixed constant.
//
// Per-hit sigma_i is edm4eic::TrackerHit::getTimeError(). Each detector's own
// TrackerHitReconstruction_factory (or equivalent) fills this from its
// digitization config. This factory reads that value directly and keeps no
// copy of its own, so it always matches what the detector plugin configured.
//
// TOF (about 30 ps resolution) dominates the weighting, so t0 is effectively
// set by the fastest detector. Silicon (about 2 microseconds resolution) is
// carried by the detector but does not help find events: only TOF and MPGD
// drive the coincidence.
//
// This factory only finds and times candidates. It does no reconstruction
// and creates no child events; the downstream EventUnfolder turns each
// surviving candidate into a JANA PhysicsEvent.
//
// Output: an edm4hep::EventHeaderCollection named "EventCandidates", one
// entry per candidate, ordered in time, plus one "EventBuilderFrameInfo"
// entry per frame carrying the frame-constant configuration (see the
// layout comments above emitCandidate() below). In short:
//   timeStamp  = round(t0_ns * 1000), t0 in ps, as an integer.
//   weight     = t0sigma_ns.
//   weights[0] = t0, in ns.
//   weights[1] = t0sigma, in ns.
//   weights[2] = flag: the number of real MC collisions inside this
//                candidate's own dt window. 0 = fake/accidental, 1 = one
//                real collision (normal), 2+ = true pile-up. Real means
//                flag > 0.

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

  // One confidence-level knob for the whole trigger side. This is a
  // separate parameter from EventUnfolder's own
  // "eventbuilder:unfolder:nsigma_window", which gates PhysicsEvent content
  // against an already-fixed t0, a different statistical step from anything
  // below.
  //
  // nsigma_window does two jobs:
  //   1. t0sigma = nsigma_window * sigma_t0: the reported candidate precision.
  //   2. The finding-stage coincidence windows dt and dt_tof, derived per
  //      frame in Process() from this value and each pool's own measured
  //      resolution. Raising nsigma_window widens both windows together, in
  //      proportion; the topology:threshold floor still protects the
  //      combined path regardless of window width.
  float m_nsigma_window = 3.0f;
  ParameterRef<float> nsigma_window{this, "nsigma_window", m_nsigma_window,
                                    "N-sigma confidence level for t0sigma AND for deriving "
                                    "the finding-stage coincidence windows (dt = nsigma_window * "
                                    "sqrt(2) * max per-hit resolution in the pool)"};

  // Topology trigger: minimum coincident hits in one staggered (theta,phi)
  // bin needed to fire. Raise it to cut accidental (fake) candidates from
  // dense machine background: the Poisson probability of N random hits
  // sharing one bin falls steeply as N grows.
  //
  // Registered by hand below, in Configure(), as "eventbuilder:topology:*",
  // instead of through ParameterRef. ParameterRef would prefix this
  // factory's own tag ("eventbuilder:trigger:*"), giving a 4-segment name
  // like "eventbuilder:trigger:topology:threshold". "topology" is its own
  // component of the eventbuilder:<component>:<parameter> namespace, at the
  // same level as trigger/coincidence/unfolder/prefilter, even though it is
  // implemented inside this factory.
  //
  // Tracklet windows and tag thresholds (see Configure() for parameter names
  // and meanings) are score-side only; they never affect the trigger
  // decision. They are runtime parameters, not derived from geometry, for
  // these reasons:
  //   dr_min/dz_min    lever-arm floors that exclude same-layer pairs. The
  //                    true values are the layer spacings (MPGD barrel r
  //                    about 550/745 mm, TOF about 640 mm; endcap planes dz
  //                    about 125 mm). A future geometry service could derive
  //                    these.
  //   dphi_*           curvature windows: dphi ~ 0.3*B*dr / (2*pT_min), plus
  //                    multiple scattering. Derivable from the field map and
  //                    layer radii, but only given a choice of pT_min, which
  //                    is a physics choice, not a geometric one.
  //   z0_max           vertex-region half-length: a machine/optics
  //                    parameter (bunch length, IR design), not detector
  //                    geometry.
  //   r_tol            projective tolerance, dominated by multiple
  //                    scattering: depends on material budget and pT
  //                    spectrum, so it is tuned, not derived.
  //   zone_nhits/b0/zdc detector-response thresholds (shower sizes), from
  //                    the trigger-study reference. Retune per campaign.
  double m_trk_dphi_barrel = 0.05;
  double m_trk_dphi_endcap = 0.10;
  double m_trk_dr_min      = 50.0;
  double m_trk_dz_min      = 50.0;
  double m_trk_z0_max      = 100.0;
  double m_trk_r_tol       = 20.0;
  int m_tag_zone_nhits  = 10;
  int m_tag_b0_min_hits = 4;
  int m_tag_zdc_min_hits = 50;

  // STAGE 1 TRIGGER: offline TIME && (TRK || CAL) thresholds on the already
  // stored n_tracklets/e_calo (weights[10]/[9]). Score only, like
  // tags:*/tracklet:* above: not applied to phys_flag or to candidate
  // emission. These exist so the threshold a run was configured with
  // travels with the data (weights[12]/[14]), instead of being hardcoded a
  // second time in offline analysis.
  int    m_trig_min_tracklets  = 1;
  double m_trig_min_cal_energy = 0.0;  // GeV; 0 = off. Not yet tuned: thresholds
                                        // e_calo, a sum over clusters, and needs a
                                        // background scan to set. Distinct from
                                        // cal_min_energy_barrel/forward below, which
                                        // floor one cluster at a time before it enters
                                        // the sum.

  // Per-cluster calo gate feeding e_calo (see the calo_te loop in
  // Process()): the same formula as CalTimeCoincidence (N sigma times
  // cluster.getTimeError(), plus a late-side asymmetric allowance, plus a
  // per-cluster minimum energy). Applied inline here, instead of consuming
  // CalTimeCoincidence's output, because that factory depends on
  // EventCandidates, which this factory is still producing; see the
  // eventbuilder README for that ordering constraint.
  double m_cal_nsigma            = 3.0;   // mirrors CalTimeCoincidence::nsigma
  double m_cal_late_nsigma_asym  = 0.5;   // mirrors CalTimeCoincidence::late_nsigma_window_asym; not yet tuned
  double m_cal_min_energy_barrel = 0.050; // GeV; near-barrel (|eta|<1) minimum detectable-photon energy
  double m_cal_min_energy_forward = 0.030; // GeV; forward/backward (1<=|eta|<4) minimum detectable-photon energy
  // Off by default. When on, this actually applies TIME && (TRK || CAL): a
  // candidate failing both thresholds above is dropped before it becomes an
  // EventCandidates row, so EventUnfolder never unfolds it. This is
  // destructive for any candidate that would have matched a real collision
  // (flags[i]==1 in Process()): a vetoed real candidate is lost and is not
  // recoverable without reprocessing with veto off.
  bool m_trig_veto = false;

  int m_threshold  = 2;
  int m_theta_bins = 12;
  int m_phi_bins   = 8;

  // Per-cell background model: equal solid angle does not mean equal
  // expected background, so each cell gets its own mu_i, following a
  // cell -> band -> global statistics ladder gated by mu_min_hits (the
  // minimum frame hits before an estimate is trusted). See Process() step 3.
  //
  // t0 estimator ("mean" or "leading"):
  //   mean    — inverse-variance weighted mean of all window hits (default).
  //             Biased slightly late (about +1 ns) by the one-sided late
  //             tail, but has the smallest breakdown tail.
  //   leading — coincidence-confirmed leading edge: anchors on the earliest
  //             precise hit pair and averages only hits compatible with it.
  //             A close-to-bias-free median, with a larger breakdown tail.
  std::string m_t0_estimator = "mean";

  // In-acceptance definition for the expected-track count stored in
  // weights[11]: stable charged signal MC with pT > pt_min and |eta| <
  // eta_abs. Defaults follow the EIC Yellow Report tracking requirements
  // (|eta| <= 3.5, minimum pT "100 MeV pi"; arXiv:2103.05419). Truth-side
  // bookkeeping only; never a trigger cut.
  double m_acc_pt_min  = 0.1; // GeV/c
  double m_acc_eta_abs = 3.5;

  // eventbuilder:statistics:null: which null model the stored significance
  // (weights[3]) and tail probability (weights[23]) are evaluated against.
  // Named as <law>_<granularity>. The trigger decision itself (count
  // thresholds, always the cell -> band -> global mu ladder) is the same in
  // every mode.
  //   "poisson_cell"   (default) exact Poisson tail on the per-cell,
  //                    leave-one-out mu ladder. Assumes independent hit
  //                    arrivals, which measured overdispersion (Fano > 1,
  //                    from conversion/shower hit families) contradicts, so
  //                    its p-values run too small.
  //   "nbinom_band"    negative-binomial tail, moment-fitted to the frame's
  //                    sideband windowed counts (per-band mean and
  //                    variance): the count law of Poisson-arriving hit
  //                    families. Honest dispersion, band-level conditioning.
  //   "nbinom_cell"    hybrid: negative-binomial tail at the per-cell
  //                    leave-one-out mean, with the band-measured
  //                    dispersion. Cell conditioning and honest dispersion
  //                    together; falls back to poisson_cell wherever
  //                    Fano <= 1.
  //   "empirical_band" rank of cmax against the frame's own sideband
  //                    windowed-count distribution (disjoint dt-wide
  //                    windows, including empty ones; add-one
  //                    p = (#{count >= c} + 1)/(N + 1)). Makes no
  //                    distributional assumption and is self-calibrating;
  //                    its smallest possible p per frame is about
  //                    1/(nphi*nwindows + 1), so pool frames for smaller p.
  // Bare aliases: "poisson" = poisson_cell, "nbinom" = nbinom_band,
  // "empirical" = empirical_band.
  //
  // weights[23] (empirical S) and weights[24] (band Fano factor) are always
  // stored, regardless of mode, so every run carries enough data to compare
  // modes offline.
  std::string m_stat_null = "poisson_cell";

  // The trigger's background ladder granularity is fixed: cell -> band ->
  // global, falling back automatically through mu_min_hits.
  static constexpr int m_mu_mode = 2;
  int m_mu_min_hits = 50;
  // Leave-one-out scoring null (see Process()): the stored S and p_tail are
  // evaluated against mu recomputed without the tested window's own hits.
  // 1 = on (removes self-masking, so a window is never scored against a
  // background estimate it contributed to); 0 = legacy frame-wide null.
  int m_mu_loo = 1;

  // Optional soft-hit rejection before event finding: synchrotron and
  // beam-gas hits are mostly low-eDep. 0 (default) = off.
  bool m_debug_calo = false;
  ParameterRef<bool> debug_calo{
      this, "debug_calo", m_debug_calo,
      "1 = dump each frame's time-aligned calo clusters (>0.1 GeV) to the log — "
      "diagnostic for cluster-time pathologies; score/trigger behavior unchanged"};

  float m_min_edep = 0.0f;
  ParameterRef<float> min_edep{
      this, "min_edep", m_min_edep,
      "minimum hit eDep (GeV) for a hit to participate in event finding; 0 = off"};

  // True-pile-up truth separation (trigger_classes, weights[25+]). flag
  // (weights[2]) counts the real collisions in a candidate's window and
  // trigger_classes_mask (weights[17]) says which classes they were, but
  // both discard the per-collision (time, stream) pairs. This stores that
  // list: count at weights[25], (time, stream) pairs after. On by default
  // (cheap: a few floats per candidate) since the .trigger_meta.txt
  // sidecar and per-collision hit attribution both need it.
  bool m_store_coincident_list = true;
  ParameterRef<bool> store_coincident_list{
      this, "store_coincident_list", m_store_coincident_list,
      "1 (default) = append the per-candidate coincident-collision "
      "(time, stream) list as trigger_classes, weights[25+] (count at "
      "[25], pairs at [26+2k]/[27+2k]). stream = generatorStatus/1000 of "
      "the collision's base status: 0 = native/tagged primary, >= 10 = a "
      "physics class (class_index = stream-10). 0 = off, weights end "
      "at [24]."};

  // Adaptive-threshold trigger mode: instead of the fixed topology:threshold
  // count, derives the count threshold per frame from that frame's own
  // background density. The expected background count per (bin x window) is
  //   mu = (hits.size() x dt / frame_span) / ncells
  // where dt is the derived combined-group coincidence window (see
  // Process()) and ncells = theta_bins x phi_bins. The trigger fires when a
  // bin's count is statistically incompatible with Poisson(mu): c >= the
  // smallest k with P(X >= k | Poisson(mu)) < p_fake.
  //
  // p_fake is the target accidental probability per (bin x window
  // position). Holding it fixed keeps the fake rate constant across
  // changing background conditions, where a fixed hit count would be too
  // loose in dense frames and too tight in quiet ones. The threshold is
  // floored at topology:threshold, so quiet frames never trigger on a bare
  // hit pair.
  //
  // The candidate's significance, S = -10*log10(p_tail) with
  // p_tail = P(X >= cmax | Poisson(mu)), is a dB-scale rarity score built
  // from the exact Poisson tail (see poissonSurvival). It is not a Gaussian
  // sigma/z-score: the operating mu is routinely below 1, well outside the
  // range where a z-score is meaningful. It is computed in both modes and
  // stored in the candidate (weights[3]), so the operating point can also
  // be chosen offline by cutting on S.
  //
  // p_fake alone selects the mode: a value in (0,1) enables the adaptive
  // per-frame threshold targeting that accidental probability; 0 (the
  // default) or 1 or more means fixed topology:threshold.
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
  // window dt is sized for MPGD (sigma about 10 ns; see the derivation in
  // Process()), which wastes TOF's much better timing resolution: a TOF
  // pair inside its own tight dt_tof is far more significant than the same
  // pair judged against dt. In adaptive mode the trigger therefore fires on
  // either group: the combined group (all fast hits within dt), or the TOF
  // group (TOF hits only, within dt_tof), each tested against its own
  // Poisson background with the same p_fake. Silicon has no trigger group:
  // it cannot time-resolve events. Both group significances are stored per
  // candidate regardless of mode. Neither dt nor dt_tof is a separate
  // parameter: both derive per frame from nsigma_window (see Process()).
  //
  // Floor on the TOF-group count threshold, the TOF-group counterpart of
  // topology:threshold. Registered next to it in Configure(), as
  // "eventbuilder:topology:threshold_tof".
  int m_threshold_tof = 2;

  // Trigger detectors: TOF (0,1) plus MPGD (2..5). Silicon (6..9) is
  // carried by the unfolder but excluded from event finding.
  static constexpr size_t kNumTriggerDet = 6;

  std::vector<std::string> m_trk_collection_names = {
      "TOFBarrelTimeAlignRecHits",          "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",         "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits", "ForwardMPGDEndcapTimeAlignRecHits",
      "SiBarrelVertexTimeAlignRecHits",     "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits",    "B0TrackerTimeAlignRecHits"};

  VariadicPodioInput<edm4eic::TrackerHit, true> m_trk_in{this, m_trk_collection_names};
  PodioInput<edm4hep::MCParticle, true> m_mc_in{this, "MCParticles"};

  // Calorimeter clusters (time-aligned, frame level). Used only to compute
  // the stored calorimeter-coincidence energy (weights[9]): the summed
  // cluster energy within the candidate's own combined window dt of t0.
  // Never part of the trigger decision, so efficiency is unaffected. A real
  // collision deposits GeV-scale synchronous calo energy, while a
  // background splash deposits only about an MeV, so this energy sum
  // separates real from fake candidates better than any count-based
  // statistic; see the eventbuilder README, "Outlook / TODO".
  //
  // 10 entries, deliberately: JOmniFactory distributes a flat wiring-tag
  // list evenly across variadic inputs, and m_trk_in above also has 10.
  // EcalLumiSpec has no plugin yet; as an optional input it resolves to
  // nullptr until one is wired, at which point it starts contributing
  // automatically.
  std::vector<std::string> m_calo_collection_names = {
      "B0ECalTimeAlignClusters", "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters", "EcalEndcapPTimeAlignClusters",
      "HcalBarrelTimeAlignClusters", "HcalEndcapPInsertTimeAlignClusters",
      "LFHCALTimeAlignClusters", "EcalFarForwardZDCTimeAlignClusters",
      "HcalFarForwardZDCTimeAlignClusters", "EcalLumiSpecTimeAlignClusters"};
  VariadicPodioInput<edm4eic::Cluster, true> m_clu_in{this, m_calo_collection_names};

  PodioOutput<edm4hep::EventHeader> m_candidates_out{this, "EventCandidates"};
  // One entry per frame (not per candidate): the frame-constant
  // configuration hoisted out of the per-candidate weights[]. Layout is
  // documented above emitFrameInfo() below. EventUnfolder copies this
  // record into every child event with eventNumber = frame*1000, so
  // offline readers join it to candidates via eventNumber / 1000.
  PodioOutput<edm4hep::EventHeader> m_frame_info_out{this, "EventBuilderFrameInfo"};

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
        "in-acceptance expected-track count (weights[11]): minimum MC pT (GeV/c). "
        "Default 0.1 = EIC Yellow Report tracking requirement ('100 MeV pi'). "
        "Truth bookkeeping only — never a trigger cut");
    pm->SetDefaultParameter("eventbuilder:findable:eta_abs", m_acc_eta_abs,
        "in-acceptance expected-track count (weights[11]): maximum MC |eta|. "
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
        "weights[23] (empirical S) and weights[24] (band Fano) are stored in every mode");
    pm->SetDefaultParameter("eventbuilder:topology:mu_min_hits", m_mu_min_hits,
        "minimum frame hits in a cell (or band) for its own mu estimate to be used; "
        "below this the ladder falls back to the band (then global) estimate");
    pm->SetDefaultParameter("eventbuilder:topology:mu_loo", m_mu_loo,
        "1 = leave-one-out scoring null (stored S/p_tail evaluated against mu computed "
        "without the tested window's own hits, removing self-masking circularity); "
        "0 = legacy frame-wide null. Trigger thresholds are frame-wide either way");
    // Tracklet-pointing windows (weights[10], score-only). Sized to absorb
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
    // Trigger-tag bitmask thresholds (weights[16], inclusion tags, score-only).
    pm->SetDefaultParameter("eventbuilder:tags:zone_nhits", m_tag_zone_nhits,
        "tag bits 0/2/4: min cluster nhits for an ECal zone region tag");
    pm->SetDefaultParameter("eventbuilder:tags:b0_min_hits", m_tag_b0_min_hits,
        "tag bit 6: min in-window B0 hits for the far-forward tag");
    pm->SetDefaultParameter("eventbuilder:tags:zdc_min_hits", m_tag_zdc_min_hits,
        "tag bit 7: min summed in-window ZDC-ECal cluster hits");
    // STAGE 1 trigger TRK/CAL thresholds (weights[12]/[14], score-only —
    // see the member declarations above for why these exist unapplied).
    pm->SetDefaultParameter("eventbuilder:trigger:min_tracklets", m_trig_min_tracklets,
        "TRK term of offline TIME && (TRK || CAL): min n_tracklets (weights[10]) "
        "to count as tracker content. Not applied to phys_flag/emission.");
    pm->SetDefaultParameter("eventbuilder:trigger:min_cal_energy", m_trig_min_cal_energy,
        "CAL term of offline TIME && (TRK || CAL): min E_calo (GeV, weights[9]) "
        "to count as calo content. 0 = off (UNMEASURED PLACEHOLDER). Not applied "
        "to phys_flag/emission.");
    pm->SetDefaultParameter("eventbuilder:trigger:cal_nsigma", m_cal_nsigma,
        "N-sigma half-width for the per-cluster gate feeding e_calo (weights[9]): "
        "gate = N x cluster.getTimeError(). Mirrors CalTimeCoincidence::nsigma.");
    pm->SetDefaultParameter("eventbuilder:trigger:cal_late_nsigma_asym", m_cal_late_nsigma_asym,
        "extra LATE acceptance for the e_calo per-cluster gate, in units of "
        "that cluster's own sigma, added to the + side only. Mirrors "
        "CalTimeCoincidence::late_nsigma_window_asym. UNMEASURED PLACEHOLDER "
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
        "matches weights[13]/[15] without affecting output.");
    // A coincidence needs at least 2 hits: fewer cannot form a cluster, and
    // t0 over 0 hits is undefined. Treat a threshold of 0 or 1 as "no
    // topology cut" and clamp it to 2, the minimum pair, so behavior and
    // the stored value (EventBuilderFrameInfo weights[5]) stay consistent.
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
    // Cylindrical position and barrel/endcap split. Used only by the
    // tracklet-pointing score at emit time (weights[10]), never by the
    // trigger. Barrel layers: TOFBarrel(0), MPGDBarrel(2),
    // OuterMPGDBarrel(3). Endcap disks: TOFEndcap(1), Backward(4)/Forward(5)
    // MPGDEndcap.
    float r_cyl, z_pos, phi_raw;
    bool is_barrel;
  };

  // Incremental per-cell windowed-count tracker: one instance per (grid,
  // cell); see tof1/tof2 in Process(). It tracks the best hit count
  // achievable in any `window`-wide span, among the hits currently active
  // in this cell (within the outer scan's [lo, hi) range), and updates that
  // count on Insert()/EvictFront() as hits enter or leave. Reading `best` is
  // O(1). A caller's combined max across all cells (for example, cmax_tof
  // in Process()) then costs O(ncells) per query. This class is generic in
  // the hit pool it tracks (today, the TOF-only pool); a future second
  // detector group with its own window width can reuse it unchanged.
  //
  // Correctness argument:
  //  - Insert(idx) treats idx as a candidate window END, looking backward:
  //    it advances back_lo (the earliest active hit still within `window`
  //    of idx) forward as needed. The window ending at idx then has count
  //    active.size() - back_lo. Recording the max of this over every
  //    insertion covers the same candidate windows as a forward scan (try
  //    every start, grow forward), just walked in the opposite direction,
  //    so both give the same true max.
  //  - EvictFront() removes the earliest active hit in this cell, mirroring
  //    the outer scan's lo advancing. If that hit was not part of the
  //    window achieving `best` (best_pos > 0), `best` is unaffected and its
  //    index only shifts down by one. If it was part of that window
  //    (best_pos == 0), `best` may no longer be achievable, so Recompute()
  //    reruns a two-pointer scan, now bounded by this cell's own active
  //    size instead of the whole window's hit count.
  //
  //  Amortized cost over one frame: O(n) for Insert, plus the sum of every
  //  Recompute()'s cost. When background spreads roughly evenly across
  //  cells, that sum stays small, since each cell's active size is a
  //  fraction of n. If nearly all hits land in one cell, Recompute() can
  //  fire often and this degrades toward O(n) per query; it is a real
  //  improvement in the typical case, not a guarantee independent of how
  //  background is distributed.
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
          continue; // unusable resolution (unset or zero); would blow up the 1/sigma^2 weight
        if (min_edep() > 0.0f && hit.getEdep() < min_edep())
          continue; // soft background hit; excluded from event finding only
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
    // dt and dt_tof are not fixed constants. Each is the widest pairwise
    // time gap that two hits of this frame's own measured resolution could
    // genuinely differ by, at nsigma_window confidence. For a pair with
    // resolutions sigma_a and sigma_b, the matched width is
    // sqrt(sigma_a^2 + sigma_b^2), which is largest when both hits sit at
    // the pool's own worst resolution. So max(sigma) times sqrt(2) is the
    // exact worst-case pairwise width, not an approximation. max_sigma_all
    // is dominated by MPGD (about 10 ns) whenever any MPGD hit is present,
    // which background alone guarantees in practice. max_sigma_tof is the
    // TOF-only subset's own resolution (about 30 ps). Both are measured
    // directly from this frame's hits, so this adds no parameter beyond
    // nsigma_window itself.
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

    // -- MC collision times (for the flag count and mc_t) -----------------------
    // One unified loop over every real-collision final-state particle.
    // Every physics class sits in its own 1000-wide band at
    // 10000 + class_index*1000, so a stable primary of class k is
    // 10000 + k*1000 + 1 and the top class ends at 35999 -- inside the
    // 16-bit ceiling DDG4 imposes on generatorStatus (see status_bands.py).
    // Native status 1 (a tagged primary from a legacy single-class run) is
    // still accepted and carries stream 0. `stream` is status/1000 of the
    // particle's own status: 0 = native primary, >= 10 = a physics class
    // (class_index = stream - 10). It is only ever used as a per-class
    // grouping key here, so it works unchanged for the legacy 10000-wide
    // layout too. Machine background (status < 10000, bands 2000-6999)
    // never enters this list. Every stable primary of one
    // collision carries that collision's frame time, so this dedupes by
    // (stream, time) within the frame-derived coincidence window to get one
    // entry per injected collision; two different classes' collisions
    // landing close in time must not merge into one entry.
    std::vector<double> collision_times;
    std::vector<int> collision_nexp; // per-collision in-acceptance stable charged
                                     // (the n_expected count, weights[11])
    std::vector<int> collision_stream;
    if (m_mc_in() != nullptr) {
      const double dedup = double(dt);
      for (const auto& mcp : *m_mc_in()) {
        const int status = mcp.getGeneratorStatus();
        const bool native_primary = (status == 1);
        const bool physics_stable = (status >= 10000 && status % 1000 == 1);
        if (!native_primary && !physics_stable)
          continue; // not a real-collision final-state particle
        const int stream = native_primary ? 0 : status / 1000;
        const double t = mcp.getTime();
        int ci = -1;
        for (size_t j = 0; j < collision_times.size(); ++j)
          if (collision_stream[j] == stream && std::abs(collision_times[j] - t) < dedup) {
            ci = int(j);
            break;
          }
        if (ci < 0) {
          collision_times.push_back(t);
          collision_nexp.push_back(0);
          collision_stream.push_back(stream);
          ci = int(collision_times.size()) - 1;
        }
        if (std::abs(mcp.getCharge()) > 0.01) {
          const auto p    = mcp.getMomentum();
          const double pt = std::hypot(double(p.x), double(p.y));
          if (pt > m_acc_pt_min &&
              std::abs(std::asinh(double(p.z) / pt)) < m_acc_eta_abs)
            ++collision_nexp[ci];
        }
      }
    }

    // Calorimeter clusters of this frame, sorted by time. Used only at
    // emit time, for the windowed energy sum (weights[24]) and the
    // trigger-tag bitmask (weights[16]). `coll` indexes
    // m_calo_collection_names (1 = barrel ECal, 2 = backward ECal, 3 =
    // forward ECal, 7 = ZDC ECal: the zones the bitmask reads). nhits is
    // the cluster hit count, the regional hit multiplicity the bitmask
    // thresholds on.
    struct CaloCl {
      double t, e;
      float sigma; // clu.getTimeError(): the calo-specific per-cluster
                   // resolution that e_calo's gate uses below, not dt
                   // (which is TOF/MPGD hit resolution; see max_sigma_all).
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
    std::sort(calo_te.begin(), calo_te.end()); // kept sorted for the debug_calo dump only
    if (debug_calo()) {
      // Opt-in diagnostic (eventbuilder:trigger:debug_calo=1): dumps the
      // frame's time-aligned clusters, so cluster-time problems (for
      // example, a cluster whose time is unset and lands at -r/c after
      // alignment) are visible without persisting Timeslice-level
      // collections.
      jout << "EB debug_calo frame=" << event_number << " nclusters=" << calo_te.size() << jendl;
      for (const auto& c : calo_te)
        if (c.e > 0.1)
          jout << "  coll=" << c.coll << " t=" << c.t << " E=" << c.e
               << " nhits=" << c.nhits << jendl;
    }

    // B0 hit times (sorted). Tracker input index 9. Not part of event
    // finding; used only by the bitmask's far-forward tag (bit 6).
    std::vector<double> b0_times;
    if (m_trk_in().size() > 9 && m_trk_in().at(9) != nullptr)
      for (const auto& hit : *m_trk_in().at(9))
        b0_times.push_back(double(hit.getTime()));
    std::sort(b0_times.begin(), b0_times.end());

    // -- 3. O(N) sliding-window coincidence + (theta,phi) topology trigger -----
    // Per-frame background density estimate, used by the adaptive threshold
    // and by the stored per-candidate significance in both modes. Signal
    // contaminates mu by only a few percent (one cluster among roughly
    // dT/dt window positions), so this is ignored.
    double mu = 0.0, mu_tof = 0.0;
    double dT = double(dt);
    if (hits.size() >= 2) {
      dT     = std::max(double(hits.back().time - hits.front().time), double(dt));
      mu     = (double(hits.size()) * double(dt) / dT) / double(ncells);
      mu_tof = (double(n_tof) * double(dt_tof) / dT) / double(ncells);
    }

    // Inhomogeneous Poisson background model. The grid keeps equal solid
    // angle, but equal solid angle does not mean equal expected background:
    // measured background is peaked toward the beam axis, varying about
    // 19x across theta rows. So the Poisson mean is allowed to vary with
    // detector position. Each cell of each staggered grid carries its own
    // local null,
    //   H_{0,i}: X_i ~ Poisson(mu_i),   mu_i = N_i * dt / dT,
    // where N_i is this frame's hit count in cell i, normalized to the
    // coincidence window, and, in adaptive mode, its own threshold k_i.
    //
    // These frame-wide mu_i values drive the trigger thresholds only, where
    // the tested window's self-contamination is small (about 2%) and does
    // not affect the firing decision. The stored significance of each
    // candidate is evaluated against a separate, leave-one-out estimate
    // instead: the same ladder recomputed with the candidate's own window
    // excluded, so no candidate is ever scored against a background
    // estimate it contributed to. A cross-frame running background, using
    // TimesliceBuffer_service, is the natural upgrade if per-frame
    // statistics ever become limiting; see the README Outlook.
    //
    // Statistics ladder (topology:mu_min_hits):
    //   N_i           >= mu_min_hits  -> per-cell mu_i
    //   N_band(theta) >= mu_min_hits  -> cos-theta-band mu/nphi
    //   otherwise                     -> global mu
    // At a 1x1 grid, every path reduces to the global mu.
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

    // Sideband/dispersion statistics (eventbuilder:statistics), computed
    // here so the adaptive k-thresholds just below can reuse the same
    // per-band dispersion that weights[23]/[24] are stored from later,
    // instead of computing it twice. This is the windowed-count
    // distribution of the frame itself: disjoint dt-wide windows tiling
    // [t_first, t_last], counted per grid-1 cell, histogrammed per
    // cos-theta band, including empty windows (they carry most of the
    // probability mass when mu is much less than 1). Correlated hit
    // families (conversions, showers) land in the sideband tail exactly as
    // they land in candidate windows, so this prices them correctly where a
    // plain Poisson tail cannot. A candidate's own hits occupy only a
    // handful of the nphi x nwin windows of its band, so no leave-one-out
    // correction is needed at this sample size. This is always computed,
    // not only in adaptive mode, since weights[23]/[24] are stored every
    // frame regardless.
    //
    // dt can be 0 when the frame has fewer than 2 fast hits (max_sigma_all
    // stays at its default of 0; see step 2 above), which would make dT/dt
    // a 0/0 division. This is guarded explicitly, rather than relying on
    // std::max(1, ...) to mask a NaN cast to int, which is undefined
    // behavior.
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

    // Per-cell applied thresholds k_i: in adaptive mode, each cell's own
    // P(X_i >= k_i) < p_fake, floored at topology:threshold; in fixed mode,
    // every cell uses the configured floor. k1/k2 (combined-group; these are
    // the ones that actually gate firing, in the scan below) use a
    // negative-binomial tail: mean is the cell/band/global mu ladder
    // (mu1[i]/mu2[i], unchanged), and variance is band_fano[band] * mean.
    // The mean keeps the ladder's resolution, but the width is corrected by
    // the separately measured local overdispersion (Fano; see weights[24]).
    // A plain Poisson k understates the true accidental rate wherever
    // background is bursty (Fano > 1, from conversion/shower hit families).
    // nbinomMinK reduces exactly to poissonMinK wherever the local Fano is
    // 1 or less, so this only ever widens the threshold, never narrows it.
    //
    // k_tof (the global scalar) stays on plain Poisson: its pool is
    // comparatively sparse, and its own dispersion has not been separately
    // measured. The applied per-cell k1[i]/k2[i] gate the combined-group
    // firing; the winning cell's k_i is stored per candidate (weights[19]).
    const bool adaptive = adaptiveMode();
    int k_tof = std::numeric_limits<int>::max(); // TOF group fires only in adaptive mode
    std::vector<int> k1(ncells, m_threshold), k2(ncells, m_threshold);
    if (adaptive) {
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

    // -- Frame-constant record (EventBuilderFrameInfo) --------------------------
    // One entry per frame, emitted once dt/dt_tof, mu_tof/k_tof and the
    // injected-collision list are all known. Layout documented above
    // emitFrameInfo() below.
    emitFrameInfo(run_number, int(collision_times.size()), mu_tof, k_tof, dt, dt_tof);

    // Candidates found by the scan, buffered so MC matching can be done
    // globally (nearest candidate per collision) once every candidate of the
    // frame is known — see the matching pass after the scan loop.
    struct Pending {
      double t0, t0sigma, S, S_tof;
      int cmax, cmax_tof, win_theta, win_phi, win_grid;
      double e_calo;
      int n_tracklets, tagbits;
      double mu_win;  // winning cell's LEAVE-ONE-OUT scoring null (weights[18]) —
                      // the mu that S/p_tail were computed against
      int k_win;      // winning cell's applied trigger threshold (weights[19])
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

      // Also track the winning cell (grid, theta bin, phi bin), stored per
      // candidate so offline diagnostics can build the empirical occupancy
      // map of where triggers actually fire (see weights[20..22]). The
      // winning cell is the most significant one under its own local null
      // H_{0,i} (the smallest P(X >= occ_i | mu_i)), not the cell with the
      // highest raw count: with non-uniform mu_i, a hot beam-axis cell
      // needs more hits than a quiet one to be equally surprising. The
      // trigger fires when any cell exceeds its own k_i. With mu_mode=0 (or
      // a 1x1 grid), mu_i is the same everywhere and this reduces to the
      // simple max-count behavior.
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
      // cmax_tof: the max over every cell's already-maintained `best`
      // (O(1) each), not a fresh rescan; see CellTracker.
      int cmax_tof = 0;
      for (int i = 0; i < ncells; ++i) {
        cmax_tof = std::max(cmax_tof, tof1[i].best);
        cmax_tof = std::max(cmax_tof, tof2[i].best);
      }
      const bool trigger = fired || (cmax_tof >= k_tof);

      if (trigger) {
        double sumw = 0.0, sumwt = 0.0;
        if (m_t0_estimator == "leading") {
          // Leading edge: sigma_min identifies the most precise detector
          // class in the window. The anchor is the earliest hit of that
          // class, since the fast population arrives first and the late
          // tail is one-sided. Only hits whose time is compatible with the
          // anchor, within their own resolution, are averaged; stragglers
          // (slow hadrons, curlers) are excluded instead of dragging t0
          // late at full weight.
          double sigma_min = std::numeric_limits<double>::infinity();
          for (size_t j = lo; j < hi; ++j)
            sigma_min = std::min(sigma_min, double(hits[j].sigma));
          // Coincidence-confirmed anchor: the earliest precise hit that has
          // at least one other precise hit within 3*(sigma_i+sigma_j). This
          // stops a lone early accidental hit from hijacking the anchor.
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
          }  // no confirmed anchor: sums stay 0, so the mean fallback below runs
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

        // S/S_tof: -10*log10(p_tail), a dB-scale significance built
        // directly from the exact one-sided Poisson tail, valid at any mu.
        // Deliberately not a z-score/sigma count: it is a monotone
        // rarity score on a log scale, well-defined even at the routinely
        // sub-1 mu this pipeline operates at, where a Gaussian z-score does
        // not apply. Computed with poissonLogSurvival10 (log-space upper
        // tail), not -log10(poissonSurvival(...)): the 1-cdf form saturates
        // at double precision (p floors around 2e-16, so S caps at about
        // 157, and exact zeros clamp to a 1e-300 floor, piling S up at
        // 3000), which would destroy the ranking of every strong candidate.
        //
        // Local null, leave-one-out: the stored significance is evaluated
        // against the winning cell's mu_i, recomputed without this window's
        // own hits (counts minus the window's contribution, off-window time
        // span dT - dt), so the tested window never contributes to the null
        // it is tested against. The frame-wide mu_win/k_win (trigger side)
        // remain whatever fired; only the scoring null is decontaminated.
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
        // Calorimeter-coincidence energy (weights[9]) and the per-zone
        // cluster tags feeding the bitmask, over the time-sorted cluster
        // list. A zone tag is set when any in-window cluster of that zone's
        // ECal has nhits >= m_tag_zone_nhits: a cluster with that many hits
        // is a shower. The ZDC tag instead thresholds the summed in-window
        // ZDC-ECal cluster hit count, since ZDC showers fragment across
        // clusters.
        double e_calo = 0.0;
        bool zone_bwd = false, zone_bar = false, zone_fwd = false;
        int zdc_nhits = 0;
        // Per-cluster gate: N sigma times this cluster's own
        // getTimeError(), plus a late-side asymmetric allowance, the same
        // formula as CalTimeCoincidence (see m_cal_nsigma's declaration
        // comment for why this is applied inline instead). The energy
        // floor is split barrel (coll==1 EcalBarrel, coll==4 HcalBarrel)
        // versus forward/backward, the same indices the zone tags below use.
        for (const auto& c : calo_te) {
          if (!(c.sigma > 0.0f))
            continue; // unusable resolution (unset or zero)
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
        // Vertex-pointing tracklets (weights[10]) from the gated fast hits
        // of this window [lo,hi). Score only; never part of the trigger
        // decision. Barrel pairs use straight-line z0 extrapolation to the
        // beamline. Endcap pairs (same side) use projective consistency,
        // r2 ~ r1*z2/z1, since the short two-plane lever arm makes z0
        // itself ill-conditioned there. Window sizes absorb solenoid
        // curvature (about 0.05-0.1 rad at 0.5 GeV), not detector
        // resolution.
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
        // Trigger bitmask (`trigger`, weights[16]): inclusion tags, each an
        // independent positive excess. This stored field is never read by
        // the trigger decision that produced the candidate. Bits 0/1 are
        // the trigger-origin bits (which hit group fired); bits 2-9 are
        // the detector-zone tags (formerly bits 0-7).
        //   bit 0: TOF-only group fired: cmax_tof >= k_tof
        //   bit 1: combined (TOF+MPGD) group fired: cmax >= k_win
        //   bit 2: backward-ECal region tag: in-window backward-ECal cluster, nhits >= 10
        //   bit 3: bit 2, plus >= 1 backward endcap tracklet
        //   bit 4: barrel-ECal region tag
        //   bit 5: bit 4, plus a barrel tracklet
        //   bit 6: forward-ECal region tag
        //   bit 7: bit 6, plus a forward endcap tracklet
        //   bit 8: far-forward B0 tag: >= 4 in-window B0 hits
        //   bit 9: ZDC tag: in-window ZDC-ECal cluster hits sum >= 50
        // Combined tags are computed offline from these bits; none are stored:
        //   cTag1 = (#set{2,4,6} >= 2) && bit8      cTag2 = (#set{3,5,7} >= 2) && bit8
        //   cTag3 = (#set{2,4,6} >= 2) && bit9      cTag4 = (#set{3,5,7} >= 2) && bit9
        //   cTag5 = (#set{2,4,6} == 3)              cTag6 = (#set{3,5,7} >= 2)
        int tagbits = 0;
        if (cmax_tof >= k_tof)         tagbits |= 1 << 0;
        if (cmax >= k_win)             tagbits |= 1 << 1;
        if (zone_bwd)                  tagbits |= 1 << 2;
        if (zone_bwd && trk_bwd > 0)   tagbits |= 1 << 3;
        if (zone_bar)                  tagbits |= 1 << 4;
        if (zone_bar && trk_barrel > 0) tagbits |= 1 << 5;
        if (zone_fwd)                  tagbits |= 1 << 6;
        if (zone_fwd && trk_fwd > 0)   tagbits |= 1 << 7;
        if (n_b0 >= m_tag_b0_min_hits)        tagbits |= 1 << 8;
        if (zdc_nhits >= m_tag_zdc_min_hits)  tagbits |= 1 << 9;

        // Defer emission: MC matching happens globally after the scan
        // (each collision goes to its nearest candidate), so no flag is
        // decided here. Matching greedily during the scan would mislabel
        // real events whenever an earlier background window also sat
        // within dt of the collision.
        pending.push_back({t0, t0sigma, S, S_tof,
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
    // band_fano, empiricalTail) were computed earlier, in step 3's "Per-cell
    // applied thresholds", so the adaptive-mode k-thresholds there could
    // reuse the same per-band dispersion that weights[23]/[24] are stored
    // from, instead of computing it twice.

    // -- 4. Global MC matching (inclusive per-candidate count), then emission --
    // flag (weights[2]) = the number of real collisions within dt of this
    // candidate's own t0: 0 = fake/accidental, 1 = normal, 2+ = true
    // pile-up. mc_t (weights[3]) = the EARLIEST matched collision's time
    // (deterministic when several overlap); n_expected (weights[11]) is
    // that earliest collision's in-acceptance track count. This replaces
    // the old exclusive nearest-match 1/2/3 scheme, which could only name
    // one collision per candidate and whose "tagged primary" branch never
    // fired under unified mixing.
    std::vector<int> flags(pending.size(), 0);
    std::vector<double> mcts(pending.size(), kNoMCTime);
    std::vector<int> nexps(pending.size(), 0);
    // Which class(es) coincided (trigger_classes_mask, weights[17]): bit N
    // = class_index N (0-25), where class_index = stream - 10 for
    // physics streams. A native primary (stream 0) sets no bit. Decode a
    // bit to a class name with the class table in AGENTS.md.
    std::vector<uint32_t> classes_mask(pending.size(), 0);
    // Full (time, stream) pair list (trigger_classes, weights[25+]); see
    // store_coincident_list. Read through the ParameterRef accessor, not
    // the bound m_ member directly (matching min_edep()/debug_calo()
    // elsewhere in this file), since JANA's parameter-strictness check
    // needs the accessor call to register this parameter as read.
    const bool store_list = store_coincident_list();
    std::vector<std::vector<std::pair<double, int>>> coincident_list(pending.size());
    for (size_t i = 0; i < pending.size(); ++i) {
      const double t0i = pending[i].t0;
      int cnt = 0;
      uint32_t smask = 0;
      double earliest = kNoMCTime;
      int earliest_nexp = 0;
      for (size_t m = 0; m < collision_times.size(); ++m) {
        if (std::abs(t0i - collision_times[m]) > double(dt))
          continue;
        ++cnt;
        const int stream = collision_stream[m];
        if (stream >= 10) {
          // stream = status / 1000 = 10 + class_index (base 10000, 1000-wide
          // bands, see status_bands.py's stream_to_class_index) -- NOT the
          // abandoned 10000-wide layout, so no extra /10 here. That stray
          // /10 collapsed classes 0-9/10-19/20-25 into three buckets
          // (0/1/2), the same class of bug already fixed on the Python side
          // this session (data_stats.py's cls_inj and cls_hits //10000
          // bugs) -- this was the matching bug on the C++ emit side.
          const int cls = stream - 10;
          if (cls >= 0 && cls < 26)
            smask |= (1u << cls);
        }
        if (earliest == kNoMCTime || collision_times[m] < earliest) {
          earliest      = collision_times[m];
          earliest_nexp = collision_nexp[m];
        }
        if (store_list)
          coincident_list[i].emplace_back(collision_times[m], stream);
      }
      flags[i]        = cnt;
      mcts[i]         = earliest;
      nexps[i]        = earliest_nexp;
      classes_mask[i] = smask;
    }
    for (size_t i = 0; i < pending.size(); ++i) {
      const Pending& c = pending[i];
      const int band     = std::min(std::max(c.win_theta, 0), ntheta - 1);
      const double p_emp = empiricalTail(c.cmax, band);
      const double s_emp = -10.0 * std::log10(p_emp);
      const double fano  = band_fano[band];
      // Mode-selected stored significance/tail. The trigger decision is
      // already made by this point; this only picks what gets stored.
      // Bare names alias to their canonical granularity: poisson->cell,
      // nbinom->band, empirical->band.
      double s_stored = c.S;
      // c.S was computed with poissonLogSurvival10 (log-space,
      // non-saturating) against the same (cmax, mu_win) pair used here.
      // Deriving p_stored from it (p = 10^(-S/10)) keeps the two fields
      // consistent; calling poissonSurvival(c.cmax, c.mu_win) directly
      // would saturate around 2e-16 and disagree with c.S for any candidate
      // significant beyond that floor.
      double p_stored = std::pow(10.0, -s_stored / 10.0);
      if (m_stat_null == "empirical" || m_stat_null == "empirical_band") {
        s_stored = s_emp;
        p_stored = p_emp;
      } else if (m_stat_null == "nbinom" || m_stat_null == "nbinom_band") {
        p_stored = nbinomSurvival(c.cmax, band_mean[band], band_var[band]);
        s_stored = -10.0 * std::log10(std::max(p_stored, 1e-300));
      } else if (m_stat_null == "nbinom_cell") {
        // Hybrid: per-cell leave-one-out mean, band-measured dispersion.
        // Falls back to the Poisson tail automatically wherever
        // Fano_band <= 1.
        p_stored = nbinomSurvival(c.cmax, c.mu_win, fano * c.mu_win);
        s_stored = -10.0 * std::log10(std::max(p_stored, 1e-300));
      }
      // Early veto (eventbuilder:trigger:veto, off by default; see the
      // member declaration above). TIME already kept this window (it is in
      // `pending`), but if neither TRK nor CAL clears its threshold, this
      // drops it before emitCandidate() creates the EventCandidates row.
      // flags[i]/mcts[i]/nexps[i], already assigned by the MC-matching pass
      // above, are simply discarded for this i. A flags[i] > 0 candidate
      // vetoed here is a real collision that is genuinely lost, not
      // recoverable offline.
      if (m_trig_veto && c.n_tracklets < m_trig_min_tracklets &&
          c.e_calo < m_trig_min_cal_energy) {
        continue;
      }
      emitCandidate(run_number, c.t0, c.t0sigma, flags[i], mcts[i], s_stored,
                    c.S_tof, p_stored, c.cmax, c.cmax_tof, c.e_calo,
                    c.n_tracklets, nexps[i], c.tagbits, classes_mask[i],
                    c.mu_win, c.k_win, c.win_theta, c.win_phi, c.win_grid,
                    s_emp, fano, coincident_list[i]);
    }
  }

  // ---------------------------------------------------------------------------
  // Sentinel for "no matched MC collision" in the candidate encoding (ns).
  static constexpr double kNoMCTime = -1.0e9;

  // Tracklet-pointing windows (weights[10]) and trigger-tag thresholds
  // (weights[16]) are runtime parameters; see Configure():
  //   eventbuilder:tracklet:{dphi_barrel,dphi_endcap,dr_min,dz_min,z0_max,r_tol}
  //   eventbuilder:tags:{zone_nhits,b0_min_hits,zdc_min_hits}
  // Defaults were tuned offline on simulated quantities; retuning on
  // reconstructed quantities is pending. All score-side; never the trigger.

  // MC matching happens globally in Process() step 4: nearest candidate per
  // collision.

  // Smallest k with P(X >= k | Poisson(mu)) < p: the adaptive-mode count threshold.
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

  // Negative-binomial upper tail P(X >= c), moment-fitted to a mean m and
  // variance v: the count law of Poisson-arriving hit families
  // (overdispersed, Fano = v/m > 1). Reduces to the Poisson tail when the
  // measured dispersion is consistent with independence (v <= m).
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

  // Smallest k with P(X >= k | NegBinom(mean, var)) < p: the honest analog
  // of poissonMinK once background is measurably overdispersed
  // (Fano = var/mean > 1). Reduces to poissonMinK exactly wherever the
  // local dispersion is consistent with independence, so this only ever
  // widens the threshold, never gives a lower k than plain Poisson would at
  // the same (mean, p).
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

  // Exact upper-tail probability P(X >= c | Poisson(mu)): the honest
  // rarity of a c-hit clump under background, valid at any mu. Stored raw
  // as weights[23], so offline cuts have the untransformed tail
  // probability. The trigger's own firing decision already uses this tail,
  // through poissonMinK. Returns 1.0 for c <= 0, and clamps to a tiny floor
  // so log10(p) downstream never hits negative infinity.
  //
  // As a probability, this saturates at double precision: 1-cdf
  // cancellation floors p around 2e-16. The stored significance S
  // therefore uses poissonLogSurvival10 below instead, which stays exact
  // far past that floor.
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

  // log10 P(X >= c | Poisson(mu)), computed entirely in log space: the
  // basis of the stored significance S = -10*log10(p_tail). Unlike
  // poissonSurvival's 1-cdf, which cancels once the tail drops below about
  // 2e-16 (capping S around 157), this sums the upper-tail series directly
  // from k=c:
  //   P(X>=c) = pmf(c) * [1 + mu/(c+1) + mu^2/((c+1)(c+2)) + ...]
  // with log pmf(c) = -mu + c*ln(mu) - lgamma(c+1). The bracket converges
  // geometrically, and the result stays exact to double precision at any
  // c. c <= 0 returns 0.0 (p=1). mu <= 0 with c > 0 returns -300, the same
  // floor poissonSurvival uses for an empty background.
  static double poissonLogSurvival10(int c, double mu) {
    if (c <= 0)
      return 0.0;
    if (!(mu > 0.0))
      return -300.0;
    // For mu >= c the tail is O(1): the direct 1-cdf has no cancellation
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

  // Candidate encoding, in EventCandidates. EventUnfolder copies every entry
  // verbatim into each child's EventBuilderInfo. Keep external readers in
  // sync with this list. Frame-constant configuration lives in
  // EventBuilderFrameInfo (see emitFrameInfo() below), not here.
  //
  //   timeStamp   = round(t0_ns * 1000), t0 in ps, as an integer.
  //   weights[0]  = t0, in ns.
  //   weights[1]  = t0sigma, in ns.
  //   weights[2]  = flag (MC only): the number of real collisions within dt
  //                 of this candidate's own t0. 0 = fake/accidental,
  //                 1 = one real collision (normal), 2+ = true pile-up.
  //                 Real means flag > 0.
  //   weights[3]  = mc_t: the EARLIEST matched collision's time, in ns
  //                 (MC only; kNoMCTime if none). When flag >= 2 the other
  //                 collisions' times are in trigger_classes (weights[25+]).
  //   weights[4]  = S: combined-group (dt) significance, S = -10*log10(p_tail).
  //                 A dB-scale rarity score, not a Gaussian sigma count.
  //   weights[5]  = S_tof: TOF-group (dt_tof) significance.
  //   weights[6]  = p_tail = P(X >= cmax | Poisson(mu)): the exact upper-tail
  //                 background probability behind weights[4]'s S.
  //   weights[7]  = cmax: the combined-group coincidence count that fired.
  //   weights[8]  = cmax_tof: the TOF-group coincidence count.
  //   weights[9]  = E_calo, in GeV: summed time-aligned cluster energy with
  //                 |t_cluster - t0| <= dt. Score only.
  //   weights[10] = n_tracklets: vertex-pointing fast-hit pairs in the fired
  //                 window. Score only.
  //   weights[11] = n_expected: in-acceptance stable charged MC particles of
  //                 the earliest matched collision. The tracking-efficiency
  //                 denominator. 0 for fakes and MC-less runs.
  //   weights[12] = trk_threshold_cfg (eventbuilder:trigger:min_tracklets).
  //   weights[13] = trk_pass (1.0/0.0): n_tracklets >= weights[12].
  //   weights[14] = cal_threshold_cfg (eventbuilder:trigger:min_cal_energy, GeV).
  //   weights[15] = cal_pass (1.0/0.0): E_calo >= weights[14].
  //   weights[16] = trigger: bitmask, bits 0-9 (see Process()): bit 0 = TOF
  //                 group fired, bit 1 = combined group fired, bits 2-9 =
  //                 detector-zone tags. Inclusion tags only.
  //   weights[17] = trigger_classes_mask: bitmask, bit N = class_index N
  //                 (0-25) had a collision within dt of this candidate.
  //                 0 = none (fake, or a native-primary-only match).
  //   weights[18] = mu: winning cell's leave-one-out background
  //                 (hits per cell per window). Per candidate.
  //   weights[19] = k: winning cell's applied count threshold. Per candidate.
  //   weights[20] = winning cell's theta bin.
  //   weights[21] = winning cell's phi bin.
  //   weights[22] = winning grid: 1 = plain, 2 = half-cell staggered.
  //   weights[23] = S_empirical: significance of cmax against this frame's
  //                 own sideband windowed-count distribution. Always stored.
  //   weights[24] = Fano factor (variance/mean of sideband windowed counts)
  //                 of the winning cell's cos-theta band. 1 = Poisson-like.
  //   weights[25+] = trigger_classes (on by default through
  //                 store_coincident_list; the vector ends at [24] when
  //                 off). weights[25] = N, the number of coincident
  //                 collisions (same value as flag). N (time, stream) pairs
  //                 follow, at weights[26+2k] and weights[27+2k] for k in
  //                 [0, N). stream = status/1000 of the collision's base
  //                 status: 0 = native primary, >= 10 = a physics class
  //                 (class_index = stream-10). This lets a consumer
  //                 attribute individual hits or MCParticles to the correct
  //                 collision when two or more genuinely overlap one
  //                 candidate; see eventbuilder.cc's MCParticle/hit gating.
  //
  // weights[4] and weights[6] always carry the significance and tail
  // probability of the active eventbuilder:statistics:null mode. The
  // default mode is "poisson".
  void emitCandidate(int64_t run_number, double t0, double t0sigma, int flag,
                     double mc_t, double signif, double signif_tof, double p_stored,
                     int cmax, int cmax_tof, double e_calo, int n_tracklets,
                     int n_expected, int tagbits, uint32_t classes_mask,
                     double mu, int k, int win_theta, int win_phi, int win_grid,
                     double s_emp, double fano,
                     const std::vector<std::pair<double, int>>& coincident_list) {
    auto cand = m_candidates_out()->create();
    cand.setRunNumber(static_cast<uint32_t>(run_number));
    cand.setEventNumber(m_candidates_out()->size() - 1);
    cand.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    cand.setWeight(t0sigma);                                             // ns (legacy)
    cand.addToWeights(t0);                                               // [0]
    cand.addToWeights(t0sigma);                                          // [1]
    cand.addToWeights(static_cast<double>(flag));                        // [2]
    cand.addToWeights(mc_t);                                             // [3]
    cand.addToWeights(signif);                                           // [4]
    cand.addToWeights(signif_tof);                                       // [5]
    cand.addToWeights(p_stored);                                         // [6]
    cand.addToWeights(double(cmax));                                     // [7]
    cand.addToWeights(double(cmax_tof));                                 // [8]
    cand.addToWeights(e_calo);                                           // [9]
    cand.addToWeights(double(n_tracklets));                              // [10]
    cand.addToWeights(double(n_expected));                               // [11]
    cand.addToWeights(double(m_trig_min_tracklets));                     // [12]
    // Precomputed pass/fail, not just the raw threshold above: a future
    // online early-veto can read this boolean directly, instead of
    // re-deriving the comparison from weights[9]/[10] vs [12]/[14].
    cand.addToWeights(n_tracklets >= m_trig_min_tracklets ? 1.0 : 0.0);  // [13]
    cand.addToWeights(m_trig_min_cal_energy);                            // [14]
    cand.addToWeights(e_calo >= m_trig_min_cal_energy ? 1.0 : 0.0);      // [15]
    cand.addToWeights(double(tagbits));                                  // [16]
    cand.addToWeights(double(classes_mask));                             // [17]
    cand.addToWeights(mu);                                               // [18]
    cand.addToWeights(double(k));                                        // [19]
    cand.addToWeights(double(win_theta));                                // [20]
    cand.addToWeights(double(win_phi));                                  // [21]
    cand.addToWeights(double(win_grid));                                 // [22]
    cand.addToWeights(s_emp);                                            // [23]
    cand.addToWeights(fano);                                             // [24]
    // [25+] trigger_classes, on by default through store_coincident_list;
    // see that parameter's own comment and the layout block above.
    if (store_coincident_list()) {
      cand.addToWeights(double(coincident_list.size()));                 // [25]
      for (const auto& [ctime, cstream] : coincident_list) {
        cand.addToWeights(ctime);                                        // [26+2k]
        cand.addToWeights(double(cstream));                              // [27+2k]
      }
    }
  }

  // Frame-constant configuration, in EventBuilderFrameInfo: one entry per
  // frame, holding the values that are identical for every candidate of the
  // frame (formerly duplicated into each candidate's weights). EventUnfolder
  // copies this record into each child event with eventNumber = frame*1000
  // (candidates carry frame*1000 + candidate), so frame = eventNumber/1000
  // is the join key on both sides.
  //
  //   weights[0] = n_physics_events: physics events injected into this frame.
  //                Denominator for pileup efficiency; the flag > 0 candidate
  //                count is the numerator.
  //   weights[1] = mu_tof: measured TOF-group background (hits per cell per window).
  //   weights[2] = k_tof: TOF-group count threshold applied (capped at 999
  //                when the TOF group is disabled).
  //   weights[3] = dt, in ns: combined-group coincidence window this frame.
  //   weights[4] = dt_tof, in ns: TOF-group coincidence window this frame.
  //   weights[5] = topology:threshold.
  //   weights[6] = topology:theta_bins.
  //   weights[7] = topology:phi_bins.
  //   weights[8] = adaptive-mode flag: 1 if p_fake is in (0,1).
  //   weights[9] = p_fake.
  void emitFrameInfo(int64_t run_number, int n_physics_events, double mu_tof,
                     int k_tof, float dt, float dt_tof) {
    auto rec = m_frame_info_out()->create();
    rec.setRunNumber(static_cast<uint32_t>(run_number));
    rec.setEventNumber(0); // EventUnfolder rewrites this to frame*1000
    rec.addToWeights(double(n_physics_events));                          // [0]
    rec.addToWeights(mu_tof);                                            // [1]
    rec.addToWeights(double(std::min(k_tof, 999)));                      // [2]
    rec.addToWeights(double(dt));                                        // [3]
    rec.addToWeights(double(dt_tof));                                    // [4]
    rec.addToWeights(double(m_threshold));                               // [5]
    rec.addToWeights(double(m_theta_bins));                              // [6]
    rec.addToWeights(double(m_phi_bins));                                // [7]
    rec.addToWeights(adaptiveMode() ? 1.0 : 0.0);                        // [8]
    rec.addToWeights(double(p_fake()));                                  // [9]
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
    // (the stagger), so a cluster on a grid-1 edge is central in grid 2.
    //
    // theta is binned in cos(theta), not raw theta: equal-width cos(theta)
    // bins subtend equal solid angle, but equal-width theta bins do not,
    // since the poles subtend far less solid angle than the equator. With
    // raw-theta bins, an isotropic background would pile up in the
    // equatorial rows, and the flat mu = N/ncells share would underestimate
    // the local rate there, inflating those cells' S. Binning in cos(theta)
    // restores the equal-share assumption for an isotropic background. u
    // runs 0..ntheta over cos(theta) in [+1,-1] (theta 0..pi), so u
    // increases with theta and the bin ordering is preserved.
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
