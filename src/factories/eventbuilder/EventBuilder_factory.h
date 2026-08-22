// Copyright 2024, Jefferson Science Associates, LLC.
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// EventBuilder_factory
// -------------------
// Streaming "event building" stage. Runs once per time-frame (JANA Timeslice
// level) over the time-aligned tracker hits and *locates* candidate physics
// events. For each candidate it produces a precise event time t0 and its
// half-width dt0:
//
//   t0  = inverse-variance weighted mean of the coincident trigger-detector
//         hit times,  t0 = Sum(t_i / sigma_i^2) / Sum(1 / sigma_i^2)
//   dt0 = nsigma * sigma_t0,  sigma_t0 = 1 / sqrt(Sum 1/sigma_i^2)
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
//   weight     = dt0_ns
//   weights[0] = t0_ns   weights[1] = dt0_ns   weights[2] = phys/fake flag
//                (1 = matched MC collision, 2 = fake)

#pragma once

#include <cmath>
#include <vector>

#include "TMath.h"

#include <extensions/jana/JOmniFactory.h>

#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4eic/TrackerHitCollection.h>

struct EventBuilder_factory : public JOmniFactory<EventBuilder_factory> {

  // -- time resolutions (ns), per detector group --------------------------------
  // (slots declared before the ParameterRefs that bind to them)
  float m_res_tof      = 0.03f;   // ns (AC-LGAD ~30 ps)
  float m_res_mpgd     = 10.0f;   // ns
  float m_res_si       = 2000.0f; // ns (MAPS integration window)
  float m_nsigma       = 3.0f;
  float m_coinc_window = 50.0f;   // ns

  ParameterRef<float> timeResolution_TOF{this, "timeResolution_TOF", m_res_tof,
                                         "time resolution of TOF detectors in ns"};
  ParameterRef<float> timeResolution_MPGD{this, "timeResolution_MPGD", m_res_mpgd,
                                          "time resolution of MPGD detectors in ns"};
  ParameterRef<float> timeResolution_Silicon{this, "timeResolution_Silicon", m_res_si,
                                             "time resolution of Silicon detectors in ns"};
  ParameterRef<float> nsigma_window{this, "nsigma_window", m_nsigma,
                                    "N-sigma half-width that defines dt0 around t0"};
  ParameterRef<float> coincidence_window{
      this, "coincidence_window", m_coinc_window,
      "time window (ns) over the fast detectors used to group coincident hits"};

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

  PodioOutput<edm4hep::EventHeader> m_candidates_out{this, "EventCandidates"};

  void Configure() {}
  void ChangeRun(int32_t /*run_nr*/) {}

  // One time-aligned trigger hit, flattened across detectors.
  struct FastHit {
    double time;
    float sigma;
    Int_t th1, ph1, th2, ph2;
  };

  void Process(int64_t run_number, uint64_t /*event_number*/) {
    // m_candidates_out() is reset to a fresh empty collection by the framework.

    // -- 1. Flatten the fast (TOF+MPGD) hits, with their (theta,phi) bins ------
    std::vector<FastHit> hits;
    for (size_t det = 0; det < kNumTriggerDet && det < m_trk_in().size(); ++det) {
      const auto* coll = m_trk_in().at(det);
      if (coll == nullptr)
        continue;
      const float sigma = detTimeRes(det);
      for (const auto& hit : *coll) {
        FastHit fh{hit.getTime(), sigma, 999, 999, 999, 999};
        thetaPhiBinCalc(hit, fh.th1, fh.ph1, fh.th2, fh.ph2);
        hits.push_back(fh);
      }
    }
    std::sort(hits.begin(), hits.end(),
              [](const FastHit& a, const FastHit& b) { return a.time < b.time; });

    // -- MC collision times (for phys/fake labelling), generatorStatus == 61 ---
    std::vector<double> mc_times;
    if (m_mc_in() != nullptr) {
      for (const auto& mcp : *m_mc_in())
        if (mcp.getGeneratorStatus() == 61)
          mc_times.push_back(mcp.getTime());
    }

    // -- 2. O(N) sliding-window coincidence + (theta,phi) topology trigger -----
    const float coinc = coincidence_window();
    Int_t cnt1[12][8] = {};
    Int_t cnt2[12][8] = {};
    size_t lo = 0, hi = 0;
    while (lo < hits.size()) {
      while (hi < hits.size() && (hits[hi].time - hits[lo].time) <= coinc) {
        ++cnt1[hits[hi].th1][hits[hi].ph1];
        ++cnt2[hits[hi].th2][hits[hi].ph2];
        ++hi;
      }

      bool trigger = false;
      for (size_t it = 0; it < 12 && !trigger; ++it)
        for (size_t ip = 0; ip < 8 && !trigger; ++ip)
          if (cnt1[it][ip] > 2 || cnt2[it][ip] > 2)
            trigger = true;

      if (trigger) {
        // inverse-variance t0 over the coincident window [lo, hi)
        double sumw = 0.0, sumwt = 0.0;
        for (size_t k = lo; k < hi; ++k) {
          const double w = 1.0 / (double(hits[k].sigma) * double(hits[k].sigma));
          sumw += w;
          sumwt += w * hits[k].time;
        }
        const double t0      = sumwt / sumw;
        const double sigmat0 = std::sqrt(1.0 / sumw);
        const double dt0     = double(nsigma_window()) * sigmat0;

        emitCandidate(run_number, t0, dt0, matchPhysics(mc_times, t0));

        // consume the cluster: clear counts and jump past it
        for (size_t k = lo; k < hi; ++k) {
          --cnt1[hits[k].th1][hits[k].ph1];
          --cnt2[hits[k].th2][hits[k].ph2];
        }
        lo = hi;
      } else {
        // slide: drop the earliest hit and advance
        --cnt1[hits[lo].th1][hits[lo].ph1];
        --cnt2[hits[lo].th2][hits[lo].ph2];
        ++lo;
        if (hi < lo)
          hi = lo;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // (not const: ParameterRef::operator() is non-const)
  float detTimeRes(size_t detIdx) {
    if (detIdx < 2)
      return timeResolution_TOF();
    if (detIdx < 6)
      return timeResolution_MPGD();
    return timeResolution_Silicon();
  }

  // Returns 1 if an MC collision (status 61) sits within the window of this t0
  // and consumes it; 2 (fake) otherwise.
  int matchPhysics(std::vector<double>& mc_times, double t0) {
    const double half = coincidence_window();
    for (auto it = mc_times.begin(); it != mc_times.end(); ++it) {
      if (*it > t0 - half && *it < t0 + half) {
        mc_times.erase(it);
        return 1;
      }
    }
    return 2;
  }

  void emitCandidate(int64_t run_number, double t0, double dt0, int phys_flag) {
    auto cand = m_candidates_out()->create();
    cand.setRunNumber(static_cast<uint32_t>(run_number));
    cand.setEventNumber(m_candidates_out()->size() - 1);
    cand.setTimeStamp(static_cast<uint64_t>(std::llround(t0 * 1000.0))); // t0 in ps
    cand.setWeight(dt0);                                                 // ns
    cand.addToWeights(t0);                                               // weights[0]
    cand.addToWeights(dt0);                                              // weights[1]
    cand.addToWeights(static_cast<double>(phys_flag));                   // weights[2]
  }

  inline void thetaPhiBinCalc(edm4eic::TrackerHit hit, Int_t& thetaID1, Int_t& phiID1,
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

    thetaID1 = hitTheta / (TMath::Pi() / 12.);
    thetaID2 = (hitTheta + TMath::Pi() / 24.) / (TMath::Pi() / 12.);
    phiID1   = hitPhi / (TMath::Pi() / 8.);
    phiID2   = (hitPhi + TMath::Pi() / 16.) / (TMath::Pi() / 8.);

    // clamp into the [0,12) x [0,8) grid (staggered bins can overflow by one)
    thetaID1 = std::min(std::max(thetaID1, 0), 11);
    thetaID2 = std::min(std::max(thetaID2, 0), 11);
    phiID1   = std::min(std::max(phiID1, 0), 7);
    phiID2   = std::min(std::max(phiID2 % 8, 0), 7);
  }
};
