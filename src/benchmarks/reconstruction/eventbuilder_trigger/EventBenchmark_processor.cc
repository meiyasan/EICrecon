// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "EventBenchmark_processor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <set>

#include <JANA/JApplication.h>
#include <JANA/JEvent.h>
#include <JANA/JEventSource.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>
#include <edm4eic/MCRecoTrackParticleAssociationCollection.h>
#include <utility>
#include <vector>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/MCParticleCollection.h>

#include <string>
#include <stdexcept>
#include <spdlog/logger.h>
#include "services/log/Log_service.h"

namespace eicrecon::eb {

namespace {
/// Fetch a collection, tolerating BOTH absence and an upstream factory that
/// throws while producing it.
///
/// The eventbuilder's own consumers declare these inputs optional, so when a
/// producer like TimeAlignment_factory<edm4eic::Cluster> throws (it can, on
/// real frames) they simply see nullptr and carry on. A JEventProcessor has
/// no such shielding: an exception escaping a tap is fatal to the whole job.
/// Measuring a run must never be able to kill it, so swallow the throw and
/// report no data -- the affected counters stay at zero and the table says so.
template <typename T>
const typename T::collection_type* tryGet(const JEvent& event, const std::string& name,
                                          bool& warned, spdlog::logger* log) {
  try {
    return event.GetCollection<T>(name, false);
  } catch (const std::exception& e) {
    if (!warned) {
      warned = true;
      if (log != nullptr)
        log->warn("benchmark: producer for '{}' threw ({}); those counters stay empty", name,
                  e.what());
    }
    return nullptr;
  }
}
} // namespace


namespace {

/// Majority match. An association with weight > 0.5 means more than half the
/// reconstructed object's hits came from that one MC particle -- the standard
/// ghost criterion. Any-shared-hit associations exist for nearly every track
/// and would fake ~100% purity.
constexpr double kMajority = 0.5;

/// Yellow Report calorimeter coverage and minimum detectable photon energy.
constexpr double kCaloEtaMax     = 4.0;
constexpr double kEFloorBarrel   = 0.050; // GeV, |eta| < 1
constexpr double kEFloorForward  = 0.030; // GeV, 1 <= |eta| < 4

/// The contributing sim hit chosen for one raw hit: the largest energy
/// deposit among the links pointing at it.
struct Truth {
  int                    cls = -1;
  edm4hep::SimTrackerHit sim;
  float                  edep = -1.0f;
};

bool isNeutrino(int pdg) {
  const int a = std::abs(pdg);
  return a == 12 || a == 14 || a == 16;
}

} // namespace

const std::vector<std::string>& EventBenchmark_processor::caloAssociationNames() {
  static const std::vector<std::string> names = {
      "B0ECalClusterAssociations",           "EcalBarrelClusterAssociations",
      "EcalEndcapNClusterAssociations",      "EcalEndcapPClusterAssociations",
      "HcalBarrelClusterAssociations",       "HcalEndcapPInsertClusterAssociations",
      "LFHCALClusterAssociations",           "EcalFarForwardZDCClusterAssociations",
      "HcalFarForwardZDCClusterAssociations"};
  return names;
}

EventBenchmark_processor::EventBenchmark_processor() {
  SetTypeName(NAME_OF_THIS);
  SetLevel(JEventLevel::PhysicsEvent);
  SetCallbackStyle(CallbackStyle::ExpertMode);
  EnableOrdering(true);
}

void EventBenchmark_processor::Init() {
  auto* app = GetApplication();
  m_log     = app->GetService<Log_service>()->logger("EventBuilderBenchmark");
  m_bench   = app->GetService<TriggerBenchmark_service>();
  m_bench->registerTap();

  // Mirror the eventbuilder's own acceptance so the STAGE 3 denominator is
  // the same population the trigger calls "findable".
  app->SetDefaultParameter("eventbuilder:findable:pt_min", m_pt_min);
  app->SetDefaultParameter("eventbuilder:findable:eta_abs", m_eta_abs);
}

void EventBenchmark_processor::processFrame(const JEvent& parent, std::uint64_t frame) {
  Totals                  d;
  std::map<int, ClassRow> per_class;
  std::vector<TriggerBenchmark_service::FoundKey> found_keys;

  d.frames = 1;
  

  // The prefilter's surviving list when the GNN ran, else the raw candidates.
  // Scoring the filtered list is the point: it is what actually reaches
  // reconstruction.
  const auto* cands = tryGet<edm4hep::EventHeader>(parent, "EventCandidatesFiltered",
                                                  m_warned_fetch, m_log.get());
  if (cands == nullptr)
    cands = tryGet<edm4hep::EventHeader>(parent, "EventCandidates", m_warned_fetch, m_log.get());
  if (cands == nullptr) {
    if (!m_warned_no_cands) {
      m_warned_no_cands = true;
      m_log->warn("no EventCandidates at Timeslice level -- is the eventbuilder plugin loaded "
                  "with -Peventbuilder=true? Benchmark will report zeros.");
    }
    return;
  }

  const auto* frame_info =
      tryGet<edm4hep::EventHeader>(parent, "EventBuilderFrameInfo", m_warned_fetch, m_log.get());
  const auto* scores =
      tryGet<edm4hep::EventHeader>(parent, "PrefilterScores", m_warned_fetch, m_log.get());

  // dt: the frame's own coincidence window, the same width EventBuilder_factory
  // used to group MCParticles into collisions. Reusing it keeps the injected
  // count here identical to the factory's own n_physics_events.
  double dt = 0.0;
  if (frame_info != nullptr && frame_info->size() > 0) {
    const auto w = (*frame_info)[0].getWeights();
    if (w.size() > fw::DT)
      dt = w[fw::DT];
  }

  // -- STAGE 1: per-candidate counters --------------------------------------
  bool frame_has_real = false;
  for (const auto& c : *cands) {
    const auto w = c.getWeights();
    if (w.size() <= w::CAL_PASS)
      continue; // pre-dates these weights; nothing to score

    if (w.size() > w::CAL_THRESHOLD_CFG) {
      d.trk_threshold = w[w::TRK_THRESHOLD_CFG];
      d.cal_threshold = w[w::CAL_THRESHOLD_CFG];
    }
    const bool real     = isReal(w[w::FLAG]);
    const bool trk_pass = w[w::TRK_PASS] > 0.5;
    const bool cal_pass = w[w::CAL_PASS] > 0.5;

    ++d.candidates;
    if (real) {
      ++d.real;
      frame_has_real = true;
      if (trk_pass)
        ++d.trk_ok;
      if (cal_pass)
        ++d.cal_ok;
      if (trk_pass || cal_pass)
        ++d.trig_ok;
    } else {
      ++d.fake;
      if (trk_pass)
        ++d.trk_fake;
      if (cal_pass)
        ++d.cal_fake;
      if (trk_pass || cal_pass)
        ++d.trig_fake;
    }

    // Classes this candidate is credited with (bitmask), and the collisions
    // it actually recovered (the variable-length trigger_classes tail).
    if (w.size() > w::TRIGGER_CLASSES_MASK) {
      auto mask = static_cast<unsigned long long>(w[w::TRIGGER_CLASSES_MASK]);
      for (int ci = 0; mask != 0ULL && ci < kNumClasses; ++ci, mask >>= 1)
        if ((mask & 1ULL) != 0ULL)
          ++per_class[ci].candidates;
    }
    if (real && w.size() > w::TRIGCLS) {
      const auto n = static_cast<std::size_t>(w[w::TRIGCLS]);
      std::set<int> seen;
      for (std::size_t k = 0; k < n; ++k) {
        const std::size_t it = w::TRIGCLS + 1 + 2 * k; // (time, stream) pairs
        if (it + 1 >= w.size())
          break;
        const int ci = classIndexOfStream(static_cast<int>(w[it + 1]));
        if (ci < 0)
          continue;
        found_keys.push_back({frame, ci, static_cast<double>(w[it]), trk_pass, cal_pass});
        seen.insert(ci);
      }
      for (int ci : seen) {
        ++per_class[ci].real_cands;
        // Cost of this candidate, charged to every class it carries. A
        // multi-class candidate is counted once per class on purpose: the
        // question is "what does a class cost me", and the work was done for
        // all of them together.
        if (m_event_secs > 0.0) {
          per_class[ci].proc_seconds += m_event_secs;
          ++per_class[ci].proc_cands;
        }
        // Same pass flags the totals use, charged to each class this real
        // candidate is credited with.
        if (trk_pass)
          ++per_class[ci].trk_ok;
        if (cal_pass)
          ++per_class[ci].cal_ok;
        if (trk_pass || cal_pass)
          ++per_class[ci].trig_ok;
      }
    }
  }

  // -- STAGE 2: GNN prefilter, one decision per frame ------------------------
  // "Passed" = the best non-background class outscores background, not a fixed
  // cut on the background probability alone.
  if (scores != nullptr && scores->size() > 0) {
    const auto p = (*scores)[0].getWeights();
    if (p.size() >= 2) {
      d.saw_gnn    = true;
      d.gnn_frames = 1;
      double best_other = p[1];
      for (std::size_t i = 2; i < p.size(); ++i)
        best_other = std::max<double>(best_other, p[i]);
      const bool passed = best_other > p[0];
      if (passed)
        ++d.gnn_passed;
      // A frame with zero real candidates is what the retired bkg class used
      // to be: 100% fake by construction. Accepting one is a false accept.
      if (!frame_has_real) {
        ++d.gnn_fake_frames;
        if (passed)
          ++d.gnn_fake_accepted;
      }
    }
  }

  // -- injected collisions, per class ---------------------------------------
  // Same dedup rule as EventBuilder_factory's collision_times loop: one entry
  // per (class, time) within the frame's own dt. Counted from the FRAME's
  // MCParticles, so collisions in frames that produced no candidate still
  // reach the denominator.
  const auto* mcps = tryGet<edm4hep::MCParticle>(parent, "MCParticles", m_warned_fetch, m_log.get());
  if (mcps != nullptr) {
    std::map<int, std::vector<double>> times;
    for (const auto& p : *mcps) {
      const int gs = p.getGeneratorStatus();
      if (!isStablePhysics(gs))
        continue;
      const int ci = classIndexOfStatus(gs);
      if (ci < 0)
        continue;
      auto&        v = times[ci];
      const double t = p.getTime();
      bool         merged = false;
      for (double e : v)
        if (std::abs(e - t) < dt) {
          merged = true;
          break;
        }
      if (!merged)
        v.push_back(t);
    }
    for (const auto& [ci, v] : times) {
      per_class[ci].injected += v.size();
      d.collisions_injected += v.size();
    }
  }

  m_bench->addFrame(d, per_class, found_keys);
  m_bench->maybeReport();
}


/// Accumulate one event's factory timings. The eventbuilder's own factories
/// run on the PARENT Timeslice, so reading only the child's graph (which is
/// what a PhysicsEvent-level tap sees) returns nothing at all -- the child
/// graph holds just the per-candidate reconstruction.
double EventBenchmark_processor::accumulateCallGraph(const JEvent& ev) {
  auto* cg = ev.GetJCallGraphRecorder();
  if (cg == nullptr || !cg->IsEnabled())
    return 0.0;
  double total = 0.0;
  for (const auto& n : cg->GetCallGraph()) {
    const double secs = std::chrono::duration<double>(n.end_time - n.start_time).count();
    total += secs;
    std::string label = n.callee_name;
    if (!n.callee_tag.empty())
      label = label.empty() ? n.callee_tag : label + ":" + n.callee_tag;
    if (label.empty())
      label = "(unnamed)";
    m_bench->addFactoryTime(label, secs);
  }
  return total;
}

void EventBenchmark_processor::fillResolutions(const JEvent& event, double t_ref) {
  if (!m_bench->res().enabled())
    return;

  // -- trackers: time and space, class taken from the per-hit truth label ----
  for (std::size_t d = 0; d < ResolutionHists::trackerNames().size(); ++d) {
    const std::string& det = ResolutionHists::trackerNames()[d];
    const auto* hits = tryGet<edm4eic::TrackerHit>(event, det + "RecHits", m_warned_fetch,
                                                   m_log.get());
    if (hits == nullptr || hits->size() == 0)
      continue;
    // raw-hit index -> (class, sim hit). The link weight is
    // generatorStatus + slot*1e6, so the class of the collision that actually
    // produced this hit is available without chasing pointers.
    const auto* links = tryGet<edm4eic::MCRecoTrackerHitLink>(event, det + "RawHitLinks",
                                                              m_warned_fetch, m_log.get());
    // A raw hit can have SEVERAL contributing sim hits -- charge shared across
    // neighbouring strips, or two tracks through one strip. std::map::emplace
    // does not overwrite, so keeping the first link seen picked an ARBITRARY
    // contributor: when it was not the dominant one the position residual was
    // measured against the wrong sim hit, piling up at the characteristic
    // strip separation (a spurious ~1.5 mm spike on MPGDBarrel). The link
    // weight is generatorStatus + slot*1e6 -- a label, not a fraction -- so it
    // cannot rank contributors; use the largest energy deposit, the usual
    // choice for "which particle made this hit".
    std::map<std::uint32_t, Truth> truth;
    if (links != nullptr) {
      for (const auto& l : *links) {
        const auto from = l.getFrom();
        const auto to   = l.getTo();
        if (!from.isAvailable())
          continue;
        const auto  status = static_cast<int>(std::fmod(l.getWeight(), 1.0e6));
        const auto  key    = static_cast<std::uint32_t>(from.getObjectID().index);
        const float edep   = to.isAvailable() ? to.getEDep() : -1.0f;
        auto        it     = truth.find(key);
        if (it == truth.end())
          truth.emplace(key, Truth{classIndexOfStatus(status), to, edep});
        else if (edep > it->second.edep)
          it->second = Truth{classIndexOfStatus(status), to, edep};
      }
    }
    for (const auto& h : *hits) {
      const auto raw = h.getRawHit();
      int  cls = -1;
      bool has_sim = false;
      edm4hep::SimTrackerHit sim;
      if (raw.isAvailable()) {
        auto it = truth.find(static_cast<std::uint32_t>(raw.getObjectID().index));
        if (it != truth.end()) {
          cls     = it->second.cls;
          sim     = it->second.sim;
          has_sim = sim.isAvailable();
        }
      }
      if (cls < 0)
        continue; // unlabelled hit (background, or links not written)
      if (has_sim) {
        const auto   r  = h.getPosition();
        // SENSOR time resolution. Two corrections are needed and both matter:
        //
        //  1. The rec hit is TIME-ALIGNED. TimeAlignment_factory::correct()
        //     subtracts |r|/c -- the straight-line, beta=1 propagation time --
        //     from every hit, so hits of one collision share a common time.
        //     Adding it back with the hit's own position undoes that term
        //     exactly (it is the same position the factory used).
        //  2. Compare against the SIM hit, not the collision time. The sim hit
        //     carries the true absolute arrival time, so what is left is the
        //     digitisation smear -- SiliconTrackerDigi writes
        //     t_sim + N(0, timeResolution) -- plus the factory's per-detector
        //     calibration offset, which is a constant and so shifts the peak
        //     without widening it. fitCore fits a free mean, so the width the
        //     table reports is the sensor resolution regardless of that offset.
        //
        // Measured directly from the digitiser's own debug output: the raw
        // timestamp reproduces sim_hit.getTime() to 1.6 ps median with a
        // 25 ps spread, and the child's SimHits carry those times exactly
        // (0.0000 ns). Referencing the collision time instead measures how
        // well the beta=1 straight-line alignment recovers t0 -- a real
        // quantity, but a one-sided few-ns one dominated by curvature and
        // beta<1, with no core to fit.
        constexpr double kInvC_ns_per_mm = 1.0 / 299.792458;
        const double     r_mm            = std::sqrt(static_cast<double>(r.x) * r.x +
                                                     static_cast<double>(r.y) * r.y +
                                                     static_cast<double>(r.z) * r.z);
        m_bench->res().fillTrackerTime(d, cls,
                                       static_cast<double>(h.getTime()) +
                                           r_mm * kInvC_ns_per_mm -
                                           static_cast<double>(sim.getTime()));
        const auto   sp = sim.getPosition();
        // Decompose instead of taking the 3D magnitude. The barrel readout is
        // a CylindricalGridPhiZ pinned to a FIXED radius per segmentation, so
        // every rec hit is projected onto a cylinder while the sim hit sits at
        // its true radius. Folding that fixed radial offset into one |rec-sim|
        // produced a spurious sharp spike (~1.5 mm on both barrel MPGDs) on
        // top of the real in-plane resolution.
        const double r_rec = std::hypot(static_cast<double>(r.x), static_cast<double>(r.y));
        const double r_sim = std::hypot(static_cast<double>(sp.x), static_cast<double>(sp.y));
        m_bench->res().fillTrackerSpaceRadial(d, cls, r_rec - r_sim);

        double dphi = std::atan2(static_cast<double>(r.y), static_cast<double>(r.x)) -
                      std::atan2(static_cast<double>(sp.y), static_cast<double>(sp.x));
        while (dphi > M_PI)
          dphi -= 2 * M_PI;
        while (dphi < -M_PI)
          dphi += 2 * M_PI;
        const double rphi = r_sim * dphi;                                  // arc length
        const double dz   = static_cast<double>(r.z) - static_cast<double>(sp.z);
        m_bench->res().fillTrackerSpace(d, cls, std::hypot(rphi, dz));
      }
    }
  }


  // -- TOF charge-sharing chain -------------------------------------------
  // What central tracking consumes: SiliconChargeSharing splits each deposit
  // across neighbouring AC-LGAD pads, LGADHitClustering turns those into the
  // Measurement2D below. Truth therefore lives in *SharedRawHitLinks and
  // *SharedSimHits, NOT the plain *RawHitLinks used above -- tracking.cc's
  // CentralTrackingRawHitLinks collector makes the same distinction.
  {
    const std::vector<std::string> meas = {"TOFBarrelClusterHits", "TOFEndcapClusterHits"};
    const std::vector<std::string> lnk  = {"TOFBarrelSharedRawHitLinks",
                                           "TOFEndcapSharedRawHitLinks"};
    for (std::size_t i = 0; i < meas.size(); ++i) {
      const auto* ms = tryGet<edm4eic::Measurement2D>(event, meas[i], m_warned_fetch, m_log.get());
      if (ms == nullptr || ms->size() == 0)
        continue;
      const auto* links =
          tryGet<edm4eic::MCRecoTrackerHitLink>(event, lnk[i], m_warned_fetch, m_log.get());
      std::map<std::uint32_t, Truth> truth;
      if (links != nullptr) {
        for (const auto& l : *links) {
          const auto from = l.getFrom();
          if (!from.isAvailable())
            continue;
          const auto  status = static_cast<int>(std::fmod(l.getWeight(), 1.0e6));
          const auto  to     = l.getTo();
          const auto  key    = static_cast<std::uint32_t>(from.getObjectID().index);
          const float edep   = to.isAvailable() ? to.getEDep() : -1.0f;
          auto        itr    = truth.find(key);
          if (itr == truth.end())
            truth.emplace(key, Truth{classIndexOfStatus(status), to, edep});
          else if (edep > itr->second.edep)
            itr->second = Truth{classIndexOfStatus(status), to, edep};
        }
      }
      const std::size_t det = 10 + i; // TOFBarrelShared / TOFEndcapShared
      for (const auto& m : *ms) {
        // A Measurement2D is a cluster of pads; score its constituent hits so
        // the residual is per hit, comparable with the other nine trackers.
        for (const auto& h : m.getHits()) {
          if (!h.isAvailable())
            continue;
          const auto raw = h.getRawHit();
          if (!raw.isAvailable())
            continue;
          auto it = truth.find(static_cast<std::uint32_t>(raw.getObjectID().index));
          if (it == truth.end() || it->second.cls < 0)
            continue;
          const int cls = it->second.cls;
          const auto sim = it->second.sim;
          if (sim.isAvailable()) {
            const auto r  = h.getPosition();
            // Sensor resolution, undoing the r/c alignment as at the first
            // fill site above.
            constexpr double kInvC_ns_per_mm = 1.0 / 299.792458;
            const double     r_mm            = std::sqrt(static_cast<double>(r.x) * r.x +
                                                         static_cast<double>(r.y) * r.y +
                                                         static_cast<double>(r.z) * r.z);
            m_bench->res().fillTrackerTime(det, cls,
                                           static_cast<double>(h.getTime()) +
                                               r_mm * kInvC_ns_per_mm -
                                               static_cast<double>(sim.getTime()));
            const auto sp = sim.getPosition();
            const double r_rec = std::hypot(static_cast<double>(r.x), static_cast<double>(r.y));
            const double r_sim = std::hypot(static_cast<double>(sp.x), static_cast<double>(sp.y));
            m_bench->res().fillTrackerSpaceRadial(det, cls, r_rec - r_sim);
            double dphi = std::atan2(static_cast<double>(r.y), static_cast<double>(r.x)) -
                          std::atan2(static_cast<double>(sp.y), static_cast<double>(sp.x));
            while (dphi > M_PI)
              dphi -= 2 * M_PI;
            while (dphi < -M_PI)
              dphi += 2 * M_PI;
            m_bench->res().fillTrackerSpace(
                det, cls,
                std::hypot(r_sim * dphi,
                           static_cast<double>(r.z) - static_cast<double>(sp.z)));
          }
        }
      }
    }
  }

  // -- calo: time, angle and energy against the matched MC particle ----------
  for (std::size_t sysi = 0; sysi < caloAssociationNames().size(); ++sysi) {
    const auto* ca = tryGet<edm4eic::MCRecoClusterParticleAssociation>(
        event, caloAssociationNames()[sysi], m_warned_fetch, m_log.get());
    if (ca == nullptr)
      continue;
    for (const auto& a : *ca) {
      if (a.getWeight() <= kMajority)
        continue;
      const auto clu = a.getRec();
      const auto mc  = a.getSim();
      if (!clu.isAvailable() || !mc.isAvailable())
        continue;
      const int cls = classIndexOfStatus(mc.getGeneratorStatus());
      if (cls < 0)
        continue;
      m_bench->res().fillCaloTime(sysi, cls, static_cast<double>(clu.getTime()) - t_ref);

      const auto mom = mc.getMomentum();
      const double pt = std::hypot(static_cast<double>(mom.x), static_cast<double>(mom.y));
      const double mc_eta = std::asinh(static_cast<double>(mom.z) / std::max(pt, 1e-9));
      const double mc_phi = std::atan2(static_cast<double>(mom.y), static_cast<double>(mom.x));
      const auto   cp     = clu.getPosition();
      const double cpt    = std::hypot(static_cast<double>(cp.x), static_cast<double>(cp.y));
      const double c_eta  = std::asinh(static_cast<double>(cp.z) / std::max(cpt, 1e-9));
      const double c_phi  = std::atan2(static_cast<double>(cp.y), static_cast<double>(cp.x));
      double dphi = std::fmod(c_phi - mc_phi + M_PI, 2 * M_PI);
      if (dphi < 0)
        dphi += 2 * M_PI;
      dphi -= M_PI;
      m_bench->res().fillCaloAngle(sysi, cls, std::hypot(c_eta - mc_eta, dphi));

      const double e_mc = std::hypot(pt, static_cast<double>(mom.z));
      if (e_mc > 0.0)
        m_bench->res().fillCaloEnergy(sysi, cls, (static_cast<double>(clu.getEnergy()) - e_mc) / e_mc);
    }
  }
}

void EventBenchmark_processor::ProcessSequential(const JEvent& event) {
  m_bench->tickClock();

  // Per-factory timing, straight from JANA's call graph: every node carries
  // the start/end timestamps of one factory call. Summing per callee gives the
  // "what is slow" breakdown without instrumenting any factory by hand.
  // Child graph here; the parent frame's graph is read once per frame in
  // processFrame(), so a frame's factories are not counted once per child.
  m_event_secs = m_bench->profiling() ? accumulateCallGraph(event) : 0.0;

  const bool has_parent = event.HasParent(JEventLevel::Timeslice);
  {
    // Every event, not just the first: a watch-mode run rotates through files
    // and the source changes underneath us.
    const JEvent* se = has_parent ? &event.GetParent(JEventLevel::Timeslice) : &event;
    if (se->GetJEventSource() != nullptr)
      m_bench->addInputFile(se->GetJEventSource()->GetResourceName());
  }
  if (!m_checked_parent) {
    m_checked_parent = true;
    // Name the dataset in the report rather than "(unknown)". A child
    // PhysicsEvent has NO source of its own -- the unfolder made it, and the
    // file belongs to the parent Timeslice -- so ask the parent first.
    const JEvent* src_ev = has_parent ? &event.GetParent(JEventLevel::Timeslice) : &event;
    if (src_ev->GetJEventSource() != nullptr)
      m_bench->addInputFile(src_ev->GetJEventSource()->GetResourceName());
    if (!has_parent) {
      m_bench->markStandalone();
      m_log->info("no Timeslice parent -- standalone mode: scoring written child events only.");
    }
  }

  // Frame-level work runs ONCE per parent frame, on its first child. Children
  // of one frame arrive together, so tracking the last frame seen is enough.
  if (has_parent) {
    const auto& parent = event.GetParent(JEventLevel::Timeslice);
    const auto  frame  = parent.GetEventNumber();
    // Remember every frame seen, not just the last one. Children of different
    // frames interleave across worker threads -- ordering is enforced within a
    // level, not between a frame and the children of the frame before it -- so
    // `frame != m_last_frame` both re-admitted a frame whose children arrived
    // in two bursts (processing it TWICE and double-counting its STAGE 3 hits)
    // and read every backwards step as a gap.
    //
    // That gap heuristic was where BLIND came from, and it was pure noise:
    // measured on one 40-frame sample, 20 and 26 blind frames on two
    // 4-thread runs of identical input, and 0 on two single-threaded runs.
    // Zero is the truth -- `children` equalled `candidates` in all four, so
    // every child was always processed. BLIND is now the frame tap's count
    // minus the distinct frames seen here, which is order-independent.
    if (m_seen_frames.insert(frame).second)
      if (m_bench->profiling())
        accumulateCallGraph(parent); // books each factory; the total is not a row
      processFrame(parent, frame);
  }
  if (!m_bench->stage3())
    return;

  const auto* mcps = tryGet<edm4hep::MCParticle>(event, "MCParticles", m_warned_fetch, m_log.get());
  if (mcps == nullptr)
    return;

  // The collisions THIS candidate is credited with, from its trigger_classes
  // tail. The child's MCParticles are gated to +-max(mc_time_window, dt), so
  // they also contain particles from neighbouring collisions the candidate was
  // never credited with -- a different physics class entirely, in general.
  // Counting those in the denominator understates the efficiency badly (a
  // factor ~2 when this filter was first added). Restrict to the credited
  // collisions so numerator and denominator describe the same population; the
  // track-matching loop below applies the same restriction, so a track built
  // from an uncredited collision's particle is dropped from BOTH sides.
  const auto* info0 =
      tryGet<edm4hep::EventHeader>(event, "EventBuilderInfo", m_warned_fetch, m_log.get());
  std::vector<std::pair<int, double>> collisions; // (class_index, time)
  double match_tol = 5.0;                          // ns, fallback
  if (info0 != nullptr && info0->size() > 0) {
    const auto wv = (*info0)[0].getWeights();
    if (wv.size() > w::TRIGCLS) {
      const auto n = static_cast<std::size_t>(wv[w::TRIGCLS]);
      for (std::size_t k = 0; k < n; ++k) {
        const std::size_t it = w::TRIGCLS + 1 + 2 * k;
        if (it + 1 >= wv.size())
          break;
        const int ci = classIndexOfStream(static_cast<int>(wv[it + 1]));
        if (ci >= 0)
          collisions.emplace_back(ci, static_cast<double>(wv[it]));
      }
    }
  }
  // Tolerance: the frame's own coincidence window, the same width the trigger
  // used to decide a particle's collision belongs to this candidate.
  if (const auto* fi =
          tryGet<edm4hep::EventHeader>(event, "EventBuilderFrameInfo", m_warned_fetch, m_log.get());
      fi != nullptr && fi->size() > 0) {
    const auto fw = (*fi)[0].getWeights();
    if (fw.size() > fw::DT && fw[fw::DT] > 0.0)
      match_tol = fw[fw::DT];
  }
  /// Does this particle belong to a collision the candidate is credited with?
  const auto belongs = [&](int ci, double t) {
    if (collisions.empty())
      return true; // no tail stored: fall back to the whole child
    for (const auto& [cci, ct] : collisions)
      if (cci == ci && std::abs(ct - t) <= match_tol)
        return true;
    return false;
  };

  Totals                  d;
  std::map<int, ClassRow> per_class;

  // -- truth populations, indexed by MCParticle index ------------------------
  std::map<std::size_t, int> charged_class; // in-acceptance stable charged -> class
  std::map<std::size_t, int> neutral_class; // in-acceptance stable neutral -> class
  std::vector<std::size_t>   neutral_idx;
  std::vector<double>        neutral_eta, neutral_phi;

  for (std::size_t i = 0; i < mcps->size(); ++i) {
    const auto p  = (*mcps)[i];
    const int  gs = p.getGeneratorStatus();
    if (!isStable(gs))
      continue;
    const auto   mom = p.getMomentum();
    const double pt  = std::hypot(static_cast<double>(mom.x), static_cast<double>(mom.y));
    const double eta = std::asinh(static_cast<double>(mom.z) / std::max(pt, 1e-9));
    const int    ci  = classIndexOfStatus(gs); // -1 for native primaries

    if (std::abs(p.getCharge()) > 0.01) {
      if (pt > m_pt_min && std::abs(eta) < m_eta_abs && belongs(ci, p.getTime()))
        charged_class[i] = ci;
    } else if (!isNeutrino(p.getPDG())) {
      // Energy taken as |p|: exact for photons, a slight underestimate for
      // massive neutrals whose rest mass alone clears these floors anyway.
      const double e     = std::hypot(pt, static_cast<double>(mom.z));
      const double floor = (std::abs(eta) < 1.0) ? kEFloorBarrel : kEFloorForward;
      if (std::abs(eta) < kCaloEtaMax && e > floor && belongs(ci, p.getTime())) {
        neutral_class[i] = ci;
        neutral_idx.push_back(i);
        neutral_eta.push_back(eta);
        neutral_phi.push_back(std::atan2(static_cast<double>(mom.y), static_cast<double>(mom.x)));
      }
    }
  }

  for (const auto& [idx, ci] : charged_class) {
    ++d.trk_expected;
    if (ci >= 0)
      ++per_class[ci].trk_expected;
  }

  // The eventbuilder also counted in-acceptance stable charged MC when it
  // built this candidate -- weights[N_EXPECTED]. It is recorded here as a
  // cross-check only; the DENOMINATOR is the loop above. See below for why
  // the two legitimately differ.
  const auto* info = tryGet<edm4hep::EventHeader>(event, "EventBuilderInfo", m_warned_fetch,
                                                  m_log.get());
  if (info != nullptr && info->size() > 0) {
    const auto w = (*info)[0].getWeights();
    if (w.size() > w::N_EXPECTED && w[w::N_EXPECTED] > 0.0) {
      // Recorded ALONGSIDE the recomputed count, never instead of it, and the
      // two are EXPECTED to differ. Both apply the same acceptance -- they
      // read the same JANA keys, eventbuilder:findable:{pt_min,eta_abs}, so
      // they cannot drift apart on cuts. What differs is how many collisions
      // each covers:
      //
      //   stored weights[N_EXPECTED]  the EARLIEST matched collision only
      //                               (EventBuilder_factory.h: nexps[i] =
      //                               earliest_nexp)
      //   recomputed trk_expected     EVERY collision this candidate is
      //                               credited with, i.e. its whole
      //                               trigger_classes tail
      //
      // So a candidate credited with two overlapping collisions contributes
      // both to the recomputed count and only the earlier one to the stored
      // count, and the recomputed value is the larger whenever pile-up is
      // credited. Measured on an 80-frame unified-mix sample: 3612 recomputed
      // against 2718 stored, a ratio of 1.33 -- pile-up, not a bug.
      //
      // The recomputed value is the one the matched-track numerator is drawn
      // from, so it is the only valid denominator here. Note what this is NOT:
      // it does not sweep in particles from UNCREDITED collisions. The
      // belongs() filter above drops those, from numerator and denominator
      // alike -- verified on a child where an upsilon collision overlapped a
      // credited dempq10to20 one and CKF reconstructed one of the upsilon
      // particles: that track was excluded from both sides, leaving 2/2.
      d.stored_nexp += static_cast<std::uint64_t>(w[w::N_EXPECTED]);
      d.used_stored_nexp = true;
    }
  }

  // -- STAGE 3 tracking ------------------------------------------------------
  const auto* trk_assoc = tryGet<edm4eic::MCRecoTrackParticleAssociation>(
      event, "CentralCKFTrackAssociations", m_warned_fetch, m_log.get());
  const auto* trks = tryGet<edm4eic::Track>(event, "CentralCKFTracks", m_warned_fetch, m_log.get());
  if (trks != nullptr)
    d.trk_reco += trks->size();
  if (trk_assoc != nullptr) {
    d.saw_stage3 = true;
    std::set<std::size_t> matched;
    std::set<std::size_t> matched_rec; // distinct RECONSTRUCTED tracks with a
                                       // majority match -> purity numerator
    // Every reconstructed track touching this class's truth at ALL, majority
    // or not: the per-class purity denominator. The shortfall against
    // matched is share-contaminated reconstruction, which unlike a ghost IS
    // attributable to a class.
    std::map<int, std::set<std::size_t>> reco_any, reco_good;
    for (const auto& a : *trk_assoc) {
      if (a.getWeight() <= kMajority)
        continue;
      {
        // Per-class purity, in tracks: every track with a majority owner of
        // this class, and the subset whose owner is credited + in acceptance.
        const auto sim_any = a.getSim();
        const auto rec_any = a.getRec();
        if (sim_any.isAvailable() && rec_any.isAvailable()) {
          const int ci_any = classIndexOfStatus(sim_any.getGeneratorStatus());
          if (ci_any >= 0) {
            const auto ridx = static_cast<std::size_t>(rec_any.getObjectID().index);
            reco_any[ci_any].insert(ridx);
            if (charged_class.count(static_cast<std::size_t>(sim_any.getObjectID().index)) != 0)
              reco_good[ci_any].insert(ridx);
          }
        }
      }
      const auto sim = a.getSim();
      if (!sim.isAvailable())
        continue;
      const auto idx = static_cast<std::size_t>(sim.getObjectID().index);
      matched.insert(idx);
      const auto rec = a.getRec();
      if (rec.isAvailable()) {
        // PURITY numerator: the track's majority owner is truth this candidate
        // is credited with AND in acceptance. Counting merely "has a majority
        // owner" is vacuous -- every track in a gated child has one, so that
        // ratio is 100% by construction and measures nothing. The shortfall
        // here is real contamination: tracks owned by another collision, by a
        // secondary, or by a particle outside acceptance.
        if (charged_class.find(idx) != charged_class.end())
          matched_rec.insert(static_cast<std::size_t>(rec.getObjectID().index));
      }
    }
    d.trk_reco_matched += matched_rec.size();
    for (const auto& [ci_t, set_t] : reco_any)
      per_class[ci_t].trk_reco_any += set_t.size();
    for (const auto& [ci_t, set_t] : reco_good)
      per_class[ci_t].trk_reco_good += set_t.size();
    for (std::size_t idx : matched) {
      auto it = charged_class.find(idx);
      if (it == charged_class.end())
        continue; // shower secondary or out of acceptance: not a denominator entry
      ++d.trk_matched;
      if (it->second >= 0)
        ++per_class[it->second].trk_matched;
    }
  }

  // -- STAGE 3 calorimetry ---------------------------------------------------
  // Truth targets are angularly GROUPED: two photons merged into one shower
  // cannot be resolved independently by any clustering algorithm, so counting
  // them separately would cap efficiency below 100% for perfect reco.
  const double dr    = m_bench->caloDr();
  const auto   n_neu = neutral_idx.size();
  std::vector<std::size_t> parent(n_neu);
  for (std::size_t i = 0; i < n_neu; ++i)
    parent[i] = i;
  const std::function<std::size_t(std::size_t)> find = [&](std::size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x         = parent[x];
    }
    return x;
  };
  for (std::size_t a = 0; a < n_neu; ++a)
    for (std::size_t b = a + 1; b < n_neu; ++b) {
      const double deta = neutral_eta[a] - neutral_eta[b];
      double       dphi = std::fmod(neutral_phi[a] - neutral_phi[b] + M_PI, 2 * M_PI);
      if (dphi < 0)
        dphi += 2 * M_PI;
      dphi -= M_PI;
      if (deta * deta + dphi * dphi < dr * dr)
        parent[find(a)] = find(b);
    }
  std::map<std::size_t, std::set<std::size_t>> groups; // root -> member indices
  for (std::size_t i = 0; i < n_neu; ++i)
    groups[find(i)].insert(neutral_idx[i]);

  std::set<std::size_t> matched_neutral;
  for (const auto& name : caloAssociationNames()) {
    // Cluster collection behind this association, for the purity denominator.
    std::string cl = name;
    cl.replace(cl.find("ClusterAssociations"), std::string("ClusterAssociations").size(),
               "Clusters");
    const auto* cls = tryGet<edm4eic::Cluster>(event, cl, m_warned_fetch, m_log.get());
    if (cls != nullptr)
      d.cal_reco += cls->size();
    const auto* ca = tryGet<edm4eic::MCRecoClusterParticleAssociation>(event, name,
                                                                      m_warned_fetch, m_log.get());
    if (ca == nullptr)
      continue;
    d.saw_stage3 = true;
    std::set<std::size_t> matched_rec_cl;
    std::map<int, std::set<std::size_t>> cl_any, cl_good;
    for (const auto& a : *ca) {
      if (a.getWeight() <= kMajority)
        continue;
      {
        const auto sim_any = a.getSim();
        const auto rec_any = a.getRec();
        if (sim_any.isAvailable() && rec_any.isAvailable()) {
          const int ci_any = classIndexOfStatus(sim_any.getGeneratorStatus());
          if (ci_any >= 0) {
            const auto ridx = static_cast<std::size_t>(rec_any.getObjectID().index);
            cl_any[ci_any].insert(ridx);
            if (neutral_class.count(static_cast<std::size_t>(sim_any.getObjectID().index)) != 0)
              cl_good[ci_any].insert(ridx);
          }
        }
      }
      const auto sim = a.getSim();
      if (!sim.isAvailable())
        continue;
      const auto idx_s = static_cast<std::size_t>(sim.getObjectID().index);
      matched_neutral.insert(idx_s);
      const auto rec = a.getRec();
      if (rec.isAvailable() && neutral_class.find(idx_s) != neutral_class.end())
        matched_rec_cl.insert(static_cast<std::size_t>(rec.getObjectID().index));
    }
    d.cal_reco_matched += matched_rec_cl.size();
    for (const auto& [ci_t, set_t] : cl_any)
      per_class[ci_t].cal_reco_any += set_t.size();
    for (const auto& [ci_t, set_t] : cl_good)
      per_class[ci_t].cal_reco_good += set_t.size();
  }
  for (const auto& [root, members] : groups) {
    // The group's class is its first member's; merged photons come from the
    // same shower and therefore the same collision.
    int ci = -1;
    for (std::size_t m : members) {
      auto it = neutral_class.find(m);
      if (it != neutral_class.end()) {
        ci = it->second;
        break;
      }
    }
    ++d.cal_expected;
    if (ci >= 0)
      ++per_class[ci].cal_expected;
    bool found = false;
    for (std::size_t m : members)
      if (matched_neutral.count(m) != 0) {
        found = true;
        break;
      }
    if (found) {
      ++d.cal_matched;
      if (ci >= 0)
        ++per_class[ci].cal_matched;
    }
  }

  // Reference every residual to the TRUE collision time (weights[MC_T]), not
  // to t0. t0 is the inverse-variance weighted mean of the very TOF/MPGD hits
  // being measured, so a t0-referenced residual is correlated with its own
  // reference and collapses toward zero whenever one precise hit dominates the
  // weighting. Fakes matched no collision and so have no truth time: they are
  // skipped rather than referenced to a meaningless mc_t of 0.
  if (info != nullptr && info->size() > 0) {
    const auto wv = (*info)[0].getWeights();
    if (wv.size() > w::MC_T && isReal(wv[w::FLAG])) {
      fillResolutions(event, wv[w::MC_T]);

      // Candidate-level timing: how well the trigger's own t0 recovers the
      // true collision time, and whether the uncertainty it quotes for that
      // t0 matches the spread actually observed.
      const double dt  = wv[w::T0] - wv[w::MC_T];
      // weights[T0SIGMA] is an N-sigma half-width, not a 1-sigma error (see
      // the service). Convert before using it as a pull denominator.
      const double nsig  = std::max(1.0, m_bench->nsigmaWindow());
      const double sig   = (wv.size() > w::T0SIGMA) ? wv[w::T0SIGMA] / nsig : 0.0;
      // Charged to every class the candidate is credited with, the same
      // attribution the per-class efficiency rows use.
      auto mask = (wv.size() > w::TRIGGER_CLASSES_MASK)
                      ? static_cast<unsigned long long>(wv[w::TRIGGER_CLASSES_MASK])
                      : 0ULL;
      for (int ci = 0; mask != 0ULL && ci < kNumClasses; ++ci, mask >>= 1) {
        if ((mask & 1ULL) == 0ULL)
          continue;
        m_bench->res().fillT0Residual(ci, dt);
        if (sig > 0.0) {
          m_bench->res().fillT0Sigma(ci, sig);
          m_bench->res().fillT0Pull(ci, dt / sig);
        }
      }
    }
  }

  d.children = 1;
  m_bench->addEvent(d, per_class);
}

void EventBenchmark_processor::Finish() {
  // Report the distinct frames this tap reached. The service subtracts it
  // from the frame tap's own count to get the genuinely blind frames -- the
  // ones that produced no child at all.
  m_bench->addChildFrames(m_seen_frames.size());
  m_bench->tapFinished();
}

} // namespace eicrecon::eb
