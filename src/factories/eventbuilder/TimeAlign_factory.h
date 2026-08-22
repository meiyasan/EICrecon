// Copyright 2024, Jefferson Science Associates, LLC.
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// Time-alignment factory, shared by trackers and calorimeter clusters.
//
// For each input collection it clones every hit/cluster, subtracts a straight-
// line propagation term r/c (so all times refer to a common collision t0), and
// sorts the collection by time. The logic is identical for tracker hits
// (edm4eic::TrackerHit) and calo clusters (edm4eic::Cluster), so it is written
// once as a template and instantiated once per family:
//
//   using TrkTimeAlign_factory = TimeAlign_factory<edm4eic::TrackerHit>;
//   using CalTimeAlign_factory = TimeAlign_factory<edm4eic::Cluster>;
//
// Register one instance per family in eventbuilder.cc, each with a single
// variadic input/output so the collection names distribute correctly. (The
// previous design put both families' variadic inputs on one factory, which made
// the wiring split a single name list across both and silently align only half
// the tracker collections.)

#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include <extensions/jana/JOmniFactory.h>

#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/ClusterCollection.h>

template <typename HitT> struct TimeAlign_factory : public JOmniFactory<TimeAlign_factory<HitT>> {
  using Base = JOmniFactory<TimeAlign_factory<HitT>>;

  // Collection names are supplied by the wiring (see eventbuilder.cc).
  typename Base::template VariadicPodioInput<HitT, true> m_in{this, {}};
  typename Base::template VariadicPodioOutput<HitT> m_out{this, {}};

  void Configure() {}
  void ChangeRun(int32_t /*run_nr*/) {}

  void Process(int64_t /*run_number*/, uint64_t /*event_number*/) {
    using MutT = decltype(std::declval<const HitT&>().clone());
    for (size_t i = 0; i < m_in().size(); ++i) {
      const auto* in = m_in().at(i);
      auto& out      = m_out().at(i);
      if (in == nullptr)
        continue;

      std::vector<MutT> sorted;
      sorted.reserve(in->size());
      for (const auto& hit : *in) {
        MutT copied        = hit.clone();
        const double x     = hit.getPosition()[0];
        const double y     = hit.getPosition()[1];
        const double z     = hit.getPosition()[2];
        const double r     = std::sqrt(x * x + y * y + z * z);
        copied.setTime(hit.getTime() - r * 0.0034); // subtract r/c propagation term
        sorted.push_back(copied);
      }

      std::sort(sorted.begin(), sorted.end(),
                [](const auto& a, const auto& b) { return a.getTime() < b.getTime(); });

      for (const auto& h : sorted)
        out->push_back(h);
    }
  }
};

using TrkTimeAlign_factory = TimeAlign_factory<edm4eic::TrackerHit>;
using CalTimeAlign_factory = TimeAlign_factory<edm4eic::Cluster>;
