// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// EventBenchmark_processor -- the PhysicsEvent-level half of the benchmark.
//
// Runs on each unfolded candidate, where the full reconstruction is available,
// and scores STAGE 3: ACTS tracks and calorimeter clusters against MC truth,
// split per physics class by the truth particle's own generatorStatus band.
//
// Also detects standalone mode. Replaying a written *.eicrecon.root gives
// PhysicsEvents with no Timeslice parent, so the frame tap never runs; the
// benchmark then reports on the children alone and says so in the header.

#pragma once

#include <cstdint>
#include <set>
#include <memory>
#include <string>
#include <vector>

#include <JANA/JEventProcessor.h>
#include <spdlog/fwd.h>

#include "TriggerBenchmark_service.h"

namespace eicrecon::eb {

class EventBenchmark_processor : public JEventProcessor {
public:
  EventBenchmark_processor();

  void Init() override;
  void ProcessSequential(const JEvent& event) override;
  void Finish() override;

private:
  std::shared_ptr<TriggerBenchmark_service> m_bench;
  std::shared_ptr<spdlog::logger>           m_log;

  double m_pt_min  = 0.1; ///< eventbuilder:findable:pt_min  (GeV/c)
  double m_eta_abs = 3.5; ///< eventbuilder:findable:eta_abs
  bool   m_checked_parent = false;
  /// Distinct parent frames whose frame-level work has been done. A set, not
  /// a last-seen value: children arrive interleaved across threads.
  std::set<std::uint64_t> m_seen_frames;
  bool          m_warned_no_cands = false;

  /// Frame-level counters (STAGE 1, STAGE 2, injected collisions), taken from
  /// the parent Timeslice on the first child of each frame.
  void processFrame(const JEvent& parent, std::uint64_t frame);

  /// Per-detector, per-class time/space/energy residuals for this candidate.
  /// t0 is the candidate's own weights[T0]; the class comes from the per-hit
  /// RawHitLink label on the tracker side and from the matched MC particle on
  /// the calo side.
  void fillResolutions(const JEvent& event, double t_ref);

  bool m_warned_fetch = false;

  /// The 9 wired calorimeter association collections. EcalLumiSpec has no
  /// plugin yet, so it is absent here rather than resolving to nullptr.
  static const std::vector<std::string>& caloAssociationNames();
};

} // namespace eicrecon::eb
