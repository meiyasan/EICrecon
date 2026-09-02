// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// eventbuilder_trigger -- trigger/eventbuilder performance benchmark.
//
//   eicrecon -Peventbuilder=true -Pplugins=eventbuilder,eventbuilder_trigger ...
//
// Prints a per-physics-class efficiency table every N frames and once at the
// end. See README.md in this directory for the column semantics.
//
// One tap, at PhysicsEvent level, reaching the parent frame for frame-level
// data. A Timeslice-level JEventProcessor would see every frame directly --
// including those that produced no candidate -- but JANA cannot currently
// route one: a Timeslice tap coexisting with the unfold/fold pair throws
// std::out_of_range from JEventPool::Ingest. See README.md, "Blind frames".

#include <JANA/JApplicationFwd.h>
#include <JANA/JApplication.h>
#include <JANA/Services/JParameterManager.h>

#include "EventBenchmark_processor.h"
#include "TriggerBenchmark_service.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  // RECORD_CALL_STACK has to be set HERE, at plugin load. Setting it from the
  // service's acquire_services() is too late: JANA has already built its event
  // pool by then and each JEvent's call-graph recorder is already stamped
  // disabled, so the graph comes back empty and the profile page is blank.
  bool profile = false;
  app->SetDefaultParameter(
      "eventbuilder:benchmark:profile", profile,
      "Per-factory timing breakdown, read from JANA's call graph. Implies "
      "RECORD_CALL_STACK=1, which costs throughput -- leave off for production runs.");
  if (profile)
    app->SetParameterValue("RECORD_CALL_STACK", true);

  app->ProvideService(std::make_shared<eicrecon::eb::TriggerBenchmark_service>(app));
  app->Add(new eicrecon::eb::EventBenchmark_processor());
}
}
