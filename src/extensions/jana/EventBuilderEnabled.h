// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde

#pragma once

#include <JANA/JApplication.h>

namespace eicrecon {

/// True when the EventBuilder workflow is switched on (-Peventbuilder=true,
/// alias -Peventbld=true).
///
/// The PODIO source reads its input at JEventLevel::Timeslice only when this is
/// set -- see JEventSourcePODIO_generator.h. Detector plugins have to register
/// their Timeslice-level "*Frame" factories under the same condition: register
/// them unconditionally and JANA builds a Timeslice topology for *every* job,
/// including ordinary PhysicsEvent reconstruction, which has no time-frames to
/// feed it and dies with
///
///   Could not find databundle with type_index=edm4hep::MCParticle
///   and tag=MCParticles
///
/// Safe to call from every plugin: RegisterParameter forwards to
/// SetDefaultParameter, which is idempotent for a matching default.
inline bool eventbuilderEnabled(JApplication* app) {
  const bool use_eventbuilder = app->RegisterParameter<bool>(
      "eventbuilder", false, "Read input as time-frames for the EventBuilder");
  const bool use_eventbld =
      app->RegisterParameter<bool>("eventbld", false, "Alias of eventbuilder");
  return use_eventbuilder || use_eventbld;
}

} // namespace eicrecon
