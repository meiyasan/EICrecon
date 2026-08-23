// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 - 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// TruthClassLabels.h
// ------------------
// THE canonical eventbuilder truth conventions: the physics-class table, the
// generatorStatus decode, and the EventBuilderInfo / EventBuilderFrameInfo
// weight indices.
//
// Everything that encodes or decodes these lived in four hand-kept copies
// (EventBuilder_factory.h, EventPrefilter_factory.h, EventBuilderTest.cpp,
// and the external analysis tooling). The rules in this directory's
// README.md require any change to be made in lockstep across all of them;
// this header exists so that in-repo "lockstep" is a single edit. See that
// README for the prose.
//
// Header-only and dependency-free on purpose: the pure-C++ unit tests link
// it without JANA, PODIO, or ROOT.

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

namespace eicrecon::eb {

// -- physics classes ---------------------------------------------------------
// Order is load-bearing: the index IS the class_index encoded in
// generatorStatus and in trigger_classes_mask. Never reorder; append only.
inline constexpr std::array<std::string_view, 26> kClassNames = {
    "ncdisq1",     "ncdisq10",    "ncdisq100",  "ncdisq1000", "ccdisq100",
    "ccdisq1000",  "ddis",        "dvcs",       "ddvcs",      "dvmp",
    "dempq3to10",  "dempq10to20", "dempq20to35", "tcs",       "jpsi",
    "jpsiphoto",   "urho",        "cohrho",     "upi0",       "upsilon",
    "mesonsf",     "cohphi",      "rho",        "spectroscopy", "omega",
    "photoprod"};

inline constexpr int kNumClasses = static_cast<int>(kClassNames.size());

// Physics classes occupy 1000-wide bands starting at 10000 (LIVE scheme):
//   base = kPhysicsStatusBase + class_index * kClassBandWidth
// so class 0 is 10000-10999 and class 25 ends at 35999, inside the 16-bit
// ceiling DDG4 imposes on generatorStatus. Machine background sits below,
// in 2000-6999. A wider 10000-wide layout is designed but NOT live; do not
// decode against it (see README.md).
inline constexpr int kPhysicsStatusBase = 10000;
inline constexpr int kClassBandWidth    = 1000;

/// A stable final-state particle of a physics collision.
constexpr bool isStablePhysics(int gs) {
  return gs >= kPhysicsStatusBase && gs % kClassBandWidth == 1;
}

/// Stable in the broad sense: native status 1, or a band-encoded stable
/// primary. Covers legacy trees (bands from 7000) as well as live ones.
constexpr bool isStable(int gs) { return gs == 1 || (gs >= 7000 && gs % kClassBandWidth == 1); }

/// class_index from a full generatorStatus, or -1 when it carries no class
/// (native primaries, Geant4 secondaries, machine background).
constexpr int classIndexOfStatus(int gs) {
  if (gs < kPhysicsStatusBase)
    return -1;
  const int ci = (gs - kPhysicsStatusBase) / kClassBandWidth;
  return ci < kNumClasses ? ci : -1;
}

/// class_index from a `stream` (= generatorStatus / 1000, as stored in the
/// trigger_classes tail). 0 = native primary, >= 10 = physics class.
constexpr int classIndexOfStream(int stream) {
  const int ci = stream - kPhysicsStatusBase / kClassBandWidth;
  return (ci >= 0 && ci < kNumClasses) ? ci : -1;
}

/// Display name for a class index; "?" when out of range so a corrupt
/// status prints instead of indexing out of bounds.
constexpr std::string_view className(int class_index) {
  return (class_index >= 0 && class_index < kNumClasses)
             ? kClassNames[static_cast<std::size_t>(class_index)]
             : "?";
}

// -- EventBuilderInfo weights (one row per candidate) ------------------------
// Index 25 (TRIGCLS) starts a VARIABLE-LENGTH tail: a count, then that many
// (time, stream) pairs. Nothing may ever be appended after it.
namespace w {
enum Index : std::size_t {
  T0 = 0,
  T0SIGMA,
  FLAG,     ///< COUNT of real MC collisions in the window; 0 = fake, >0 = real
  MC_T,     ///< time of the EARLIEST matched collision
  S,
  S_TOF,
  P_TAIL,
  CMAX,
  CMAX_TOF,
  E_CALO,
  N_TRACKLETS,
  N_EXPECTED,
  TRK_THRESHOLD_CFG,
  TRK_PASS,
  CAL_THRESHOLD_CFG,
  CAL_PASS,
  TRIGGER_BITS,
  TRIGGER_CLASSES_MASK, ///< bit N set = class_index N had a collision in window
  MU,
  K,
  WIN_THETA,
  WIN_PHI,
  WIN_GRID,
  S_EMPIRICAL,
  FANO,
  TRIGCLS ///< variable-length tail; MUST stay last
};
} // namespace w

// -- EventBuilderFrameInfo weights (one row per frame) -----------------------
// eventNumber = frame*1000; join to candidates via eventNumber / 1000.
namespace fw {
enum Index : std::size_t {
  N_PHYSICS_EVENTS = 0,
  MU_TOF,
  K_TOF,
  DT,
  DT_TOF,
  TOPOLOGY_THRESHOLD,
  THETA_BINS,
  PHI_BINS,
  ADAPTIVE_FLAG,
  P_FAKE
};
} // namespace fw

/// A candidate is "real" when it matched at least one injected collision.
constexpr bool isReal(double flag_weight) { return flag_weight > 0.0; }

} // namespace eicrecon::eb
