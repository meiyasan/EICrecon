// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Wouter Deconinck

#pragma once

#include <string>
#include <vector>

namespace eicrecon {

struct CalorimeterHitRecoConfig {

  // digitization settings
  unsigned int capADC{1};
  double dyRangeADC{1};
  unsigned int pedMeanADC{0};
  double pedSigmaADC{0};
  double resolutionTDC{1};
  double corrMeanScale{1};

  // Hit timing resolution model: sigma_t(E) = sqrt((timeErrorScale/sqrt(E[GeV]))^2 +
  // timeErrorOffset^2) — the standard stochastic (photostatistics) + constant
  // (light-collection spread, TDC quantization, electronics jitter) calorimeter
  // timing form. Both default to 0 (timeError = 0), matching prior behavior for
  // any detector not yet measured/configured (see each detector plugin for
  // per-detector values fit from truth-matched hit residuals, 2026-07).
  double timeErrorScale{0}; // stochastic term [ns * sqrt(GeV)]
  double timeErrorOffset{0}; // constant term [ns]

  // zero suppression
  double thresholdFactor{0};
  double thresholdValue{0};

  // sampling fraction
  std::string sampFrac{"1.0"};

  // readout fields
  std::string readout{""};
  std::string layerField{""};
  std::string sectorField{""};

  // name of detelment or fields to find the local detector (for global->local transform)
  // if nothing is provided, the lowest level DetElement (from cellID) will be used
  std::string localDetElement{""};
  std::vector<std::string> localDetFields{};
  std::string maskPos{""};
  std::vector<std::string> maskPosFields{};
};

} // namespace eicrecon
