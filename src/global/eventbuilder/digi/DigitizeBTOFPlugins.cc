// Copyright 2022, Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024, Dmitry Kalinkin

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
#include <edm4eic/EDM4eicVersion.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <TMath.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/SimTrackerHit.h>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "algorithms/digi/SiliconChargeSharingConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/CFDROCDigitization_factory.h"
#include "factories/digi/PulseCombiner_factory.h"
#include "factories/digi/PulseGeneration_factory.h"
#include "factories/digi/SiliconChargeSharing_factory.h"
#include "factories/digi/SiliconPulseDiscretization_factory.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/reco/LGADHitCalibration_factory.h"
#include "factories/tracking/LGADHitClustering_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

// extern "C" {
void InitPlugin_digiBTOF(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // Digitization
  app->Add((new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "TOFBarrelRawHitDigi", {"EventHeader", "TOFBarrelHits"},
      {"TOFBarrelRawHitDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "TOFBarrelRawHitLinkDigi",
#endif
       "TOFBarrelRawHitAssociationDigi"},
      {
          .threshold      = 6.0 * dd4hep::keV,
          .timeResolution = 0.025, // [ns]
      },
      app))->SetLevel(JEventLevel::Timeslice));

  // Convert raw digitized hits into hits with geometry info (ready for tracking)
  app->Add((new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "TOFBarrelRecHitDigi", {"TOFBarrelRawHitDigi"}, // Input data collection tags
      {"TOFBarrelRecHitDigi"},                        // Output data tag
      {},
      app))->SetLevel(JEventLevel::Timeslice)); // Hit reco default config for factories

  // Convert raw digitized hits into calibrated hits
  // time walk correction is still TBD
  app->Add((new JOmniFactoryGeneratorT<LGADHitCalibration_factory>(
      "TOFBarrelCalibratedHitDigi", {"TOFBarrelADCTDCDigi"}, // Input data collection tags
      {"TOFBarrelCalibratedHitDigi"},                       // Output data tag
      {},
      app))->SetLevel(JEventLevel::Timeslice)); // Hit reco default config for factories

  // cluster all hits in a sensor into one hit location
  // Currently it's just a simple weighted average
  // More sophisticated algorithm TBD
  app->Add((new JOmniFactoryGeneratorT<LGADHitClustering_factory>(
      "TOFBarrelClusterHitDigi", {"TOFBarrelCalibratedHitDigi"}, // Input data collection tags
      {"TOFBarrelClusterHitDigi"},                               // Output data tag
      {}, app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<SiliconChargeSharing_factory>(
      "TOFBarrelSharedHitDigi", {"TOFBarrelHitDigi"}, {"TOFBarrelSharedHitDigi"},
      {
          .sigma_mode     = SiliconChargeSharingConfig::ESigmaMode::rel,
          .sigma_sharingx = 1,
          .sigma_sharingy = 0.5,
          .min_edep       = 0.0 * edm4eic::unit::GeV,
          .readout        = "TOFBarrelHits",
      },
      app))->SetLevel(JEventLevel::Timeslice));

  // calculation of the extreme values for Landau distribution can be found on lin 514-520 of
  // https://root.cern.ch/root/html524/src/TMath.cxx.html#fsokrB Landau reaches minimum for mpv =
  // 0 and sigma = 1 at x = -0.22278
  const double x_when_landau_min = -0.22278;
  const double landau_min        = TMath::Landau(x_when_landau_min, 0, 1, true);
  const double sigma_analog      = 0.293951 * edm4eic::unit::ns;
  const double Vm                = 3e-4 * dd4hep::GeV;
  const double adc_range         = 256;
  // gain is set such that pulse reaches a height of adc_range when EDep = Vm
  // gain is negative as LGAD voltage is always negative
  const double gain = -adc_range / Vm / landau_min * sigma_analog;
  const int offset  = 3;
  app->Add((new JOmniFactoryGeneratorT<PulseGeneration_factory<edm4hep::SimTrackerHit>>(
      "LGADPulseGenerationDigi", {"TOFBarrelSharedHitDigi"}, {"TOFBarrelSmoothPulseDigi"},
      {
          .pulse_shape_function = "LandauPulse",
          .pulse_shape_params   = {gain, sigma_analog, offset},
          .ignore_thres         = 0.05 * adc_range,
          .timestep             = 0.01 * edm4eic::unit::ns,
      },
      app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<PulseCombiner_factory>(
      "TOFBarrelPulseCombinerDigi", {"TOFBarrelSmoothPulseDigi"}, {"TOFBarrelCombinedPulseDigi"},
      {
          .minimum_separation = 25 * edm4eic::unit::ns,
      },
      app))->SetLevel(JEventLevel::Timeslice));

  double risetime = 0.45 * edm4eic::unit::ns;
  app->Add((new JOmniFactoryGeneratorT<SiliconPulseDiscretization_factory>(
      "TOFBarrelPulseDigi", {"TOFBarrelCombinedPulseDigi"}, {"TOFBarrelPulseDigi"},
      {
          .EICROC_period = 25 * edm4eic::unit::ns,
          .local_period  = 25 * edm4eic::unit::ns / 1024,
          .global_offset = -offset * sigma_analog + risetime,
      },
      app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CFDROCDigitization_factory>(
      "CFDROCDigitizationDigi", {"TOFBarrelPulseDigi"}, {"TOFBarrelADCTDCDigi"}, {}, app))->SetLevel(JEventLevel::Timeslice));
}
// } // extern "C"
