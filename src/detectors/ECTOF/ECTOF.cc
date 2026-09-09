// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025, Dmitry Romanov,  Wouter Deconinck, Kolja Kauder, Barak Schmookler, Honey Khindri, Dmitry Kalinkin

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
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
#include "extensions/jana/EventBuilderEnabled.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/EICROCDigitization_factory.h"
#include "factories/digi/PulseCombiner_factory.h"
#include "factories/digi/PulseGeneration_factory.h"
#include "factories/digi/SiliconChargeSharing_factory.h"
#include "factories/digi/SiliconPulseDiscretization_factory.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/tracking/LGADHitClustering_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // cluster all hits in a sensor into one hit location
  // Currently it's just a simple weighted average
  // More sophisticated algorithm TBD
  app->Add(new JOmniFactoryGeneratorT<LGADHitClustering_factory>(
      "TOFEndcapClusterHits", {"TOFEndcapSharedRecHits"}, // Input data collection tags
      {"TOFEndcapClusterHits"},                           // Output data tag
      {
          .readout = "TOFEndcapHits",
          .useAve  = true,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<SiliconChargeSharing_factory>(
      "TOFEndcapSharedHits", {"TOFEndcapHits"}, {"TOFEndcapSharedHits"},
      {

          .sigma_mode     = SiliconChargeSharingConfig::ESigmaMode::rel,
          .sigma_sharingx = 0.5,
          .sigma_sharingy = 0.5,
          .min_edep       = 6 * dd4hep::keV,
          .readout        = "TOFEndcapHits",
      },
      app));

  // temporary steps to bypass pulse digitization and jump right from ChargeSharing to clusters
  // Avoid efficiency loss until we can simulate hardware accurately
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "TOFEndcapSharedRawHits", {"EventHeader", "TOFEndcapSharedHits"},
      {"TOFEndcapSharedRawHits", "TOFEndcapSharedRawHitLinks", "TOFEndcapSharedRawHitAssociations"},
      {
          .threshold      = 0.0,
          .timeResolution = 0.025, // [ns]
      },
      app));

  // Convert raw digitized hits into hits with geometry info (ready for tracking)
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "TOFEndcapSharedRecHits", {"TOFEndcapSharedRawHits"}, // Input data collection tags
      {"TOFEndcapSharedRecHits"},                           // Output data tag
      {
          // AC-LGAD: encode the REAL ~25 ps resolution (see BTOF.cc note).
          .timeResolution = 0.025, // [ns]
      },
      app)); // Hit reco default config for factories

  const double x_when_landau_min = -0.22278;
  const double landau_min        = TMath::Landau(x_when_landau_min, 0, 1, true);
  const double sigma_analog      = 0.293951 * edm4eic::unit::ns;
  const double Vm                = 3e-4 * dd4hep::GeV;
  const double adc_range         = 256;

  const double gain = -adc_range / Vm / landau_min * sigma_analog;
  const int offset  = 3;
  app->Add(new JOmniFactoryGeneratorT<PulseGeneration_factory<edm4hep::SimTrackerHit>>(
      "TOFEndcapSmoothPulses", {"TOFEndcapSharedHits"}, {"TOFEndcapSmoothPulses"},
      {
          .pulse_shape_function = "LandauPulse",
          .pulse_shape_params   = {gain, sigma_analog, offset},
          .ignore_thres         = 0.05 * adc_range,
          .timestep             = 0.01 * edm4eic::unit::ns,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<PulseCombiner_factory>(
      "TOFEndcapCombinedPulses", {"TOFEndcapSmoothPulses"}, {"TOFEndcapCombinedPulses"},
      {
          .minimum_separation = 25 * edm4eic::unit::ns,
      },
      app));

  double risetime = 0.45 * edm4eic::unit::ns;
  app->Add(new JOmniFactoryGeneratorT<SiliconPulseDiscretization_factory>(
      "TOFEndcapPulses", {"TOFEndcapCombinedPulses"}, {"TOFEndcapPulses"},
      {
          .EICROC_period = 25 * edm4eic::unit::ns,
          .local_period  = 25 * edm4eic::unit::ns / 1024,
          .global_offset = -offset * sigma_analog + risetime,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<EICROCDigitization_factory>(
      "TOFEndcapADCTDC", {"TOFEndcapPulses"}, {"TOFEndcapADCTDC"}, {}, app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the
  // frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in
  // sync with the chain above.
  if (eicrecon::eventbuilderEnabled(app)) {
    app->Add((new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
        "TOFEndcapRawHitFrame", {"EventHeader", "TOFEndcapHits"},
        {"TOFEndcapRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "TOFEndcapRawHitLinkFrame",
#endif
         "TOFEndcapRawHitAssociationFrame"},
        {
            .threshold      = 6.0 * dd4hep::keV,
            .timeResolution = 0.025,
        },
        app))->SetLevel(JEventLevel::Timeslice));

    // Convert raw digitized hits into hits with geometry info (ready for tracking)
    app->Add((new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
        "TOFEndcapRecHitFrame", {"TOFEndcapRawHitFrame"}, // Input data collection tags
        {"TOFEndcapRecHitFrame"},                         // Output data tag
        {
            .timeResolution = 0.025,
        },
        app))->SetLevel(JEventLevel::Timeslice));

    // Frame-level mirror of the event-level "bypass" ClusterHits chain above
    // (ChargeSharing -> Digi -> HitReco -> LGADHitClustering), so the
    // eventbuilder's unfolder can gate TOF Measurement2D per candidate the same
    // way it gates hits and clusters. Keep configs in sync with the chain above.
    app->Add((new JOmniFactoryGeneratorT<SiliconChargeSharing_factory>(
        "TOFEndcapSharedHitFrame", {"TOFEndcapHits"}, {"TOFEndcapSharedHitFrame"},
        {
            .sigma_mode     = SiliconChargeSharingConfig::ESigmaMode::rel,
            .sigma_sharingx = 0.5,
            .sigma_sharingy = 0.5,
            .min_edep       = 6 * dd4hep::keV,
            .readout        = "TOFEndcapHits",
        },
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
        "TOFEndcapSharedRawHitFrame", {"EventHeader", "TOFEndcapSharedHitFrame"},
        {"TOFEndcapSharedRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "TOFEndcapSharedRawHitLinkFrame",
#endif
         "TOFEndcapSharedRawHitAssociationFrame"},
        {
            .threshold      = 0.0,
            .timeResolution = 0.025, // [ns]
        },
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
        "TOFEndcapSharedRecHitFrame", {"TOFEndcapSharedRawHitFrame"},
        {"TOFEndcapSharedRecHitFrame"},
        {
            .timeResolution = 0.025, // [ns] — keep in sync with the event level
        },
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<LGADHitClustering_factory>(
        "TOFEndcapClusterHitFrame", {"TOFEndcapSharedRecHitFrame"},
        {"TOFEndcapClusterHitFrame"},
        {
            .readout = "TOFEndcapHits",
            .useAve  = true,
        },
        app))->SetLevel(JEventLevel::Timeslice));
  }

  //   app->Add(new JOmniFactoryGeneratorT<SiliconChargeSharing_factory>(
  //       "TOFEndcapSharedHitFrame", {"TOFEndcapHitFrame"}, {"TOFEndcapSharedHitFrame"},
  //       {

  //           .sigma_mode     = SiliconChargeSharingConfig::ESigmaMode::rel,
  //           .sigma_sharingx = 1,
  //           .sigma_sharingy = 1,
  //           .min_edep       = 0.001 * dd4hep::keV,
  //           .readout        = "TOFEndcapHits",
  //       },
  //       app));

  //   const double x_when_landau_min = -0.22278;
  //   const double landau_min        = TMath::Landau(x_when_landau_min, 0, 1, true);
  //   const double sigma_analog      = 0.293951 * edm4eic::unit::ns;
  //   const double Vm                = 3e-4 * dd4hep::GeV;
  //   const double adc_range         = 256;

  //   const double gain = -adc_range / Vm / landau_min * sigma_analog;
  //   const int offset  = 3;
  //   app->Add(new JOmniFactoryGeneratorT<PulseGeneration_factory<edm4hep::SimTrackerHit>>(
  //       "TOFEndcapSmoothPulseFrame", {"TOFEndcapSharedHitFrame"}, {"TOFEndcapSmoothPulseFrame"},
  //       {
  //           .pulse_shape_function = "LandauPulse",
  //           .pulse_shape_params   = {gain, sigma_analog, offset},
  //           .ignore_thres         = 0.05 * adc_range,
  //           .timestep             = 0.01 * edm4eic::unit::ns,
  //       },
  //       app));

  //   app->Add(new JOmniFactoryGeneratorT<PulseCombiner_factory>(
  //       "TOFEndcapCombinedPulseFrame", {"TOFEndcapSmoothPulseFrame"}, {"TOFEndcapCombinedPulseFrame"},
  //       {
  //           .minimum_separation = 25 * edm4eic::unit::ns,
  //       },
  //       app));

  //   double risetime = 0.45 * edm4eic::unit::ns;
  //   app->Add(new JOmniFactoryGeneratorT<SiliconPulseDiscretization_factory>(
  //       "TOFEndcapPulseFrame", {"TOFEndcapCombinedPulseFrame"}, {"TOFEndcapPulseFrame"},
  //       {
  //           .EICROC_period = 25 * edm4eic::unit::ns,
  //           .local_period  = 25 * edm4eic::unit::ns / 1024,
  //           .global_offset = -offset * sigma_analog + risetime,
  //       },
  //       app));

  //   app->Add(new JOmniFactoryGeneratorT<EICROCDigitization_factory>(
  //       "TOFEndcapADCTDCFrame", {"TOFEndcapPulseFrame"}, {"TOFEndcapADCTDCFrame"}, {}, app));
}
} // extern "C"
