// Copyright 2022, Dmitry Romanov, Minjung Kim, Joshua Sobaljic, Shujie Li
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <edm4eic/RawTrackerHit.h>
#include <memory>
#include <string>
#include <vector>

#include "extensions/jana/EventBuilderEnabled.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/RandomNoisePixel_factory.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/meta/CollectionCollector_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // Digitization
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "SiEndcapTrackerRawHits", {"EventHeader", "TrackerEndcapHits"},
      {"SiEndcapTrackerRawHits", "SiEndcapTrackerRawHitLinks", "SiEndcapTrackerRawHitAssociations"},
      {
          .threshold = 0.54 * dd4hep::keV,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<RandomNoisePixel_factory>(
      "SiEndcapTrackerNoiseRawHits", {"EventHeader"}, {"SiEndcapTrackerNoiseRawHits"},
      {.addNoise                       = true,
       .noise_rate_per_pixel_per_event = 2.0e-7,
       .readout_name                   = "TrackerEndcapHits"},
      app));
  app->Add(new JOmniFactoryGeneratorT<CollectionCollector_factory<edm4eic::RawTrackerHit>>(
      "SiEndcapTrackerRawHitsWithNoise", {"SiEndcapTrackerRawHits", "SiEndcapTrackerNoiseRawHits"},
      {"SiEndcapTrackerRawHitsWithNoise"}, {}, app));
  // Convert raw digitized hits into hits with geometry info (ready for tracking)
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "SiEndcapTrackerRecHits", {"SiEndcapTrackerRawHitsWithNoise"}, {"SiEndcapTrackerRecHits"},
      {}, // default config
      app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the
  // frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in
  // sync with the chain above.
  if (eicrecon::eventbuilderEnabled(app)) {
    app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
        JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>::TypedWiring{
            .m_tag                 = "SiEndcapTrackerRawHitFrame",
            .m_default_input_tags  = {"EventHeader", "TrackerEndcapHits"},
            .m_default_output_tags = {"SiEndcapTrackerRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
                                      "SiEndcapTrackerRawHitLinkFrame",
#endif
                                      "SiEndcapTrackerRawHitAssociationFrame"},
            .m_default_cfg =
                {
                    .threshold = 0.54 * dd4hep::keV,
                },
            .level = JEventLevel::Timeslice},
        app));

    // Convert raw digitized hits into hits with geometry info (ready for tracking)
    app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
        JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>::TypedWiring{
            .m_tag                 = "SiEndcapTrackerRecHitFrame",
            .m_default_input_tags  = {"SiEndcapTrackerRawHitFrame"},
            .m_default_output_tags = {"SiEndcapTrackerRecHitFrame"},
            .m_default_cfg =
                {
                    .timeResolution = 10,
                },
            .level = JEventLevel::Timeslice},
        app));
  }
}
} // extern "C"
