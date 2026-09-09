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
      "SiBarrelRawHits", {"EventHeader", "SiBarrelHits"},
      {"SiBarrelRawHits", "SiBarrelRawHitLinks", "SiBarrelRawHitAssociations"},
      {
          .threshold = 0.54 * dd4hep::keV,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<RandomNoisePixel_factory>(
      "SiBarrelNoiseRawHits", {"EventHeader"}, {"SiBarrelNoiseRawHits"},
      {.addNoise = true, .noise_rate_per_pixel_per_event = 2.0e-7, .readout_name = "SiBarrelHits"},
      app));

  app->Add(new JOmniFactoryGeneratorT<CollectionCollector_factory<edm4eic::RawTrackerHit>>(
      "SiBarrelRawHitsWithNoise",                  // Name of the combiner instance
      {"SiBarrelRawHits", "SiBarrelNoiseRawHits"}, // Inputs: original + noise-only
      {"SiBarrelRawHitsWithNoise"},                // Output: merged collection
      {},                                          // default config
      app));

  // Convert raw digitized hits into hits with geometry info (ready for tracking)
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "SiBarrelTrackerRecHits", {"SiBarrelRawHitsWithNoise"}, {"SiBarrelTrackerRecHits"},
      {}, // default config
      app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the
  // frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in
  // sync with the chain above.
  if (eicrecon::eventbuilderEnabled(app)) {
    app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
        JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>::TypedWiring{
            .m_tag                 = "SiBarrelRawHitFrame",
            .m_default_input_tags  = {"EventHeader", "SiBarrelHits"},
            .m_default_output_tags = {"SiBarrelRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
                                      "SiBarrelRawHitLinkFrame",
#endif
                                      "SiBarrelRawHitAssociationFrame"},
            .m_default_cfg =
                {
                    .threshold = 0.54 * dd4hep::keV,
                },
            .level = JEventLevel::Timeslice},
        app));

    // Convert raw digitized hits into hits with geometry info (ready for tracking)
    app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
        JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>::TypedWiring{
            .m_tag                 = "SiBarrelTrackerRecHitFrame",
            .m_default_input_tags  = {"SiBarrelRawHitFrame"},
            .m_default_output_tags = {"SiBarrelTrackerRecHitFrame"},
            .m_default_cfg =
                {
                    .timeResolution = 10,
                },
            .level = JEventLevel::Timeslice},
        app));
  }
}
} // extern "C"
