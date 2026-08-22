// Copyright 2022, Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <string>
#include <vector>
#include <edm4eic/unit_system.h>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // Digitization
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "B0TrackerRawHits", {"EventHeader", "B0TrackerHits"},
      {"B0TrackerRawHits", "B0TrackerRawHitLinks", "B0TrackerRawHitAssociations"},
      {
          .threshold      = 10.0 * dd4hep::keV,
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  // Convert raw digitized hits into hits with geometry info (ready for tracking)
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "B0TrackerRecHits", {"B0TrackerRawHits"}, {"B0TrackerRecHits"},
      {
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the
  // frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in
  // sync with the chain above.
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>::TypedWiring{
          .m_tag                 = "B0TrackerRawHitFrame",
          .m_default_input_tags  = {"EventHeader", "B0TrackerHits"},
          .m_default_output_tags = {"B0TrackerRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
                                    "B0TrackerRawHitLinkFrame",
#endif
                                    "B0TrackerRawHitAssociationFrame"},
          .m_default_cfg =
              {
                  .threshold      = 10.0 * dd4hep::keV,
                  .timeResolution = 8,
              },
          .level = JEventLevel::Timeslice},
      app));

  // Convert raw digitized hits into hits with geometry info (ready for tracking)
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>::TypedWiring{
          .m_tag                 = "B0TrackerRecHitFrame",
          .m_default_input_tags  = {"B0TrackerRawHitFrame"},
          .m_default_output_tags = {"B0TrackerRecHitFrame"},
          .m_default_cfg =
              {
                  .timeResolution = 8,
              },
          .level = JEventLevel::Timeslice},
      app));
}
} // extern "C"
