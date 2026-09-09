// Copyright 2023, Alex Jentsch
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <string>
#include <vector>
#include <edm4eic/unit_system.h>

#include "algorithms/fardetectors/MatrixTransferStaticConfig.h"
#include "extensions/jana/EventBuilderEnabled.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/fardetectors/MatrixTransferStatic_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  using namespace eicrecon;

  //Digitized hits, especially for thresholds
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "ForwardOffMTrackerRawHits", {"EventHeader", "ForwardOffMTrackerHits"},
      {"ForwardOffMTrackerRawHits", "ForwardOffMTrackerRawHitLinks",
       "ForwardOffMTrackerRawHitAssociations"},
      {
          .threshold      = 10.0 * dd4hep::keV,
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "ForwardOffMTrackerRecHits", {"ForwardOffMTrackerRawHits"}, {"ForwardOffMTrackerRecHits"},
      {
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<MatrixTransferStatic_factory>(
      "ForwardOffMRecParticles", {"MCParticles", "ForwardOffMTrackerRecHits"},
      {"ForwardOffMRecParticles"},
      {
          .matrix_configs = {{
              .nomMomentum = 130.0,

              .aX =
                  {
                      {2.08344, 5.37571},
                      {0.188756, -2.90941},
                  },

              .aY =
                  {
                      {-0.977013, -35.7785},
                      {-0.0812252, -2.86315},
                  },

              .local_x_offset       = -1032.2,
              .local_y_offset       = 0.00462829,
              .local_x_slope_offset = -59.7363,
              .local_y_slope_offset = -0.0030213,

          }},

          .hit1minZ = 25490.0,
          .hit1maxZ = 25512.0,
          .hit2minZ = 27012.0,
          .hit2maxZ = 27035.0,

          .readout = "ForwardOffMTrackerRecHits",
      },
      app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the
  // frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in
  // sync with the chain above.
  if (eicrecon::eventbuilderEnabled(app)) {
    app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
        JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>::TypedWiring{
            .m_tag                 = "ForwardOffMTrackerRawHitFrame",
            .m_default_input_tags  = {"EventHeader", "ForwardOffMTrackerHits"},
            .m_default_output_tags = {"ForwardOffMTrackerRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
                                      "ForwardOffMTrackerRawHitLinkFrame",
#endif
                                      "ForwardOffMTrackerRawHitAssociationFrame"},
            .m_default_cfg =
                {
                    .threshold      = 10.0 * dd4hep::keV,
                    .timeResolution = 8,
                },
            .level = JEventLevel::Timeslice},
        app));

    app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
        JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>::TypedWiring{
            .m_tag                 = "ForwardOffMTrackerRecHitFrame",
            .m_default_input_tags  = {"ForwardOffMTrackerRawHitFrame"},
            .m_default_output_tags = {"ForwardOffMTrackerRecHitFrame"},
            .m_default_cfg =
                {
                    .timeResolution = 8,
                },
            .level = JEventLevel::Timeslice},
        app));

    app->Add(new JOmniFactoryGeneratorT<MatrixTransferStatic_factory>(
        JOmniFactoryGeneratorT<MatrixTransferStatic_factory>::TypedWiring{
            .m_tag                 = "ForwardOffMRecParticleFrame",
            .m_default_input_tags  = {"MCParticles", "ForwardOffMTrackerRecHitFrame"},
            .m_default_output_tags = {"ForwardOffMRecParticleFrame"},
            .m_default_cfg =
                {
                    .matrix_configs = {{
                        .nomMomentum = 130.0,

                        .aX =
                            {
                                {1.61591, 12.6786},
                                {0.184206, -2.907},
                            },

                        .aY =
                            {
                                {-0.789385, -28.5578},
                                {-0.0721796, -2.8763},
                            },

                        .local_x_offset       = -881.631,
                        .local_y_offset       = -0.00552173,
                        .local_x_slope_offset = -59.7386,
                        .local_y_slope_offset = -0.00360656,

                    }},

                    .hit1minZ = 22490.0,
                    .hit1maxZ = 22512.0,
                    .hit2minZ = 24512.0,
                    .hit2maxZ = 24535.0,

                    .readout = "ForwardOffMTrackerRecHitFrame",
                },
            .level = JEventLevel::Timeslice},
        app));
  }
}
}
