// Copyright 2024, Jefferson Science Associates, LLC.
// Subject to the terms in the LICENSE file found in the top-level directory.
// kuma edit

#include <JANA/JApplication.h>
#include "extensions/jana/JOmniFactoryGeneratorT.h"

#include "factories/eventbuilder/TimeAlign_factory.h" // defines Trk/Cal TimeAlign_factory aliases
#include "factories/eventbuilder/EventBuilder_factory.h"
#include "EventUnfolder.h"


void InitPlugin_digiBTOF(JApplication* app);
void InitPlugin_digiMPGD(JApplication* app);
void InitPlugin_digiBVTX(JApplication* app);
void InitPlugin_digiBTRK(JApplication* app);
void InitPlugin_digiECTRK(JApplication* app);
void InitPlugin_digiECTOF(JApplication* app);
void InitPlugin_digiB0TRK(JApplication* app);
// void InitPlugin_digiDIRC(JApplication* app);
// void InitPlugin_digiDRICH(JApplication* app);
void InitPlugin_digiFOFFMTRK(JApplication* app);
// void InitPlugin_digiPFRICH(JApplication* app);
void InitPlugin_digiLOWQ2(JApplication* app);

void InitPlugin_digiB0ECAL(JApplication* app);
void InitPlugin_digiBEMC(JApplication* app);
void InitPlugin_digiEEMC(JApplication* app);
void InitPlugin_digiFEMC(JApplication* app);
void InitPlugin_digiECHAL(JApplication* app);
void InitPlugin_digiBHCAL(JApplication* app);
void InitPlugin_digiFHCAL(JApplication* app);
void InitPlugin_digiFOFFMTRK(JApplication* app);
// void InitPlugin_digiLUMISPECCAL(JApplication* app);
// void InitPlugin_digiZDC(JApplication* app);


extern "C" {
void InitPlugin(JApplication* app) {

  // if (!app->RegisterParameter<bool>("split_timeframes", false, "Enable timeframe splitting")) {
  //   return;
  // }

  // This is the plugin initialization function that JANA will call.
  // std::vector<std::string> m_simtrackerhit_collection_names_aligned = {
  //     "B0TrackerHits_aligned",         "BackwardMPGDEndcapHits_aligned",
  //     "DIRCBarHits_aligned",           "DRICHHits_aligned",
  //     "ForwardMPGDEndcapHits_aligned", "ForwardOffMTrackerHits_aligned",
  //     "ForwardRomanPotHits_aligned",   "LumiSpecTrackerHits_aligned",
  //     "MPGDBarrelHits_aligned",        "OuterMPGDBarrelHits_aligned",
  //     "RICHEndcapNHits_aligned",       "SiBarrelHits_aligned",
  //     "TOFBarrelHits_aligned",         "TOFEndcapHits_aligned",
  //     "TaggerTrackerHits_aligned",     "TrackerEndcapHits_aligned",
  //     "VertexBarrelHits_aligned"};

  std::vector<std::string> m_simtrackerhit_collection_names_aligned = {
      "TOFBarrelTimeAlignRecHits",
      "TOFEndcapTimeAlignRecHits",
      "MPGDBarrelTimeAlignRecHits",
      "OuterMPGDBarrelTimeAlignRecHits",
      "BackwardMPGDEndcapTimeAlignRecHits",
      "ForwardMPGDEndcapTimeAlignRecHits",
      "SiBarrelVertexTimeAlignRecHits",
      "SiBarrelTrackerTimeAlignRecHits",
      "SiEndcapTrackerTimeAlignRecHits",     
      "B0TrackerTimeAlignRecHits"
    };
    // "TaggerTrackerTimeAlignRecHits",
    // "DIRCBarTimeAlignRecHits",
    //   "DRICHTimeAlignRecHits",
    //   "ForwardOffMTrackerTimeAlignRecHits",
    //   "ForwardRomanPotTimeAlignRecHits",
    //   "LumiSpecTrackerTimeAlignRecHits",
    //   "RICHEndcapNTimeAlignRecHits"

  std::vector<std::string> m_simtrackerhit_collection_names = {
    "TOFBarrelRecHitDigi",
    "TOFEndcapRecHitDigi",
    "MPGDBarrelRecHitDigi",
    "OuterMPGDBarrelRecHitDigi",
    "BackwardMPGDEndcapRecHitDigi",
    "ForwardMPGDEndcapRecHitDigi",
    "SiBarrelVertexRecHitDigi",
    "SiBarrelTrackerRecHitDigi",
    "SiEndcapTrackerRecHitDigi",
    "B0TrackerRecHitDigi"
    };   
    // "TaggerTrackerRecHitDigi",
    // "DIRCBarRecHitDigi",
    // "DRICHRecHitDigi",
    // "ForwardOffMTrackerRecHitDigi",
    // "ForwardRomanPotRecHitDigi",
    // "LumiSpecTrackerRecHitDigi",
    // "RICHEndcapNRecHitDigi"


  std::vector<std::string> m_simcalocluster_collection_names_aligned = {
      "B0ECalTimeAlignClusters",
      "EcalBarrelTimeAlignClusters",
      "EcalEndcapNTimeAlignClusters",
      "EcalEndcapPTimeAlignClusters"
    };
    // "EcalFarForwardZDCTimeAlignClusters",
    //   "EcalLumiSpecTimeAlignClusters",
    //   "HcalBarrelTimeAlignClusters",
    //   "HcalEndcapNTimeAlignClusters",
    //   "HcalEndcapPInsertTimeAlignClusters",
    //   "HcalFarForwardZDCTimeAlignClusters",
    //   "LFHCALTimeAlignClusters"

  std::vector<std::string> m_simcalocluster_collection_names = {
    "B0ECalClusterDigi",
    "EcalBarrelClusterDigi",
    "EcalEndcapNClusterDigi",
    "EcalEndcapPClusterDigi"
    };

    // "EcalFarForwardZDCClusterDigi",
    // "EcalLumiSpecClusterDigi",
    // "HcalBarrelClusterDigi",
    // "HcalEndcapNClusterDigi",
    // "HcalEndcapPInsertClusterDigi",
    // "HcalFarForwardZDCClusterDigi",
    // "LFHCALClusterDigi",
    // "EcalBarrelImagingClusterDigi",
    // "EcalBarrelScFiClusterDigi",
    // "EcalEndcapNImagingClusterDigi",
    // "EcalEndcapPImagingClusterDigi",
    // "EcalFarForwardZDCImagingClusterDigi",
    // "EcalLumiSpecImagingClusterDigi"

  InitJANAPlugin(app);

  app->Add(new JOmniFactoryGeneratorT<TrkTimeAlign_factory>(
      JOmniFactoryGeneratorT<TrkTimeAlign_factory>::TypedWiring{
          .m_tag                 = "timeAlignment",
          .m_default_input_tags  = m_simtrackerhit_collection_names,
          .m_default_output_tags = m_simtrackerhit_collection_names_aligned,
          .level                 = JEventLevel::Timeslice,
      },
      app));

    app->Add(new JOmniFactoryGeneratorT<CalTimeAlign_factory>(
      JOmniFactoryGeneratorT<CalTimeAlign_factory>::TypedWiring{
          .m_tag                 = "CalTimeAlignment",
          .m_default_input_tags  = m_simcalocluster_collection_names,
          .m_default_output_tags = m_simcalocluster_collection_names_aligned,
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // EventBuilder: finds candidate physics events in a time-frame and assigns
  // each a (t0, dt0). Runs at the Timeslice (time-frame) level. Reads the
  // time-aligned tracker hits plus MCParticles (for phys/fake labelling).
  std::vector<std::string> m_eventbuilder_input_tags = m_simtrackerhit_collection_names_aligned;
  m_eventbuilder_input_tags.push_back("MCParticles");
  app->Add(new JOmniFactoryGeneratorT<EventBuilder_factory>(
      JOmniFactoryGeneratorT<EventBuilder_factory>::TypedWiring{
          .m_tag                 = "eventBuilder",
          .m_default_input_tags  = m_eventbuilder_input_tags,
          .m_default_output_tags = {"EventCandidates"},
          .level                 = JEventLevel::Timeslice,
      },
      app));

  // Unfolder that materialises each candidate as a PhysicsEvent (no analysis).
  app->Add(new EventUnfolder());

  // app->Add(new JOmniFactoryGeneratorT<HitChecker>(
  //     jana::components::JOmniFactoryGeneratorT<HitChecker>::TypedWiring{
  //         .tag                   = "hitChecker",
  //         .level                 = JEventLevel::Timeslice,
  //         .variadic_input_names  = m_simtrackerhit_collection_names,
  //         .variadic_output_names = m_simtrackerhit_collection_names_aligned}));
  // app->Add(new JOmniFactoryGeneratorT<HitChecker>(jana::components::JOmniFactoryGeneratorT<HitChecker>::TypedWiring
  // {.tag          = "timeslice_hit_checker",
  //  .level        = JEventLevel::PhysicsEvent,
  //  .input_names  = {"TOFBarrelRecHits"},
  //  .output_names = {"hitChecker_TS"}}));


    InitPlugin_digiBTOF(app);
    InitPlugin_digiMPGD(app);
    InitPlugin_digiBVTX(app);
    InitPlugin_digiBTRK(app);
    InitPlugin_digiECTRK(app);
    InitPlugin_digiECTOF(app);
    InitPlugin_digiB0TRK(app);
    // InitPlugin_digiDIRC(app);
    // InitPlugin_digiDRICH(app);
    // InitPlugin_digiFOFFMTRK(app);
    // InitPlugin_digiPFRICH(app);
    // InitPlugin_digiLOWQ2(app);

    InitPlugin_digiB0ECAL(app);
    InitPlugin_digiBEMC(app);
    InitPlugin_digiEEMC(app);
    InitPlugin_digiFEMC(app);
    // InitPlugin_digiECHAL(app);
    // InitPlugin_digiBHCAL(app);
    // InitPlugin_digiFHCAL(app);
    InitPlugin_digiFOFFMTRK(app);
    // InitPlugin_digiLUMISPECCAL(app);
    // InitPlugin_digiZDC(app);

}
} // "C"
