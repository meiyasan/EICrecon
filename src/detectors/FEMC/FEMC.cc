// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2021 - 2025, Chao Peng, Sylvester Joosten, Whitney Armstrong, David Lawrence, Friederike Bock, Wouter Deconinck, Kolja Kauder, Sebouh Paul, Akio Ogawa

#include <DD4hep/Detector.h>
#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <format>
#include <spdlog/logger.h>
#include <cmath>
#include <gsl/pointers>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "algorithms/calorimetry/CalorimeterHitDigiConfig.h"
#include "extensions/jana/EventBuilderEnabled.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"
#include "factories/calorimetry/TrackClusterMergeSplitter_factory.h"
#include "services/geometry/dd4hep/DD4hep_service.h"
#include "services/log/Log_service.h"

extern "C" {
void InitPlugin(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);

  auto log_service = app->GetService<Log_service>();
  auto mLog        = log_service->logger("FEMC");

  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) EcalEndcapP_capADC =
      16384; //16384, assuming 14 bits. For approximate HGCROC resolution use 65536
  decltype(CalorimeterHitDigiConfig::dyRangeADC) EcalEndcapP_dyRangeADC   = 100 * dd4hep::GeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) EcalEndcapP_pedMeanADC   = 200;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) EcalEndcapP_pedSigmaADC = 2.4576;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) EcalEndcapP_resolutionTDC =
      10 * dd4hep::picosecond;
  // Hit timing resolution sigma_t(E) = sqrt((a/sqrt(E[GeV]))^2 + b^2), fit from
  // truth-matched RecHit residuals (2026-07); b is this detector's TDC-quantization
  // floor (resolutionTDC/sqrt(12)), not necessarily the full constant term.
  decltype(CalorimeterHitRecoConfig::timeErrorScale) EcalEndcapP_timeErrorScale = 4.0584;
  decltype(CalorimeterHitRecoConfig::timeErrorOffset) EcalEndcapP_timeErrorOffset = 0.0029;
  const double EcalEndcapP_sampFrac = 0.029043; // updated with ratio to ScFi model
  decltype(CalorimeterHitDigiConfig::corrMeanScale) EcalEndcapP_corrMeanScale =
      std::format("{}", 1.0 / EcalEndcapP_sampFrac); //only used for ScFi model
  const double EcalEndcapP_nPhotonPerGeV          = 1500;
  const double EcalEndcapP_PhotonCollectionEff    = 0.285;
  const unsigned long long EcalEndcapP_totalPixel = 4 * 159565ULL;

  int EcalEndcapP_homogeneousFlag = 0;
  try {
    auto detector               = app->GetService<DD4hep_service>()->detector();
    EcalEndcapP_homogeneousFlag = detector->constant<int>("EcalEndcapP_Homogeneous_ScFi");
    if (EcalEndcapP_homogeneousFlag <= 1) {
      mLog->info("Homogeneous geometry loaded");
    } else {
      mLog->info("ScFi geometry loaded");
    }
  } catch (...) {
    // Variable not present apply legacy homogeneous geometry implementation
    EcalEndcapP_homogeneousFlag = 0;
  };

  if (EcalEndcapP_homogeneousFlag <= 1) {
    app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
        "EcalEndcapPRawHits", {"EventHeader", "EcalEndcapPHits"},
        {"EcalEndcapPRawHits", "EcalEndcapPRawHitLinks", "EcalEndcapPRawHitAssociations"},
        {
            .eRes = {0.11333 * sqrt(dd4hep::GeV), 0.03,
                     0.0 * dd4hep::GeV}, // (11.333% / sqrt(E)) \oplus 3%
            .tRes = 0.0,
            .threshold =
                0.0, // 15MeV threshold for a single tower will be applied on ADC at Reco below
            .readoutType = "sipm",
            .lightYield  = EcalEndcapP_nPhotonPerGeV / EcalEndcapP_PhotonCollectionEff,
            .photonDetectionEfficiency = EcalEndcapP_PhotonCollectionEff,
            .numEffectiveSipmPixels    = EcalEndcapP_totalPixel,
            .capADC                    = EcalEndcapP_capADC,
            .capTime                   = 100, // given in ns, 4 samples in HGCROC
            .dyRangeADC                = EcalEndcapP_dyRangeADC,
            .pedMeanADC                = EcalEndcapP_pedMeanADC,
            .pedSigmaADC               = EcalEndcapP_pedSigmaADC,
            .resolutionTDC             = EcalEndcapP_resolutionTDC,
            .corrMeanScale             = "1.0",
            .readout                   = "EcalEndcapPHits",
        },
        app // TODO: Remove me once fixed
        ));
  } else if (EcalEndcapP_homogeneousFlag == 2) {
    app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
        "EcalEndcapPRawHits", {"EventHeader", "EcalEndcapPHits"},
        {"EcalEndcapPRawHits", "EcalEndcapPRawHitLinks", "EcalEndcapPRawHitAssociations"},
        {
            .eRes = {0.0, 0.022, 0.0}, // just constant term 2.2% based on MC data comparison
            .tRes = 0.0,
            .threshold =
                0.0, // 15MeV threshold for a single tower will be applied on ADC at Reco below
            .readoutType = "sipm",
            .lightYield =
                EcalEndcapP_nPhotonPerGeV / EcalEndcapP_PhotonCollectionEff / EcalEndcapP_sampFrac,
            .photonDetectionEfficiency = EcalEndcapP_PhotonCollectionEff,
            .numEffectiveSipmPixels    = EcalEndcapP_totalPixel,
            .capADC                    = EcalEndcapP_capADC,
            .capTime                   = 100, // given in ns, 4 samples in HGCROC
            .dyRangeADC                = EcalEndcapP_dyRangeADC,
            .pedMeanADC                = EcalEndcapP_pedMeanADC,
            .pedSigmaADC               = EcalEndcapP_pedSigmaADC,
            .resolutionTDC             = EcalEndcapP_resolutionTDC,
            .corrMeanScale             = EcalEndcapP_corrMeanScale,
            .readout                   = "EcalEndcapPHits",
            .fields                    = {"fiber_x", "fiber_y"},
        },
        app // TODO: Remove me once fixed
        ));
  }

  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "EcalEndcapPRecHits", {"EcalEndcapPRawHits"}, {"EcalEndcapPRecHits"},
      {
          .capADC          = EcalEndcapP_capADC,
          .dyRangeADC      = EcalEndcapP_dyRangeADC,
          .pedMeanADC      = EcalEndcapP_pedMeanADC,
          .pedSigmaADC     = EcalEndcapP_pedSigmaADC,
          .resolutionTDC   = EcalEndcapP_resolutionTDC,
          .timeErrorScale = EcalEndcapP_timeErrorScale,
          .timeErrorOffset = EcalEndcapP_timeErrorOffset,
          // thresholdADC = thresholdFactor * pedSigmaADC + thresholdValue (see
          // CalorimeterHitReco.cc). The flat thresholdValue=3 this replaces sat
          // at only ~1.2 sigma above EcalEndcapP_pedSigmaADC=2.4576 -- far too
          // close to the noise floor. A hit whose true energy is below that
          // floor still crosses threshold on noise alone about as often as
          // not, and then gets a reconstructed energy set by the noise
          // fluctuation instead of its own (tiny) deposit: confirmed 2026-08,
          // res_check.py's EcalEndcapPRecHits ΔE/E showed a median
          // Ereco/Esim of 5.29x overall, but 29x for the low-true-energy half
          // of the sample and a sane 1.46x for the high-energy half -- two
          // different populations, not one calibration offset. EEMC's own
          // threshold (thresholdValue=4.0, pedSigmaADC=1) is a real 4-sigma
          // cut; thresholdFactor here matches that same 4-sigma margin
          // instead of a flat count that assumed EEMC's much smaller sigma.
          .thresholdFactor = 4.0,
          .thresholdValue = 0.0,
          .sampFrac = "1.00", // already taken care in DIGI code above
          .readout  = "EcalEndcapPHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "EcalEndcapPTruthProtoClusters", {"EcalEndcapPRecHits", "EcalEndcapPRawHitLinks"},
      {"EcalEndcapPTruthProtoClusters"},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "EcalEndcapPIslandProtoClusters", {"EcalEndcapPRecHits"}, {"EcalEndcapPIslandProtoClusters"},
      {.adjacencyMatrix{},
       .peakNeighbourhoodMatrix{},
       .readout{},
       .sectorDist = 5.0 * dd4hep::cm,
       .localDistXY{},
       .localDistXZ{},
       .localDistYZ{},
       .globalDistRPhi{},
       .globalDistEtaPhi{},
       .dimScaledLocalDistXY          = {1.5, 1.5},
       .splitCluster                  = false,
       .minClusterHitEdep             = 0.0 * dd4hep::MeV,
       .minClusterCenterEdep          = 60.0 * dd4hep::MeV,
       .transverseEnergyProfileMetric = "dimScaledLocalDistXY",
       .transverseEnergyProfileScale  = 1.,
       .transverseEnergyProfileScaleUnits{}},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapPTruthClustersWithoutShapes",
      {
          "EcalEndcapPTruthProtoClusters", // edm4eic::ProtoClusterCollection
          "EcalEndcapPRawHitLinks",        // edm4eic::MCRecoCalorimeterHitLink
          "EcalEndcapPRawHitAssociations"  // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapPTruthClustersWithoutShapes", "EcalEndcapPTruthClusterLinksWithoutShapes",
       "EcalEndcapPTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = true},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapPTruthClusters",
      {"EcalEndcapPTruthClustersWithoutShapes", "EcalEndcapPTruthClusterLinksWithoutShapes"},
      {"EcalEndcapPTruthClusters", "EcalEndcapPTruthClusterLinks",
       "EcalEndcapPTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapPClustersWithoutShapes",
      {
          "EcalEndcapPIslandProtoClusters", // edm4eic::ProtoClusterCollection
          "EcalEndcapPRawHitLinks",         // edm4eic::MCRecoCalorimeterHitLink
          "EcalEndcapPRawHitAssociations"   // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapPClustersWithoutShapes", "EcalEndcapPClusterLinksWithoutShapes",
       "EcalEndcapPClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapPClusters",
      {"EcalEndcapPClustersWithoutShapes", "EcalEndcapPClusterLinksWithoutShapes"},
      {"EcalEndcapPClusters", "EcalEndcapPClusterLinks", "EcalEndcapPClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  app->Add(new JOmniFactoryGeneratorT<TrackClusterMergeSplitter_factory>(
      "EcalEndcapPSplitMergeProtoClusters",
      {"EcalEndcapPTrackClusterMatches", "EcalEndcapPClusters", "CalorimeterTrackProjections"},
      {"EcalEndcapPSplitMergeProtoClusters", "EcalEndcapPTrackSplitMergeProtoClusterLinks"},
      {.minSigCut                    = -2.0,
       .avgEP                        = 1.0,
       .sigEP                        = 0.10,
       .drAdd                        = 0.30,
       .surfaceToUse                 = 1,
       .transverseEnergyProfileScale = 1.0},
      app // TODO: remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapPSplitMergeClustersWithoutShapes",
      {
          "EcalEndcapPSplitMergeProtoClusters", // edm4eic::ProtoClusterCollection
          "EcalEndcapPRawHitLinks",             // edm4eic::MCRecoCalorimeterHitLink
          "EcalEndcapPRawHitAssociations" // edm4hep::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapPSplitMergeClustersWithoutShapes",
       "EcalEndcapPSplitMergeClusterLinksWithoutShapes",
       "EcalEndcapPSplitMergeClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapPSplitMergeClusters",
      {"EcalEndcapPSplitMergeClustersWithoutShapes",
       "EcalEndcapPSplitMergeClusterLinksWithoutShapes"},
      {"EcalEndcapPSplitMergeClusters", "EcalEndcapPSplitMergeClusterLinks",
       "EcalEndcapPSplitMergeClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in sync with the chain above.
  // Note: unlike the PhysicsEvent-level chain above, this mirror always uses the
  // homogeneous (SiPM-on-tile) defaults regardless of EcalEndcapP_homogeneousFlag,
  // and it omits TrackClusterMergeSplitter (requires ACTS CalorimeterTrackProjections,
  // only available at PhysicsEvent level). Add ScFi support here if/when that geometry
  // is used at Timeslice level.
  if (eicrecon::eventbuilderEnabled(app)) {
    app->Add((new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
        "EcalEndcapPRawHitFrame", {"EventHeader", "EcalEndcapPHits"},
        {"EcalEndcapPRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "EcalEndcapPRawHitLinkFrame",
#endif
         "EcalEndcapPRawHitAssociationFrame"},
        {
            .eRes                      = {0.11333 * sqrt(dd4hep::GeV), 0.03, 0.0 * dd4hep::GeV},
            .tRes                      = 0.0,
            .threshold                 = 0.0,
            .readoutType               = "sipm",
            .lightYield                = EcalEndcapP_nPhotonPerGeV / EcalEndcapP_PhotonCollectionEff,
            .photonDetectionEfficiency = EcalEndcapP_PhotonCollectionEff,
            .numEffectiveSipmPixels    = EcalEndcapP_totalPixel,
            .capADC                    = EcalEndcapP_capADC,
            .capTime                   = 100,
            .dyRangeADC                = EcalEndcapP_dyRangeADC,
            .pedMeanADC                = EcalEndcapP_pedMeanADC,
            .pedSigmaADC               = EcalEndcapP_pedSigmaADC,
            .resolutionTDC             = EcalEndcapP_resolutionTDC,
            .corrMeanScale             = "1.0",
            .readout                   = "EcalEndcapPHits",
        },
        app))->SetLevel(JEventLevel::Timeslice));
    app->Add((new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
        "EcalEndcapPRecHitFrame", {"EcalEndcapPRawHitFrame"}, {"EcalEndcapPRecHitFrame"},
        {
            .capADC          = EcalEndcapP_capADC,
            .dyRangeADC      = EcalEndcapP_dyRangeADC,
            .pedMeanADC      = EcalEndcapP_pedMeanADC,
            .pedSigmaADC     = EcalEndcapP_pedSigmaADC,
            .resolutionTDC   = EcalEndcapP_resolutionTDC,
            .timeErrorScale = EcalEndcapP_timeErrorScale,
            .timeErrorOffset = EcalEndcapP_timeErrorOffset,
            // Same 4-sigma threshold as EcalEndcapPRecHits above -- see that
            // block's comment for why the flat thresholdValue=3 this replaces
            // let noise-floor hits through.
            .thresholdFactor = 4.0,
            .thresholdValue = 0.0,
            .sampFrac        = "1.00",
            .readout         = "EcalEndcapPHits",
        },
        app))->SetLevel(JEventLevel::Timeslice));
    app->Add((new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
        "EcalEndcapPTruthProtoClusterFrame", {"EcalEndcapPRecHitFrame", "EcalEndcapPHits"},
        {"EcalEndcapPTruthProtoClusterFrame"},
        app))->SetLevel(JEventLevel::Timeslice));
    app->Add((new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
        "EcalEndcapPIslandProtoClusterFrame", {"EcalEndcapPRecHitFrame"},
        {"EcalEndcapPIslandProtoClusterFrame"},
        {
            .adjacencyMatrix{},
            .peakNeighbourhoodMatrix{},
            .readout{},
            .sectorDist                    = 5.0 * dd4hep::cm,
            .localDistXY{},
            .localDistXZ{},
            .localDistYZ{},
            .globalDistRPhi{},
            .globalDistEtaPhi{},
            .dimScaledLocalDistXY          = {1.5, 1.5},
            .splitCluster                  = false,
            .minClusterHitEdep             = 0.0 * dd4hep::MeV,
            .minClusterCenterEdep          = 60.0 * dd4hep::MeV,
            .transverseEnergyProfileMetric = "dimScaledLocalDistXY",
            .transverseEnergyProfileScale  = 1.,
            .transverseEnergyProfileScaleUnits{},
        },
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
        "EcalEndcapPTruthClustersWithoutShapeFrame",
        {
            "EcalEndcapPTruthProtoClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
            "EcalEndcapPRawHitLinkFrame",
#endif
            "EcalEndcapPRawHitAssociationFrame"
        },
        {"EcalEndcapPTruthClustersWithoutShapeFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "EcalEndcapPTruthClusterLinksWithoutShapeFrame",
#endif
         "EcalEndcapPTruthClusterAssociationsWithoutShapeFrame"},
        {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = true},
        app))->SetLevel(JEventLevel::Timeslice));
    app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
        "EcalEndcapPTruthClusterFrame",
        {"EcalEndcapPTruthClustersWithoutShapeFrame", "EcalEndcapPTruthClusterLinksWithoutShapeFrame"},
        {"EcalEndcapPTruthClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "EcalEndcapPTruthClusterLinkFrame",
#endif
         "EcalEndcapPTruthClusterAssociationFrame"},
        {.energyWeight = "log", .logWeightBase = 6.2}, app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
        "EcalEndcapPClustersWithoutShapeFrame",
        {
            "EcalEndcapPIslandProtoClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
            "EcalEndcapPRawHitLinkFrame",
#endif
            "EcalEndcapPRawHitAssociationFrame"
        },
        {"EcalEndcapPClustersWithoutShapeFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "EcalEndcapPClusterLinksWithoutShapeFrame",
#endif
         "EcalEndcapPClusterAssociationsWithoutShapeFrame"},
        {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
        app))->SetLevel(JEventLevel::Timeslice));
    app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
        "EcalEndcapPClusterFrame",
        {"EcalEndcapPClustersWithoutShapeFrame", "EcalEndcapPClusterLinksWithoutShapeFrame"},
        {"EcalEndcapPClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "EcalEndcapPClusterLinkFrame",
#endif
         "EcalEndcapPClusterAssociationFrame"},
        {.energyWeight = "log", .logWeightBase = 3.6}, app))->SetLevel(JEventLevel::Timeslice));
  }
}
}
