// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)

// Copyright (C) 2022 - 2025 Chao Peng, Sylvester Joosten, Whitney Armstrong, David Lawrence, Friederike Bock, Wouter Deconinck, Kolja Kauder, Sebouh Paul, Akio Ogawa

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/JApplication.h>
#include <edm4eic/EDM4eicVersion.h>
#include <cmath>
#include <string>
#include <variant>
#include <vector>

#include "algorithms/calorimetry/CalorimeterHitDigiConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"

// Note: the standard detectors/FEMC/FEMC.cc detects EcalEndcapP_Homogeneous_ScFi via
// DD4hep_service to switch between homogeneous-SiPM-on-tile and ScFi parameters.
// This plugin always uses the homogeneous (SiPM-on-tile) defaults, which are correct
// for the standard ePIC geometry. Add ScFi support if/when that geometry is used.
//
// TrackClusterMergeSplitter is also omitted: it requires CalorimeterTrackProjections
// produced by ACTS at PhysicsEvent level, which is not available at Timeslice level.

// extern "C" {
void InitPlugin_FEMC(JApplication* app) {

  using namespace eicrecon;

  // Register all factories below at the Timeslice (time-frame) level.
#define ADD_TS(theApp, ...) theApp->Add((__VA_ARGS__)->SetLevel(JEventLevel::Timeslice))

  InitJANAPlugin(app);

  // Homogeneous SiPM-on-tile parameters (standard ePIC geometry).
  decltype(CalorimeterHitDigiConfig::capADC) EcalEndcapP_capADC             = 16384;
  decltype(CalorimeterHitDigiConfig::dyRangeADC) EcalEndcapP_dyRangeADC     = 100 * dd4hep::GeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) EcalEndcapP_pedMeanADC     = 200;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) EcalEndcapP_pedSigmaADC   = 2.4576;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) EcalEndcapP_resolutionTDC =
      10 * dd4hep::picosecond;
  const double EcalEndcapP_nPhotonPerGeV         = 1500;
  const double EcalEndcapP_PhotonCollectionEff   = 0.285;
  const unsigned long long EcalEndcapP_totalPixel = 4 * 159565ULL;

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "EcalEndcapPRawHitDigi", {"EventHeader", "EcalEndcapPHits"},
      {"EcalEndcapPRawHitDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapPRawHitLinkDigi",
#endif
       "EcalEndcapPRawHitAssociationDigi"},
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
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "EcalEndcapPRecHitDigi", {"EcalEndcapPRawHitDigi"}, {"EcalEndcapPRecHitDigi"},
      {
          .capADC          = EcalEndcapP_capADC,
          .dyRangeADC      = EcalEndcapP_dyRangeADC,
          .pedMeanADC      = EcalEndcapP_pedMeanADC,
          .pedSigmaADC     = EcalEndcapP_pedSigmaADC,
          .resolutionTDC   = EcalEndcapP_resolutionTDC,
          .thresholdFactor = 0.0,
          .thresholdValue  = 3, // ≈ 15.25 MeV
          .sampFrac        = "1.00",
          .readout         = "EcalEndcapPHits",
      },
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "EcalEndcapPTruthProtoClusterDigi", {"EcalEndcapPRecHitDigi", "EcalEndcapPHits"},
      {"EcalEndcapPTruthProtoClusterDigi"},
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "EcalEndcapPIslandProtoClusterDigi", {"EcalEndcapPRecHitDigi"},
      {"EcalEndcapPIslandProtoClusterDigi"},
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
      app));

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapPTruthClustersWithoutShapeDigi",
      {
          "EcalEndcapPTruthProtoClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapPRawHitLinkDigi",
#endif
          "EcalEndcapPRawHitAssociationDigi"
      },
      {"EcalEndcapPTruthClustersWithoutShapeDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapPTruthClusterLinksWithoutShapeDigi",
#endif
       "EcalEndcapPTruthClusterAssociationsWithoutShapeDigi"},
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = true},
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapPTruthClusterDigi",
      {"EcalEndcapPTruthClustersWithoutShapeDigi", "EcalEndcapPTruthClusterAssociationsWithoutShapeDigi"},
      {"EcalEndcapPTruthClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapPTruthClusterLinkDigi",
#endif
       "EcalEndcapPTruthClusterAssociationDigi"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapPClustersWithoutShapeDigi",
      {
          "EcalEndcapPIslandProtoClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapPRawHitLinkDigi",
#endif
          "EcalEndcapPRawHitAssociationDigi"
      },
      {"EcalEndcapPClustersWithoutShapeDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapPClusterLinksWithoutShapeDigi",
#endif
       "EcalEndcapPClusterAssociationsWithoutShapeDigi"},
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapPClusterDigi",
      {"EcalEndcapPClustersWithoutShapeDigi", "EcalEndcapPClusterAssociationsWithoutShapeDigi"},
      {"EcalEndcapPClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapPClusterLinkDigi",
#endif
       "EcalEndcapPClusterAssociationDigi"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));
}
// }
