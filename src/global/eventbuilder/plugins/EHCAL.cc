// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)

// Copyright (C) 2022 - 2025 Sylvester Joosten, Chao, Chao Peng, Whitney Armstrong, David Lawrence, Friederike Bock, Nathan Brei, Wouter Deconinck, Dmitry Kalinkin, Derek Anderson

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <edm4eic/EDM4eicVersion.h>
#include <JANA/Utils/JTypeInfo.h>
#include <string>
#include <variant>
#include <vector>

#include "algorithms/calorimetry/CalorimeterHitDigiConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterHitsMerger_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"

// extern "C" {
void InitPlugin_EHCAL(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);
  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) HcalEndcapN_capADC =
      32768; // assuming 15 bit ADC like FHCal
  decltype(CalorimeterHitDigiConfig::dyRangeADC) HcalEndcapN_dyRangeADC =
      200 * dd4hep::MeV; // to be verified with simulations
  decltype(CalorimeterHitDigiConfig::pedMeanADC) HcalEndcapN_pedMeanADC   = 10;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) HcalEndcapN_pedSigmaADC = 2;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) HcalEndcapN_resolutionTDC =
      10 * dd4hep::picosecond;

  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "HcalEndcapNRawHitDigi", {"EventHeader", "HcalEndcapNHits"},
      {"HcalEndcapNRawHitDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "HcalEndcapNRawHitLinkDigi",
#endif
       "HcalEndcapNRawHitAssociationDigi"},
      {
          .eRes{},
          .tRes          = 0.0 * dd4hep::ns,
          .capADC        = HcalEndcapN_capADC,
          .capTime       = 100, // given in ns, 4 samples in HGCROC
          .dyRangeADC    = HcalEndcapN_dyRangeADC,
          .pedMeanADC    = HcalEndcapN_pedMeanADC,
          .pedSigmaADC   = HcalEndcapN_pedSigmaADC,
          .resolutionTDC = HcalEndcapN_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "HcalEndcapNHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "HcalEndcapNRecHitDigi", {"HcalEndcapNRawHitDigi"}, {"HcalEndcapNRecHitDigi"},
      {
          .capADC          = HcalEndcapN_capADC,
          .dyRangeADC      = HcalEndcapN_dyRangeADC,
          .pedMeanADC      = HcalEndcapN_pedMeanADC,
          .pedSigmaADC     = HcalEndcapN_pedSigmaADC,
          .resolutionTDC   = HcalEndcapN_resolutionTDC,
          .thresholdFactor = 0.0,
          .thresholdValue =
              41.0, // 0.1875 MeV deposition out of 200 MeV max (per layer) --> adc = 10 + 0.1875 / 200 * 32768 == 41
          .sampFrac =
              "0.0095", // from latest study - implement at level of reco hits rather than clusters
          .readout = "HcalEndcapNHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitsMerger_factory>(
      "HcalEndcapNMergedHitDigi", {"HcalEndcapNRecHitDigi"}, {"HcalEndcapNMergedHitDigi"},
      {.readout = "HcalEndcapNHits", .fieldTransformations = {"layer:4", "slice:0"}},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "HcalEndcapNTruthProtoClusterDigi", {"HcalEndcapNMergedHitDigi", "HcalEndcapNHits"},
      {"HcalEndcapNTruthProtoClusterDigi"},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "HcalEndcapNIslandProtoClusterDigi", {"HcalEndcapNMergedHitDigi"},
      {"HcalEndcapNIslandProtoClusterDigi"},
      {
          .adjacencyMatrix{},
          .peakNeighbourhoodMatrix{},
          .readout{},
          .sectorDist  = 5.0 * dd4hep::cm,
          .localDistXY = {15 * dd4hep::cm, 15 * dd4hep::cm},
          .localDistXZ{},
          .localDistYZ{},
          .globalDistRPhi{},
          .globalDistEtaPhi{},
          .dimScaledLocalDistXY{},
          .splitCluster                  = true,
          .minClusterHitEdep             = 0.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "globalDistEtaPhi",
          .transverseEnergyProfileScale  = 1.,
          .transverseEnergyProfileScaleUnits{},
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "HcalEndcapNTruthClustersWithoutShapeDigi",
      {
          "HcalEndcapNTruthProtoClusterDigi", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "HcalEndcapNRawHitLinkDigi", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "HcalEndcapNRawHitAssociationDigi" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"HcalEndcapNTruthClustersWithoutShapeDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "HcalEndcapNTruthClusterLinksWithoutShapeDigi",
#endif
       "HcalEndcapNTruthClusterAssociationsWithoutShapeDigi"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "HcalEndcapNTruthClusterDigi",
      {"HcalEndcapNTruthClustersWithoutShapeDigi", "HcalEndcapNTruthClusterAssociationsWithoutShapeDigi"},
      {"HcalEndcapNTruthClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "HcalEndcapNTruthClusterLinkDigi",
#endif
       "HcalEndcapNTruthClusterAssociationDigi"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "HcalEndcapNClustersWithoutShapeDigi",
      {
          "HcalEndcapNIslandProtoClusterDigi", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "HcalEndcapNRawHitLinkDigi", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "HcalEndcapNRawHitAssociationDigi" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"HcalEndcapNClustersWithoutShapeDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "HcalEndcapNClusterLinksWithoutShapeDigi",
#endif
       "HcalEndcapNClusterAssociationsWithoutShapeDigi"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 6.2,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "HcalEndcapNClusterDigi",
      {"HcalEndcapNClustersWithoutShapeDigi", "HcalEndcapNClusterAssociationsWithoutShapeDigi"},
      {"HcalEndcapNClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "HcalEndcapNClusterLinkDigi",
#endif
       "HcalEndcapNClusterAssociationDigi"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));
  // TrackClusterMergeSplitter (and its SplitMerge cluster chain) is not
  // registered at Timeslice level: it requires ACTS CalorimeterTrackProjections
  // (plus track-cluster matches), which only exist at PhysicsEvent level.
}
// }
