// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024 David Lawrence, Derek Anderson, Wouter Deconinck

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <string>
#include <variant>
#include <vector>

#include "algorithms/calorimetry/CalorimeterHitDigiConfig.h"
#include "algorithms/calorimetry/CalorimeterIslandClusterConfig.h"
#include "extensions/jana/EventBuilderEnabled.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterHitsMerger_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"
#include "factories/calorimetry/TrackClusterMergeSplitter_factory.h"

extern "C" {

void InitPlugin(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);

  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) HcalBarrel_capADC         = 65536; //65536,  16bit ADC
  decltype(CalorimeterHitDigiConfig::dyRangeADC) HcalBarrel_dyRangeADC = 1.0 * dd4hep::GeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) HcalBarrel_pedMeanADC = 300;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) HcalBarrel_pedSigmaADC = 2;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) HcalBarrel_resolutionTDC =
      1 * dd4hep::picosecond;
  // Hit timing resolution sigma_t(E) = sqrt((a/sqrt(E[GeV]))^2 + b^2), fit from
  // truth-matched RecHit residuals (2026-07); b is this detector's TDC-quantization
  // floor (resolutionTDC/sqrt(12)), not necessarily the full constant term.
  decltype(CalorimeterHitRecoConfig::timeErrorScale) HcalBarrel_timeErrorScale  = 1.4793;
  decltype(CalorimeterHitRecoConfig::timeErrorOffset) HcalBarrel_timeErrorOffset = 0.0003;

  // Set default adjacency matrix. Magic constants:
  //  320 - number of tiles per row
  decltype(CalorimeterIslandClusterConfig::adjacencyMatrix) HcalBarrel_adjacencyMatrix =
      "("
      // check for vertically adjacent tiles
      "  ( (abs(eta_1 - eta_2) == 1) && (abs(phi_1 - phi_2) == 0) ) ||"
      // check for horizontally adjacent tiles
      "  ( (abs(eta_1 - eta_2) == 0) && (abs(phi_1 - phi_2) == 1) ) ||"
      // check for horizontally adjacent tiles at wraparound
      "  ( (abs(eta_1 - eta_2) == 0) && (abs(phi_1 - phi_2) == (320 - 1)) )"
      ") == 1";

  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "HcalBarrelRawHits", {"EventHeader", "HcalBarrelHits"},
      {"HcalBarrelRawHits", "HcalBarrelRawHitLinks", "HcalBarrelRawHitAssociations"},
      {
          .eRes          = {},
          .tRes          = 0.0 * dd4hep::ns,
          .threshold     = 0.0, // Use ADC cut instead
          .capADC        = HcalBarrel_capADC,
          .capTime       = 100, // given in ns, 4 samples in HGCROC
          .dyRangeADC    = HcalBarrel_dyRangeADC,
          .pedMeanADC    = HcalBarrel_pedMeanADC,
          .pedSigmaADC   = HcalBarrel_pedSigmaADC,
          .resolutionTDC = HcalBarrel_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "HcalBarrelHits",
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "HcalBarrelRecHits", {"HcalBarrelRawHits"}, {"HcalBarrelRecHits"},
      {
          .capADC          = HcalBarrel_capADC,
          .dyRangeADC      = HcalBarrel_dyRangeADC,
          .pedMeanADC      = HcalBarrel_pedMeanADC,
          .pedSigmaADC     = HcalBarrel_pedSigmaADC, // not used; relying on energy cut
          .resolutionTDC   = HcalBarrel_resolutionTDC,
          .timeErrorScale  = HcalBarrel_timeErrorScale,
          .timeErrorOffset = HcalBarrel_timeErrorOffset,
          .thresholdFactor = 0.0,     // not used; relying on flat ADC cut
          .thresholdValue  = 33,      // pedSigmaADC + thresholdValue = half-MIP (333 ADC)
          .sampFrac        = "0.033", // average, from sPHENIX simulations
          .readout         = "HcalBarrelHits",
          .layerField      = "",
          .sectorField     = "",
      },
      app // TODO: Remove me once fixed
      ));

  // --------------------------------------------------------------------
  // If needed, merge adjacent phi tiles into towers. By default,
  // NO merging will be done. This can be changed at runtime.
  // --------------------------------------------------------------------
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitsMerger_factory>(
      "HcalBarrelMergedHits", {"HcalBarrelRecHits"}, {"HcalBarrelMergedHits"},
      {.readout = "HcalBarrelHits", .fieldTransformations = {"phi:phi"}},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "HcalBarrelTruthProtoClusters", {"HcalBarrelRecHits", "HcalBarrelRawHitLinks"},
      {"HcalBarrelTruthProtoClusters"},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "HcalBarrelIslandProtoClusters", {"HcalBarrelRecHits"}, {"HcalBarrelIslandProtoClusters"},
      {.adjacencyMatrix = HcalBarrel_adjacencyMatrix,
       .peakNeighbourhoodMatrix{},
       .readout    = "HcalBarrelHits",
       .sectorDist = 5.0 * dd4hep::cm,
       .localDistXY{},
       .localDistXZ{},
       .localDistYZ{},
       .globalDistRPhi{},
       .globalDistEtaPhi{},
       .dimScaledLocalDistXY{},
       .splitCluster                  = false,
       .minClusterHitEdep             = 5.0 * dd4hep::MeV,
       .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
       .transverseEnergyProfileMetric = "globalDistEtaPhi",
       .transverseEnergyProfileScale  = 1.,
       .transverseEnergyProfileScaleUnits{}},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "HcalBarrelClustersWithoutShapes",
      {
          "HcalBarrelIslandProtoClusters", // edm4eic::ProtoClusterCollection
          "HcalBarrelRawHitLinks",         // edm4eic::MCRecoCalorimeterHitLink
          "HcalBarrelRawHitAssociations"   // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"HcalBarrelClustersWithoutShapes", "HcalBarrelClusterLinksWithoutShapes",
       "HcalBarrelClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "HcalBarrelClusters",
      {"HcalBarrelClustersWithoutShapes", "HcalBarrelClusterLinksWithoutShapes"},
      {"HcalBarrelClusters", "HcalBarrelClusterLinks", "HcalBarrelClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "HcalBarrelTruthClustersWithoutShapes",
      {
          "HcalBarrelTruthProtoClusters", // edm4eic::ProtoClusterCollection
          "HcalBarrelRawHitLinks",        // edm4eic::MCRecoCalorimeterHitLink
          "HcalBarrelRawHitAssociations"  // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"HcalBarrelTruthClustersWithoutShapes", "HcalBarrelTruthClusterLinksWithoutShapes",
       "HcalBarrelTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "HcalBarrelTruthClusters",
      {"HcalBarrelTruthClustersWithoutShapes", "HcalBarrelTruthClusterLinksWithoutShapes"},
      {"HcalBarrelTruthClusters", "HcalBarrelTruthClusterLinks",
       "HcalBarrelTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  app->Add(new JOmniFactoryGeneratorT<TrackClusterMergeSplitter_factory>(
      "HcalBarrelSplitMergeProtoClusters",
      {"HcalBarrelTrackClusterMatches", "HcalBarrelClustersWithoutShapes",
       "CalorimeterTrackProjections"},
      {"HcalBarrelSplitMergeProtoClusters", "HcalBarrelTrackSplitMergeProtoClusterLinks"},
      {.minSigCut                    = -2.0,
       .avgEP                        = 0.50,
       .sigEP                        = 0.25,
       .drAdd                        = 0.40,
       .surfaceToUse                 = 1,
       .transverseEnergyProfileScale = 1.0},
      app // TODO: remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "HcalBarrelSplitMergeClustersWithoutShapes",
      {
          "HcalBarrelSplitMergeProtoClusters", // edm4eic::ProtoClusterCollection
          "HcalBarrelRawHitLinks",             // edm4eic::MCRecoCalorimeterHitLink
          "HcalBarrelRawHitAssociations"       // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"HcalBarrelSplitMergeClustersWithoutShapes", "HcalBarrelSplitMergeClusterLinksWithoutShapes",
       "HcalBarrelSplitMergeClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "HcalBarrelSplitMergeClusters",
      {"HcalBarrelSplitMergeClustersWithoutShapes",
       "HcalBarrelSplitMergeClusterLinksWithoutShapes"},
      {"HcalBarrelSplitMergeClusters", "HcalBarrelSplitMergeClusterLinks",
       "HcalBarrelSplitMergeClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the
  // frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep
  // in sync with the chain above.
  if (eicrecon::eventbuilderEnabled(app)) {
    app->Add((new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
        "HcalBarrelRawHitFrame", {"EventHeader", "HcalBarrelHits"},
        {"HcalBarrelRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "HcalBarrelRawHitLinkFrame",
#endif
         "HcalBarrelRawHitAssociationFrame"},
        {
            .eRes          = {},
            .tRes          = 0.0 * dd4hep::ns,
            .threshold     = 0.0,
            .capADC        = HcalBarrel_capADC,
            .capTime       = 100,
            .dyRangeADC    = HcalBarrel_dyRangeADC,
            .pedMeanADC    = HcalBarrel_pedMeanADC,
            .pedSigmaADC   = HcalBarrel_pedSigmaADC,
            .resolutionTDC = HcalBarrel_resolutionTDC,
            .corrMeanScale = "1.0",
            .readout       = "HcalBarrelHits",
        },
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
        "HcalBarrelRecHitFrame", {"HcalBarrelRawHitFrame"}, {"HcalBarrelRecHitFrame"},
        {
            .capADC          = HcalBarrel_capADC,
            .dyRangeADC      = HcalBarrel_dyRangeADC,
            .pedMeanADC      = HcalBarrel_pedMeanADC,
            .pedSigmaADC     = HcalBarrel_pedSigmaADC,
            .resolutionTDC   = HcalBarrel_resolutionTDC,
            .timeErrorScale  = HcalBarrel_timeErrorScale,
            .timeErrorOffset = HcalBarrel_timeErrorOffset,
            .thresholdFactor = 0.0,
            .thresholdValue  = 33,
            .sampFrac        = "0.033",
            .readout         = "HcalBarrelHits",
            .layerField      = "",
            .sectorField     = "",
        },
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterHitsMerger_factory>(
        "HcalBarrelMergedHitFrame", {"HcalBarrelRecHitFrame"}, {"HcalBarrelMergedHitFrame"},
        {.readout = "HcalBarrelHits", .fieldTransformations = {"phi:phi"}},
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
        "HcalBarrelTruthProtoClusterFrame", {"HcalBarrelRecHitFrame", "HcalBarrelHits"},
        {"HcalBarrelTruthProtoClusterFrame"},
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
        "HcalBarrelIslandProtoClusterFrame", {"HcalBarrelRecHitFrame"},
        {"HcalBarrelIslandProtoClusterFrame"},
        {.adjacencyMatrix = HcalBarrel_adjacencyMatrix,
         .peakNeighbourhoodMatrix{},
         .readout    = "HcalBarrelHits",
         .sectorDist = 5.0 * dd4hep::cm,
         .localDistXY{},
         .localDistXZ{},
         .localDistYZ{},
         .globalDistRPhi{},
         .globalDistEtaPhi{},
         .dimScaledLocalDistXY{},
         .splitCluster                  = false,
         .minClusterHitEdep             = 5.0 * dd4hep::MeV,
         .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
         .transverseEnergyProfileMetric = "globalDistEtaPhi",
         .transverseEnergyProfileScale  = 1.,
         .transverseEnergyProfileScaleUnits{}},
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
        "HcalBarrelClustersWithoutShapeFrame",
        {
            "HcalBarrelIslandProtoClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
            "HcalBarrelRawHitLinkFrame",
#endif
            "HcalBarrelRawHitAssociationFrame"
        },
        {"HcalBarrelClustersWithoutShapeFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "HcalBarrelClusterLinksWithoutShapeFrame",
#endif
         "HcalBarrelClusterAssociationsWithoutShapeFrame"},
        {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
        app))->SetLevel(JEventLevel::Timeslice));

    app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
        "HcalBarrelClusterFrame",
        {"HcalBarrelClustersWithoutShapeFrame", "HcalBarrelClusterLinksWithoutShapeFrame"},
        {"HcalBarrelClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
         "HcalBarrelClusterLinkFrame",
#endif
         "HcalBarrelClusterAssociationFrame"},
        {.energyWeight = "log", .logWeightBase = 6.2}, app))->SetLevel(JEventLevel::Timeslice));
  }
}
}
