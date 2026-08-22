// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Whitney Armstrong, Sylvester Joosten, Chao Peng, David Lawrence, Wouter Deconinck, Kolja Kauder, Nathan Brei, Dmitry Kalinkin, Derek Anderson, Michael Pitt

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <cmath>
#include <string>
#include <variant>
#include <vector>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);

  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "B0ECalRawHits", {"EventHeader", "B0ECalHits"},
      {"B0ECalRawHits", "B0ECalRawHitLinks", "B0ECalRawHitAssociations"},
      {
          // The stochastic term is set using light yield in PbOW4 of N_photons = 145.75 / GeV / mm, for 6x6 mm2 sensors with PDE=0.18 (a=1/sqrt(145.75*36*0.18))
          .eRes          = {0.0326 * sqrt(dd4hep::GeV), 0.00, 0.0 * dd4hep::GeV},
          .tRes          = 0.0 * dd4hep::ns,
          .threshold     = 5.0 * dd4hep::MeV,
          .capADC        = 16384,
          .dyRangeADC    = 170 * dd4hep::GeV,
          .pedMeanADC    = 100,
          .pedSigmaADC   = 1,
          .resolutionTDC = 1e-11,
          .corrMeanScale = "1.0",
          .readout       = "B0ECalHits",
      },
      app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "B0ECalRecHits", {"B0ECalRawHits"}, {"B0ECalRecHits"},
      {
          .capADC          = 16384,
          .dyRangeADC      = 170. * dd4hep::GeV,
          .pedMeanADC      = 100,
          .pedSigmaADC     = 1,
          .resolutionTDC   = 1e-11,
          .thresholdFactor = 0.0,
          .thresholdValue  = 1.0, // using threshold of 10 photons = 10 MeV = 1 ADC
          .sampFrac        = "0.998",
          .readout         = "B0ECalHits",
          .sectorField     = "sector",
      },
      app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "B0ECalTruthProtoClusters", {"B0ECalRecHits", "B0ECalRawHitLinks"},
      {"B0ECalTruthProtoClusters"}, app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "B0ECalIslandProtoClusters", {"B0ECalRecHits"}, {"B0ECalIslandProtoClusters"},
      {
          .adjacencyMatrix{},
          .peakNeighbourhoodMatrix{},
          .readout{},
          .sectorDist = 5.0 * dd4hep::cm,
          .localDistXY{},
          .localDistXZ{},
          .localDistYZ{},
          .globalDistRPhi{},
          .globalDistEtaPhi{},
          .dimScaledLocalDistXY          = {1.8, 1.8},
          .splitCluster                  = false,
          .minClusterHitEdep             = 1.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "globalDistEtaPhi",
          .transverseEnergyProfileScale  = 1.,
          .transverseEnergyProfileScaleUnits{},
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "B0ECalClustersWithoutShapes",
      {
          "B0ECalIslandProtoClusters", // edm4eic::ProtoClusterCollection
          "B0ECalRawHitLinks",         // edm4eic::MCRecoCalorimeterHitLink
          "B0ECalRawHitAssociations"   // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalClustersWithoutShapes", // edm4eic::Cluster
       "B0ECalClusterLinksWithoutShapes",
       "B0ECalClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
      app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalClusters", {"B0ECalClustersWithoutShapes", "B0ECalClusterLinksWithoutShapes"},
      {"B0ECalClusters", "B0ECalClusterLinks", "B0ECalClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "B0ECalTruthClustersWithoutShapes",
      {
          "B0ECalTruthProtoClusters", // edm4eic::ProtoClusterCollection
          "B0ECalRawHitLinks",        // edm4eic::MCRecoCalorimeterHitLink
          "B0ECalRawHitAssociations"  // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalTruthClustersWithoutShapes", // edm4eic::Cluster
       "B0ECalTruthClusterLinksWithoutShapes",
       "B0ECalTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalTruthClusters",
      {"B0ECalTruthClustersWithoutShapes", "B0ECalTruthClusterLinksWithoutShapes"},
      {"B0ECalTruthClusters", "B0ECalTruthClusterLinks", "B0ECalTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in sync with the chain above.
  app->Add((new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "B0ECalRawHitFrame", {"EventHeader", "B0ECalHits"},
      {"B0ECalRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalRawHitLinkFrame",
#endif
       "B0ECalRawHitAssociationFrame"},
      {
          // The stochastic term is set using light yield in PbOW4 of N_photons = 145.75 / GeV / mm, for 6x6 mm2 sensors with PDE=0.18 (a=1/sqrt(145.75*36*0.18))
          .eRes          = {0.0326 * sqrt(dd4hep::GeV), 0.00, 0.0 * dd4hep::GeV},
          .tRes          = 0.0 * dd4hep::ns,
          .threshold     = 5.0 * dd4hep::MeV,
          .capADC        = 16384,
          .dyRangeADC    = 170 * dd4hep::GeV,
          .pedMeanADC    = 100,
          .pedSigmaADC   = 1,
          .resolutionTDC = 1e-11,
          .corrMeanScale = "1.0",
          .readout       = "B0ECalHits",
      },
      app))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "B0ECalRecHitFrame", {"B0ECalRawHitFrame"}, {"B0ECalRecHitFrame"},
      {
          .capADC          = 16384,
          .dyRangeADC      = 170. * dd4hep::GeV,
          .pedMeanADC      = 100,
          .pedSigmaADC     = 1,
          .resolutionTDC   = 1e-11,
          .thresholdFactor = 0.0,
          .thresholdValue  = 1.0, // using threshold of 10 photons = 10 MeV = 1 ADC
          .sampFrac        = "0.998",
          .readout         = "B0ECalHits",
          .sectorField     = "sector",
      },
      app))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "B0ECalTruthProtoClusterFrame", {"B0ECalRecHitFrame", "B0ECalHits"}, {"B0ECalTruthProtoClusterFrame"},
      app))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "B0ECalIslandProtoClusterFrame", {"B0ECalRecHitFrame"}, {"B0ECalIslandProtoClusterFrame"},
      {
          .adjacencyMatrix{},
          .peakNeighbourhoodMatrix{},
          .readout{},
          .sectorDist = 5.0 * dd4hep::cm,
          .localDistXY{},
          .localDistXZ{},
          .localDistYZ{},
          .globalDistRPhi{},
          .globalDistEtaPhi{},
          .dimScaledLocalDistXY          = {1.8, 1.8},
          .splitCluster                  = false,
          .minClusterHitEdep             = 1.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "globalDistEtaPhi",
          .transverseEnergyProfileScale  = 1.,
          .transverseEnergyProfileScaleUnits{},
      },
      app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "B0ECalClustersWithoutShapeFrame",
      {
          "B0ECalIslandProtoClusterFrame", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "B0ECalRawHitLinkFrame", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "B0ECalRawHitAssociationFrame" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalClustersWithoutShapeFrame", // edm4eic::Cluster
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalClusterLinksWithoutShapeFrame",
#endif
       "B0ECalClusterAssociationsWithoutShapeFrame"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
      app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalClusterFrame", {"B0ECalClustersWithoutShapeFrame", "B0ECalClusterAssociationsWithoutShapeFrame"},
      {"B0ECalClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalClusterLinkFrame",
#endif
       "B0ECalClusterAssociationFrame"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "B0ECalTruthClustersWithoutShapeFrame",
      {
          "B0ECalTruthProtoClusterFrame", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "B0ECalRawHitLinkFrame", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "B0ECalRawHitAssociationFrame" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalTruthClustersWithoutShapeFrame", // edm4eic::Cluster
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalTruthClusterLinksWithoutShapeFrame",
#endif
       "B0ECalTruthClusterAssociationsWithoutShapeFrame"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalTruthClusterFrame",
      {"B0ECalTruthClustersWithoutShapeFrame", "B0ECalTruthClusterAssociationsWithoutShapeFrame"},
      {"B0ECalTruthClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalTruthClusterLinkFrame",
#endif
       "B0ECalTruthClusterAssociationFrame"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app))->SetLevel(JEventLevel::Timeslice));
}
}
