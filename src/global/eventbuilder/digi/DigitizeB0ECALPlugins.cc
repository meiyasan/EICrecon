// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Whitney Armstrong, Sylvester Joosten, Chao Peng, David Lawrence, Wouter Deconinck, Kolja Kauder, Nathan Brei, Dmitry Kalinkin, Derek Anderson, Michael Pitt

#include <edm4eic/EDM4eicVersion.h>
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

// extern "C" {
void InitPlugin_digiB0ECAL(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);

  app->Add((new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "B0ECalRawHitDigi", {"EventHeader", "B0ECalHits"},
      {"B0ECalRawHitDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalRawHitLinkDigi",
#endif
       "B0ECalRawHitAssociationDigi"},
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
      "B0ECalRecHitDigi", {"B0ECalRawHitDigi"}, {"B0ECalRecHitDigi"},
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
      "B0ECalTruthProtoClusterDigi", {"B0ECalRecHitDigi", "B0ECalHits"}, {"B0ECalTruthProtoClusterDigi"},
      app))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "B0ECalIslandProtoClusterDigi", {"B0ECalRecHitDigi"}, {"B0ECalIslandProtoClusterDigi"},
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
      "B0ECalClustersWithoutShapeDigi",
      {
          "B0ECalIslandProtoClusterDigi", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "B0ECalRawHitLinkDigi", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "B0ECalRawHitAssociationDigi" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalClustersWithoutShapeDigi", // edm4eic::Cluster
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalClusterLinksWithoutShapeDigi",
#endif
       "B0ECalClusterAssociationsWithoutShapeDigi"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
      app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalClusterDigi", {"B0ECalClustersWithoutShapeDigi", "B0ECalClusterAssociationsWithoutShapeDigi"},
      {"B0ECalClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalClusterLinkDigi",
#endif
       "B0ECalClusterAssociationDigi"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "B0ECalTruthClustersWithoutShapeDigi",
      {
          "B0ECalTruthProtoClusterDigi", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "B0ECalRawHitLinkDigi", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "B0ECalRawHitAssociationDigi" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalTruthClustersWithoutShapeDigi", // edm4eic::Cluster
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalTruthClusterLinksWithoutShapeDigi",
#endif
       "B0ECalTruthClusterAssociationsWithoutShapeDigi"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalTruthClusterDigi",
      {"B0ECalTruthClustersWithoutShapeDigi", "B0ECalTruthClusterAssociationsWithoutShapeDigi"},
      {"B0ECalTruthClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "B0ECalTruthClusterLinkDigi",
#endif
       "B0ECalTruthClusterAssociationDigi"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app))->SetLevel(JEventLevel::Timeslice));
}
// }
