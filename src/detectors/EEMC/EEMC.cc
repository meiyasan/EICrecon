// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Sylvester Joosten, Chao, Chao Peng, Whitney Armstrong, Thomas Britton, David Lawrence, Dhevan Gangadharan, Wouter Deconinck, Dmitry Kalinkin, Derek Anderson

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
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
#include "factories/calorimetry/CalorimeterParticleIDPostML_factory.h"
#include "factories/calorimetry/CalorimeterParticleIDPreML_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"
#include "factories/calorimetry/TrackClusterMergeSplitter_factory.h"
#include "factories/meta/ONNXInference_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);

  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) EcalEndcapN_capADC         = 16384; //65536,  16bit ADC
  decltype(CalorimeterHitDigiConfig::dyRangeADC) EcalEndcapN_dyRangeADC = 20.0 * dd4hep::GeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) EcalEndcapN_pedMeanADC = 20;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) EcalEndcapN_pedSigmaADC = 1;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) EcalEndcapN_resolutionTDC =
      10 * dd4hep::picosecond;
  // Hit timing resolution sigma_t(E) = sqrt((a/sqrt(E[GeV]))^2 + b^2), fit from
  // truth-matched RecHit residuals (2026-07); b is this detector's TDC-quantization
  // floor (resolutionTDC/sqrt(12)), not necessarily the full constant term.
  decltype(CalorimeterHitRecoConfig::timeErrorScale) EcalEndcapN_timeErrorScale = 0.2257;
  decltype(CalorimeterHitRecoConfig::timeErrorOffset) EcalEndcapN_timeErrorOffset = 0.0029;
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "EcalEndcapNRawHits", {"EventHeader", "EcalEndcapNHits"},
      {"EcalEndcapNRawHits", "EcalEndcapNRawHitLinks", "EcalEndcapNRawHitAssociations"},
      {
          .eRes        = {0.0 * sqrt(dd4hep::GeV), 0.0, 0.0 * dd4hep::GeV},
          .tRes        = 0.0 * dd4hep::ns,
          .threshold   = 0.0 * dd4hep::MeV, // Use ADC cut instead
          .readoutType = "sipm",
          // 18. pe/MeV is measured with PMT at 25% QE
          .lightYield = 18. / 0.25 / dd4hep::MeV,
          // Based on slide 6 of https://indico.bnl.gov/event/29076/contributions/110749/attachments/63706/109457/Calo_meeting_Jun25_Updated.pdf
          // Geometric factor for 16 of 3x3 mm^2 sensors covering 20x20 mm^2 area for sensor with 28% QE
          .photonDetectionEfficiency = (16 * (3. * 3.) / (20. * 20.)) * 0.28,
          // S14160-3015PS, 16 sensors per cell
          .numEffectiveSipmPixels = 39984ULL * 16,
          .capADC                 = EcalEndcapN_capADC,
          .dyRangeADC             = EcalEndcapN_dyRangeADC,
          .pedMeanADC             = EcalEndcapN_pedMeanADC,
          .pedSigmaADC            = EcalEndcapN_pedSigmaADC,
          .resolutionTDC          = EcalEndcapN_resolutionTDC,
          .corrMeanScale          = "1.0",
          .readout                = "EcalEndcapNHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "EcalEndcapNRecHits", {"EcalEndcapNRawHits"}, {"EcalEndcapNRecHits"},
      {
          .capADC          = EcalEndcapN_capADC,
          .dyRangeADC      = EcalEndcapN_dyRangeADC,
          .pedMeanADC      = EcalEndcapN_pedMeanADC,
          .pedSigmaADC     = EcalEndcapN_pedSigmaADC,
          .resolutionTDC   = EcalEndcapN_resolutionTDC,
          .timeErrorScale = EcalEndcapN_timeErrorScale,
          .timeErrorOffset = EcalEndcapN_timeErrorOffset,
          .thresholdFactor = 0.0,
          .thresholdValue  = 4.0, // (20. GeV / 16384) * 4 ~= 5 MeV
          .sampFrac        = "0.96",
          .readout         = "EcalEndcapNHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "EcalEndcapNTruthProtoClusters", {"EcalEndcapNRecHits", "EcalEndcapNRawHitLinks"},
      {"EcalEndcapNTruthProtoClusters"},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "EcalEndcapNIslandProtoClusters", {"EcalEndcapNRecHits"}, {"EcalEndcapNIslandProtoClusters"},
      {
          .adjacencyMatrix         = "(abs(row_1 - row_2) + abs(column_1 - column_2)) == 1",
          .peakNeighbourhoodMatrix = "max(abs(row_1 - row_2), abs(column_1 - column_2)) == 1",
          .readout                 = "EcalEndcapNHits",
          .sectorDist              = 5.0 * dd4hep::cm,
          .localDistXY{},
          .localDistXZ{},
          .localDistYZ{},
          .globalDistRPhi{},
          .globalDistEtaPhi{},
          .dimScaledLocalDistXY{},
          .splitCluster                  = true,
          .minClusterHitEdep             = 1.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "globalDistEtaPhi",
          .transverseEnergyProfileScale  = 0.08,
          .transverseEnergyProfileScaleUnits{},
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapNTruthClustersWithoutShapes",
      {
          "EcalEndcapNTruthProtoClusters", // edm4eic::ProtoClusterCollection
          "EcalEndcapNRawHitLinks",        // edm4eic::MCRecoCalorimeterHitLink
          "EcalEndcapNRawHitAssociations"  // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapNTruthClustersWithoutShapes", "EcalEndcapNTruthClusterLinksWithoutShapes",
       "EcalEndcapNTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 4.6, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapNTruthClusters",
      {"EcalEndcapNTruthClustersWithoutShapes", "EcalEndcapNTruthClusterLinksWithoutShapes"},
      {"EcalEndcapNTruthClusters", "EcalEndcapNTruthClusterLinks",
       "EcalEndcapNTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 4.6}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapNClustersWithoutPIDAndShapes",
      {
          "EcalEndcapNIslandProtoClusters", // edm4eic::ProtoClusterCollection
          "EcalEndcapNRawHitLinks",         // edm4eic::MCRecoCalorimeterHitLink
          "EcalEndcapNRawHitAssociations"   // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapNClustersWithoutPIDAndShapes",             // edm4eic::Cluster
       "EcalEndcapNClusterLinksWithoutPIDAndShapes",         // edm4eic::MCRecoClusterParticleLink
       "EcalEndcapNClusterAssociationsWithoutPIDAndShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapNClustersWithoutPID",
      {"EcalEndcapNClustersWithoutPIDAndShapes", "EcalEndcapNClusterLinksWithoutPIDAndShapes"},
      {"EcalEndcapNClustersWithoutPID", "EcalEndcapNClusterLinksWithoutPID",
       "EcalEndcapNClusterAssociationsWithoutPID"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterParticleIDPreML_factory>(
      "EcalEndcapNParticleIDPreML",
      {
          "EcalEndcapNClustersWithoutPID",
          "EcalEndcapNClusterLinksWithoutPID",
      },
      {
          "EcalEndcapNParticleIDInput_features",
          "EcalEndcapNParticleIDTarget",
      },
      app));
  app->Add(new JOmniFactoryGeneratorT<ONNXInference_factory>(
      "EcalEndcapNParticleIDInference",
      {
          "EcalEndcapNParticleIDInput_features",
      },
      {
          "EcalEndcapNParticleIDOutput_label",
          "EcalEndcapNParticleIDOutput_probability_tensor",
      },
      {
          .modelPath = "calibrations/onnx/EcalEndcapN_pi_rejection.onnx",
      },
      app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterParticleIDPostML_factory>(
      "EcalEndcapNParticleIDPostML",
      {
          "EcalEndcapNClustersWithoutPID",
          "EcalEndcapNClusterLinksWithoutPID",
          "EcalEndcapNParticleIDOutput_probability_tensor",
      },

      {
          "EcalEndcapNClusters",
          "EcalEndcapNClusterLinks",
          "EcalEndcapNClusterAssociations",
          "EcalEndcapNClusterParticleIDs",
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<TrackClusterMergeSplitter_factory>(
      "EcalEndcapNSplitMergeProtoClusters",
      {"EcalEndcapNTrackClusterMatches", "EcalEndcapNClustersWithoutPID",
       "CalorimeterTrackProjections"},
      {"EcalEndcapNSplitMergeProtoClusters", "EcalEndcapNTrackSplitMergeProtoClusterLinks"},
      {.minSigCut                    = -1.0,
       .avgEP                        = 1.0,
       .sigEP                        = 0.10,
       .drAdd                        = 0.08,
       .surfaceToUse                 = 1,
       .transverseEnergyProfileScale = 1.0},
      app // TODO: remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapNSplitMergeClustersWithoutShapes",
      {"EcalEndcapNSplitMergeProtoClusters",
       "EcalEndcapNRawHitLinks", // edm4eic::MCRecoCalorimeterHitLink
       "EcalEndcapNRawHitAssociations"},
      {"EcalEndcapNSplitMergeClustersWithoutShapes",
       "EcalEndcapNSplitMergeClusterLinksWithoutShapes",
       "EcalEndcapNSplitMergeClusterAssociationsWithoutShapes"},
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapNSplitMergeClusters",

      {"EcalEndcapNSplitMergeClustersWithoutShapes",
       "EcalEndcapNSplitMergeClusterLinksWithoutShapes"},
      {"EcalEndcapNSplitMergeClusters", "EcalEndcapNSplitMergeClusterLinks",
       "EcalEndcapNSplitMergeClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  // Timeslice-level mirror of the chain above for eventbuilder ("Frame" suffix marks the frame-level variant, avoiding collision with the PhysicsEvent-level names above). Keep in sync with the chain above.
  // TrackClusterMergeSplitter (and its SplitMerge cluster chain) is not
  // registered at Timeslice level: it requires ACTS CalorimeterTrackProjections
  // (plus track-cluster matches), which only exist at PhysicsEvent level.
  app->Add((new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "EcalEndcapNRawHitFrame", {"EventHeader", "EcalEndcapNHits"},
      {"EcalEndcapNRawHitFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNRawHitLinkFrame",
#endif
       "EcalEndcapNRawHitAssociationFrame"},
      {
          .eRes        = {0.0 * sqrt(dd4hep::GeV), 0.0, 0.0 * dd4hep::GeV},
          .tRes        = 0.0 * dd4hep::ns,
          .threshold   = 0.0 * dd4hep::MeV, // Use ADC cut instead
          .readoutType = "sipm",
          // 18. pe/MeV is measured with PMT at 25% QE
          .lightYield = 18. / 0.25 / dd4hep::MeV,
          // Based on slide 6 of https://indico.bnl.gov/event/29076/contributions/110749/attachments/63706/109457/Calo_meeting_Jun25_Updated.pdf
          // Geometric factor for 16 of 3x3 mm^2 sensors covering 20x20 mm^2 area for sensor with 28% QE
          .photonDetectionEfficiency = (16 * (3. * 3.) / (20. * 20.)) * 0.28,
          // S14160-3015PS, 16 sensors per cell
          .numEffectiveSipmPixels = 39984ULL * 16,
          .capADC                 = EcalEndcapN_capADC,
          .dyRangeADC             = EcalEndcapN_dyRangeADC,
          .pedMeanADC             = EcalEndcapN_pedMeanADC,
          .pedSigmaADC            = EcalEndcapN_pedSigmaADC,
          .resolutionTDC          = EcalEndcapN_resolutionTDC,
          .corrMeanScale          = "1.0",
          .readout                = "EcalEndcapNHits",
      },
      app // TODO: Remove me once fixed
      ))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "EcalEndcapNRecHitFrame", {"EcalEndcapNRawHitFrame"}, {"EcalEndcapNRecHitFrame"},
      {
          .capADC          = EcalEndcapN_capADC,
          .dyRangeADC      = EcalEndcapN_dyRangeADC,
          .pedMeanADC      = EcalEndcapN_pedMeanADC,
          .pedSigmaADC     = EcalEndcapN_pedSigmaADC,
          .resolutionTDC   = EcalEndcapN_resolutionTDC,
          .timeErrorScale = EcalEndcapN_timeErrorScale,
          .timeErrorOffset = EcalEndcapN_timeErrorOffset,
          .thresholdFactor = 0.0,
          .thresholdValue  = 4.0, // (20. GeV / 16384) * 4 ~= 5 MeV
          .sampFrac        = "0.96",
          .readout         = "EcalEndcapNHits",
      },
      app // TODO: Remove me once fixed
      ))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "EcalEndcapNTruthProtoClusterFrame", {"EcalEndcapNRecHitFrame", "EcalEndcapNHits"},
      {"EcalEndcapNTruthProtoClusterFrame"},
      app // TODO: Remove me once fixed
      ))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "EcalEndcapNIslandProtoClusterFrame", {"EcalEndcapNRecHitFrame"}, {"EcalEndcapNIslandProtoClusterFrame"},
      {
          .adjacencyMatrix         = "(abs(row_1 - row_2) + abs(column_1 - column_2)) == 1",
          .peakNeighbourhoodMatrix = "max(abs(row_1 - row_2), abs(column_1 - column_2)) == 1",
          .readout                 = "EcalEndcapNHits",
          .sectorDist              = 5.0 * dd4hep::cm,
          .localDistXY{},
          .localDistXZ{},
          .localDistYZ{},
          .globalDistRPhi{},
          .globalDistEtaPhi{},
          .dimScaledLocalDistXY{},
          .splitCluster                  = true,
          .minClusterHitEdep             = 1.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "globalDistEtaPhi",
          .transverseEnergyProfileScale  = 0.08,
          .transverseEnergyProfileScaleUnits{},
      },
      app // TODO: Remove me once fixed
      ))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapNTruthClustersWithoutShapeFrame",
      {
          "EcalEndcapNTruthProtoClusterFrame", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapNRawHitLinkFrame", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "EcalEndcapNRawHitAssociationFrame" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapNTruthClustersWithoutShapeFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNTruthClusterLinksWithoutShapeFrame",
#endif
       "EcalEndcapNTruthClusterAssociationsWithoutShapeFrame"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 4.6, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapNTruthClusterFrame",
      {"EcalEndcapNTruthClustersWithoutShapeFrame", "EcalEndcapNTruthClusterLinksWithoutShapeFrame"},
      {"EcalEndcapNTruthClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNTruthClusterLinkFrame",
#endif
       "EcalEndcapNTruthClusterAssociationFrame"},
      {.energyWeight = "log", .logWeightBase = 4.6}, app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapNClustersWithoutPIDAndShapeFrame",
      {
          "EcalEndcapNIslandProtoClusterFrame", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapNRawHitLinkFrame", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "EcalEndcapNRawHitAssociationFrame" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapNClustersWithoutPIDAndShapeFrame", // edm4eic::Cluster
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNClusterLinksWithoutPIDAndShapeFrame", // edm4eic::MCRecoClusterParticleLink
#endif
       "EcalEndcapNClusterAssociationsWithoutPIDAndShapeFrame"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapNClustersWithoutPIDFrame",
      {"EcalEndcapNClustersWithoutPIDAndShapeFrame",
       "EcalEndcapNClusterLinksWithoutPIDAndShapeFrame"},
      {"EcalEndcapNClustersWithoutPIDFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNClusterLinksWithoutPIDFrame",
#endif
       "EcalEndcapNClusterAssociationsWithoutPIDFrame"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app))->SetLevel(JEventLevel::Timeslice));

  app->Add((new JOmniFactoryGeneratorT<CalorimeterParticleIDPreML_factory>(
      "EcalEndcapNParticleIDPreMLFrame",
      {
          "EcalEndcapNClustersWithoutPIDFrame",
          "EcalEndcapNClusterLinksWithoutPIDFrame",
      },
      {
          "EcalEndcapNParticleIDInput_Tffeatures",
          "EcalEndcapNParticleIDTargetFrame",
      },
      app))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<ONNXInference_factory>(
      "EcalEndcapNParticleIDInferenceFrame",
      {
          "EcalEndcapNParticleIDInput_Tffeatures",
      },
      {
          "EcalEndcapNParticleIDOutput_Tflabel",
          "EcalEndcapNParticleIDOutput_probability_Tftensor",
      },
      {
          .modelPath = "calibrations/onnx/EcalEndcapN_pi_rejection.onnx",
      },
      app))->SetLevel(JEventLevel::Timeslice));
  app->Add((new JOmniFactoryGeneratorT<CalorimeterParticleIDPostML_factory>(
      "EcalEndcapNParticleIDPostMLFrame",
      {
          "EcalEndcapNClustersWithoutPIDFrame",
          "EcalEndcapNClusterLinksWithoutPIDFrame",
          "EcalEndcapNParticleIDOutput_probability_Tftensor",
      },

      {
          "EcalEndcapNClusterFrame",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapNClusterLinkFrame",
#endif
          "EcalEndcapNClusterAssociationFrame",
          "EcalEndcapNClusterParticleIDFrame",
      },
      app))->SetLevel(JEventLevel::Timeslice));
}
}
