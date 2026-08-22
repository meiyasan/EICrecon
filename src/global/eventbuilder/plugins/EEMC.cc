// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)

// Copyright (C) 2022 - 2025 Sylvester Joosten, Chao, Chao Peng, Whitney Armstrong, Thomas Britton, David Lawrence, Dhevan Gangadharan, Wouter Deconinck, Dmitry Kalinkin, Derek Anderson

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <edm4eic/EDM4eicVersion.h>
#include <JANA/Utils/JTypeInfo.h>
#include <cmath>
#include <string>
#include <variant>
#include <vector>

#include "algorithms/calorimetry/CalorimeterHitDigiConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterParticleIDPostML_factory.h"
#include "factories/calorimetry/CalorimeterParticleIDPreML_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"
#include "factories/meta/ONNXInference_factory.h"

// extern "C" {
void InitPlugin_EEMC(JApplication* app) {

  using namespace eicrecon;

  // Register all factories below at the Timeslice (time-frame) level.
#define ADD_TS(theApp, ...) theApp->Add((__VA_ARGS__)->SetLevel(JEventLevel::Timeslice))

  InitJANAPlugin(app);

  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) EcalEndcapN_capADC         = 16384; //65536,  16bit ADC
  decltype(CalorimeterHitDigiConfig::dyRangeADC) EcalEndcapN_dyRangeADC = 20.0 * dd4hep::GeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) EcalEndcapN_pedMeanADC = 20;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) EcalEndcapN_pedSigmaADC = 1;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) EcalEndcapN_resolutionTDC =
      10 * dd4hep::picosecond;
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "EcalEndcapNRawHitDigi", {"EventHeader", "EcalEndcapNHits"},
      {"EcalEndcapNRawHitDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNRawHitLinkDigi",
#endif
       "EcalEndcapNRawHitAssociationDigi"},
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
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "EcalEndcapNRecHitDigi", {"EcalEndcapNRawHitDigi"}, {"EcalEndcapNRecHitDigi"},
      {
          .capADC          = EcalEndcapN_capADC,
          .dyRangeADC      = EcalEndcapN_dyRangeADC,
          .pedMeanADC      = EcalEndcapN_pedMeanADC,
          .pedSigmaADC     = EcalEndcapN_pedSigmaADC,
          .resolutionTDC   = EcalEndcapN_resolutionTDC,
          .thresholdFactor = 0.0,
          .thresholdValue  = 4.0, // (20. GeV / 16384) * 4 ~= 5 MeV
          .sampFrac        = "0.96",
          .readout         = "EcalEndcapNHits",
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "EcalEndcapNTruthProtoClusterDigi", {"EcalEndcapNRecHitDigi", "EcalEndcapNHits"},
      {"EcalEndcapNTruthProtoClusterDigi"},
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "EcalEndcapNIslandProtoClusterDigi", {"EcalEndcapNRecHitDigi"}, {"EcalEndcapNIslandProtoClusterDigi"},
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

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapNTruthClustersWithoutShapeDigi",
      {
          "EcalEndcapNTruthProtoClusterDigi", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapNRawHitLinkDigi", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "EcalEndcapNRawHitAssociationDigi" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapNTruthClustersWithoutShapeDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNTruthClusterLinksWithoutShapeDigi",
#endif
       "EcalEndcapNTruthClusterAssociationsWithoutShapeDigi"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 4.6, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapNTruthClusterDigi",
      {"EcalEndcapNTruthClustersWithoutShapeDigi", "EcalEndcapNTruthClusterAssociationsWithoutShapeDigi"},
      {"EcalEndcapNTruthClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNTruthClusterLinkDigi",
#endif
       "EcalEndcapNTruthClusterAssociationDigi"},
      {.energyWeight = "log", .logWeightBase = 4.6}, app));

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalEndcapNClustersWithoutPIDAndShapeDigi",
      {
          "EcalEndcapNIslandProtoClusterDigi", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapNRawHitLinkDigi", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "EcalEndcapNRawHitAssociationDigi" // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"EcalEndcapNClustersWithoutPIDAndShapeDigi", // edm4eic::Cluster
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNClusterLinksWithoutPIDAndShapeDigi", // edm4eic::MCRecoClusterParticleLink
#endif
       "EcalEndcapNClusterAssociationsWithoutPIDAndShapeDigi"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalEndcapNClustersWithoutPIDDigi",
      {"EcalEndcapNClustersWithoutPIDAndShapeDigi",
       "EcalEndcapNClusterAssociationsWithoutPIDAndShapeDigi"},
      {"EcalEndcapNClustersWithoutPIDDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalEndcapNClusterLinksWithoutPIDDigi",
#endif
       "EcalEndcapNClusterAssociationsWithoutPIDDigi"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  // TrackClusterMergeSplitter (and its SplitMerge cluster chain) is not
  // registered at Timeslice level: it requires ACTS CalorimeterTrackProjections
  // (plus track-cluster matches), which only exist at PhysicsEvent level.

  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterParticleIDPreML_factory>(
      "EcalEndcapNParticleIDPreMLDigi",
      {
          "EcalEndcapNClustersWithoutPIDDigi",
          "EcalEndcapNClusterAssociationsWithoutPIDDigi",
      },
      {
          "EcalEndcapNParticleIDInput_Tffeatures",
          "EcalEndcapNParticleIDTargetDigi",
      },
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<ONNXInference_factory>(
      "EcalEndcapNParticleIDInferenceDigi",
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
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterParticleIDPostML_factory>(
      "EcalEndcapNParticleIDPostMLDigi",
      {
          "EcalEndcapNClustersWithoutPIDDigi",
          "EcalEndcapNClusterAssociationsWithoutPIDDigi",
          "EcalEndcapNParticleIDOutput_probability_Tftensor",
      },

      {
          "EcalEndcapNClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalEndcapNClusterLinkDigi",
#endif
          "EcalEndcapNClusterAssociationDigi",
          "EcalEndcapNClusterParticleIDDigi",
      },
      app));
}
// }
