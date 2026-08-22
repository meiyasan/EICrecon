// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026, Marco Meyer-Conde (ARL, Tokyo City University)
//                     Takuya Kumaoka (QNSI, The University of Tokyo)

// Copyright (C) 2022 - 2025 Whitney Armstrong, Sylvester Joosten, Chao Peng, David Lawrence, Thomas Britton, Wouter Deconinck, Maria Zurek, Akshaya Vijay, Nathan Brei, Dmitry Kalinkin, Derek Anderson, Minho Kim

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <edm4eic/EDM4eicVersion.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "algorithms/calorimetry/CalorimeterHitDigiConfig.h"
#include "algorithms/calorimetry/ImagingTopoClusterConfig.h"
#include "algorithms/calorimetry/SimCalorimeterHitProcessorConfig.h"
#include "algorithms/digi/CALOROCDigitizationConfig.h"
#include "algorithms/digi/PulseCombinerConfig.h"
#include "algorithms/digi/PulseGenerationConfig.h"
#include "algorithms/digi/PulseNoiseConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/EnergyPositionClusterMerger_factory.h"
#include "factories/calorimetry/ImagingClusterReco_factory.h"
#include "factories/calorimetry/ImagingTopoCluster_factory.h"
#include "factories/calorimetry/SimCalorimeterHitProcessor_factory.h"
#include "factories/calorimetry/TruthEnergyPositionClusterMerger_factory.h"
#include "factories/digi/PulseCombiner_factory.h"
#include "factories/digi/PulseGeneration_factory.h"
#include "factories/digi/PulseNoise_factory.h"

#if EDM4EIC_VERSION_MAJOR > 8 || (EDM4EIC_VERSION_MAJOR == 8 && EDM4EIC_VERSION_MINOR >= 7)
#include "factories/digi/CALOROCDigitization_factory.h"
#endif

// extern "C" {
void InitPlugin_BEMC(JApplication* app) {

  using namespace eicrecon;

  // Register all factories below at the Timeslice (time-frame) level.
#define ADD_TS(theApp, ...) theApp->Add((__VA_ARGS__)->SetLevel(JEventLevel::Timeslice))

  InitJANAPlugin(app);

  // Make sure left and right use the same value
  decltype(SimCalorimeterHitProcessorConfig::attenuationParameters) EcalBarrelScFi_attPars = {
      0.416212, 747.39875 * edm4eic::unit::mm, 7521.88383 * edm4eic::unit::mm};
  decltype(SimCalorimeterHitProcessorConfig::hitMergeFields) EcalBarrelScFi_hitMergeFields = {
      "fiber", "z"};
  decltype(SimCalorimeterHitProcessorConfig::contributionMergeFields)
      EcalBarrelScFi_contributionMergeFields = {"fiber"};
  decltype(SimCalorimeterHitProcessorConfig::inversePropagationSpeed)
      EcalBarrelScFi_inversePropagationSpeed = {(1. / 160) * edm4eic::unit::ns / edm4eic::unit::mm};
  decltype(SimCalorimeterHitProcessorConfig::fixedTimeDelay) EcalBarrelScFi_fixedTimeDelay = {
      2 * edm4eic::unit::ns};
  decltype(SimCalorimeterHitProcessorConfig::timeWindow) EcalBarrelScFi_timeWindow = {
      100 * edm4eic::unit::ns};

  decltype(PulseGenerationConfig::pulse_shape_function) EcalBarrelScFi_pulse_shape_function = {
      "LandauPulse"};
  decltype(PulseGenerationConfig::pulse_shape_params) EcalBarrelScFi_pulse_shape_params = {
      5.0, 10 * edm4eic::unit::ns};
  decltype(PulseGenerationConfig::ignore_thres) EcalBarrelScFi_ignore_thres = {5.0e-5};
  decltype(PulseGenerationConfig::timestep) EcalBarrelScFi_timestep = {0.5 * edm4eic::unit::ns};

  decltype(PulseCombinerConfig::combine_field) EcalBarrelScFi_combine_field           = {"grid"};
  decltype(PulseCombinerConfig::minimum_separation) EcalBarrelScFi_minimum_separation = {
      100 * edm4eic::unit::ns};
  decltype(PulseNoiseConfig::poles) EcalBarrelScFi_poles                  = {2};
  decltype(PulseNoiseConfig::variance) EcalBarrelScFi_variance            = {0.5};
  decltype(PulseNoiseConfig::alpha) EcalBarrelScFi_alpha                  = {0};
  decltype(PulseNoiseConfig::scale) EcalBarrelScFi_scale                  = {5.4e-5};
  decltype(PulseNoiseConfig::pedestal) EcalBarrelScFi_pedestal            = {1.6e-4};
  decltype(CALOROCDigitizationConfig::adc_phase) EcalBarrelScFi_adc_phase = {10 *
                                                                             edm4eic::unit::ns};
  decltype(CALOROCDigitizationConfig::toa_thres) EcalBarrelScFi_toa_thres = {4.0e-4};
  decltype(CALOROCDigitizationConfig::tot_thres) EcalBarrelScFi_tot_thres = {8.0e-4};
  decltype(CALOROCDigitizationConfig::dyRangeSingleGainADC) EcalBarrelScFi_dyRangeSingleGainADC = {
      1.0e-3};
  decltype(CALOROCDigitizationConfig::dyRangeHighGainADC) EcalBarrelScFi_dyRangeHighGainADC = {
      1.0e-3};
  decltype(CALOROCDigitizationConfig::dyRangeLowGainADC) EcalBarrelScFi_dyRangeLowGainADC = {
      1.5e-2};

  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) EcalBarrelScFi_capADC = 16384; //16384,  14bit ADC
  decltype(CalorimeterHitDigiConfig::dyRangeADC) EcalBarrelScFi_dyRangeADC   = 1500 * dd4hep::MeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) EcalBarrelScFi_pedMeanADC   = 100;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) EcalBarrelScFi_pedSigmaADC = 1;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) EcalBarrelScFi_resolutionTDC =
      10 * dd4hep::picosecond;
  ADD_TS(app, new JOmniFactoryGeneratorT<SimCalorimeterHitProcessor_factory>(
      "EcalBarrelScFiPAttenuatedHitDigi", {"EcalBarrelScFiHits"},
      {"EcalBarrelScFiPAttenuatedHitDigi", "EcalBarrelScFiPAttenuatedHitContributionDigi"},
      {
          .attenuationParameters            = EcalBarrelScFi_attPars,
          .readout                          = "EcalBarrelScFiHits",
          .attenuationReferencePositionName = "EcalBarrel_LightGuide_PositivePosZ",
          .hitMergeFields                   = EcalBarrelScFi_hitMergeFields,
          .contributionMergeFields          = EcalBarrelScFi_contributionMergeFields,
          .inversePropagationSpeed          = EcalBarrelScFi_inversePropagationSpeed,
          .fixedTimeDelay                   = EcalBarrelScFi_fixedTimeDelay,
          .timeWindow                       = EcalBarrelScFi_timeWindow,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<SimCalorimeterHitProcessor_factory>(
      "EcalBarrelScFiNAttenuatedHitDigi", {"EcalBarrelScFiHits"},
      {"EcalBarrelScFiNAttenuatedHitDigi", "EcalBarrelScFiNAttenuatedHitContributionDigi"},
      {
          .attenuationParameters            = EcalBarrelScFi_attPars,
          .readout                          = "EcalBarrelScFiHits",
          .attenuationReferencePositionName = "EcalBarrel_LightGuide_NegativePosZ",
          .hitMergeFields                   = EcalBarrelScFi_hitMergeFields,
          .contributionMergeFields          = EcalBarrelScFi_contributionMergeFields,
          .inversePropagationSpeed          = EcalBarrelScFi_inversePropagationSpeed,
          .fixedTimeDelay                   = EcalBarrelScFi_fixedTimeDelay,
          .timeWindow                       = EcalBarrelScFi_timeWindow,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<PulseGeneration_factory<edm4hep::SimCalorimeterHit>>(
      "EcalBarrelScFiPPulseDigi", {"EcalBarrelScFiPAttenuatedHitDigi"}, {"EcalBarrelScFiPPulseDigi"},
      {
          .pulse_shape_function = EcalBarrelScFi_pulse_shape_function,
          .pulse_shape_params   = EcalBarrelScFi_pulse_shape_params,
          .ignore_thres         = EcalBarrelScFi_ignore_thres,
          .timestep             = EcalBarrelScFi_timestep,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<PulseGeneration_factory<edm4hep::SimCalorimeterHit>>(
      "EcalBarrelScFiNPulseDigi", {"EcalBarrelScFiNAttenuatedHitDigi"}, {"EcalBarrelScFiNPulseDigi"},
      {
          .pulse_shape_function = EcalBarrelScFi_pulse_shape_function,
          .pulse_shape_params   = EcalBarrelScFi_pulse_shape_params,
          .ignore_thres         = EcalBarrelScFi_ignore_thres,
          .timestep             = EcalBarrelScFi_timestep,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<PulseCombiner_factory>(
      "EcalBarrelScFiPCombinedPulseDigi", {"EcalBarrelScFiPPulseDigi"}, {"EcalBarrelScFiPCombinedPulseDigi"},
      {
          .minimum_separation = EcalBarrelScFi_minimum_separation,
          .readout            = "EcalBarrelScFiHits",
          .combine_field      = EcalBarrelScFi_combine_field,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<PulseCombiner_factory>(
      "EcalBarrelScFiNCombinedPulseDigi", {"EcalBarrelScFiNPulseDigi"}, {"EcalBarrelScFiNCombinedPulseDigi"},
      {
          .minimum_separation = EcalBarrelScFi_minimum_separation,
          .readout            = "EcalBarrelScFiHits",
          .combine_field      = EcalBarrelScFi_combine_field,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<PulseNoise_factory>(
      "EcalBarrelScFiPCombinedPulsesWithNoiseDigi", {"EventHeader", "EcalBarrelScFiPCombinedPulseDigi"},
      {"EcalBarrelScFiPCombinedPulsesWithNoiseDigi"},
      {
          .poles    = EcalBarrelScFi_poles,
          .variance = EcalBarrelScFi_variance,
          .alpha    = EcalBarrelScFi_alpha,
          .scale    = EcalBarrelScFi_scale,
          .pedestal = EcalBarrelScFi_pedestal,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<PulseNoise_factory>(
      "EcalBarrelScFiNCombinedPulsesWithNoiseDigi", {"EventHeader", "EcalBarrelScFiNCombinedPulseDigi"},
      {"EcalBarrelScFiNCombinedPulsesWithNoiseDigi"},
      {
          .poles    = EcalBarrelScFi_poles,
          .variance = EcalBarrelScFi_variance,
          .alpha    = EcalBarrelScFi_alpha,
          .scale    = EcalBarrelScFi_scale,
          .pedestal = EcalBarrelScFi_pedestal,
      },
      app // TODO: Remove me once fixed
      ));
#if EDM4EIC_VERSION_MAJOR > 8 || (EDM4EIC_VERSION_MAJOR == 8 && EDM4EIC_VERSION_MINOR >= 7)
  ADD_TS(app, new JOmniFactoryGeneratorT<CALOROCDigitization_factory>(
      "EcalBarrelScFiPCALOROCHitDigi", {"EcalBarrelScFiPCombinedPulsesWithNoiseDigi"},
      {"EcalBarrelScFiPCALOROCHitDigi"},
      {
          .adc_phase            = EcalBarrelScFi_adc_phase,
          .toa_thres            = EcalBarrelScFi_toa_thres,
          .tot_thres            = EcalBarrelScFi_tot_thres,
          .dyRangeSingleGainADC = EcalBarrelScFi_dyRangeSingleGainADC,
          .dyRangeHighGainADC   = EcalBarrelScFi_dyRangeHighGainADC,
          .dyRangeLowGainADC    = EcalBarrelScFi_dyRangeLowGainADC,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CALOROCDigitization_factory>(
      "EcalBarrelScFiNCALOROCHitDigi", {"EcalBarrelScFiNCombinedPulsesWithNoiseDigi"},
      {"EcalBarrelScFiNCALOROCHitDigi"},
      {
          .adc_phase            = EcalBarrelScFi_adc_phase,
          .toa_thres            = EcalBarrelScFi_toa_thres,
          .tot_thres            = EcalBarrelScFi_tot_thres,
          .dyRangeSingleGainADC = EcalBarrelScFi_dyRangeSingleGainADC,
          .dyRangeHighGainADC   = EcalBarrelScFi_dyRangeHighGainADC,
          .dyRangeLowGainADC    = EcalBarrelScFi_dyRangeLowGainADC,
      },
      app // TODO: Remove me once fixed
      ));
#endif
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "EcalBarrelScFiRawHitDigi", {"EventHeader", "EcalBarrelScFiHits"},
      {"EcalBarrelScFiRawHitDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelScFiRawHitLinkDigi",
#endif
       "EcalBarrelScFiRawHitAssociationDigi"},
      {
          .eRes          = {0.0 * sqrt(dd4hep::GeV), 0.0, 0.0 * dd4hep::GeV},
          .tRes          = 0.0 * dd4hep::ns,
          .threshold     = 0.0 * dd4hep::keV, // threshold is set in ADC in reco
          .capADC        = EcalBarrelScFi_capADC,
          .dyRangeADC    = EcalBarrelScFi_dyRangeADC,
          .pedMeanADC    = EcalBarrelScFi_pedMeanADC,
          .pedSigmaADC   = EcalBarrelScFi_pedSigmaADC,
          .resolutionTDC = EcalBarrelScFi_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "EcalBarrelScFiHits",
          .fields        = {"fiber", "z"},
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "EcalBarrelScFiRecHitDigi", {"EcalBarrelScFiRawHitDigi"}, {"EcalBarrelScFiRecHitDigi"},
      {
          .capADC          = EcalBarrelScFi_capADC,
          .dyRangeADC      = EcalBarrelScFi_dyRangeADC,
          .pedMeanADC      = EcalBarrelScFi_pedMeanADC,
          .pedSigmaADC     = EcalBarrelScFi_pedSigmaADC, // not needed; use only thresholdValue
          .resolutionTDC   = EcalBarrelScFi_resolutionTDC,
          .thresholdFactor = 0.0, // use only thresholdValue
          .thresholdValue  = 5.0, // 16384 ADC counts/1500 MeV * 0.5 MeV (desired threshold) = 5.46
          .sampFrac        = "0.09285755",
          .readout         = "EcalBarrelScFiHits",
          .layerField      = "layer",
          .sectorField     = "sector",
          .localDetFields  = {"system", "sector"},
          // here we want to use grid center position (XY) but keeps the z information from fiber-segment
          // TODO: a more realistic way to get z is to reconstruct it from timing
          .maskPos       = "xy",
          .maskPosFields = {"fiber", "z"},
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "EcalBarrelScFiProtoClusterDigi", {"EcalBarrelScFiRecHitDigi"}, {"EcalBarrelScFiProtoClusterDigi"},
      {
          .adjacencyMatrix{},
          .peakNeighbourhoodMatrix{},
          .readout{},
          .sectorDist = 50. * dd4hep::mm,
          .localDistXY{},
          .localDistXZ = {80 * dd4hep::mm, 80 * dd4hep::mm},
          .localDistYZ{},
          .globalDistRPhi{},
          .globalDistEtaPhi{},
          .dimScaledLocalDistXY{},
          .splitCluster         = false,
          .minClusterHitEdep    = 5.0 * dd4hep::MeV,
          .minClusterCenterEdep = 100.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric{},
          .transverseEnergyProfileScale{},
          .transverseEnergyProfileScaleUnits{},
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "EcalBarrelScFiClustersWithoutShapeDigi",
      {
          "EcalBarrelScFiProtoClusterDigi", // edm4eic::ProtoClusterCollection
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
          "EcalBarrelScFiRawHitLinkDigi", // edm4eic::MCRecoCalorimeterHitLink
#endif
          "EcalBarrelScFiRawHitAssociationDigi" // edm4eic::MCRecoCalorimeterHitAssociation
      },
      {"EcalBarrelScFiClustersWithoutShapeDigi", // edm4eic::Cluster
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelScFiClusterLinksWithoutShapeDigi",
#endif
       "EcalBarrelScFiClusterAssociationsWithoutShapeDigi"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalBarrelScFiClusterDigi",
      {"EcalBarrelScFiClustersWithoutShapeDigi", "EcalBarrelScFiClusterAssociationsWithoutShapeDigi"},
      {"EcalBarrelScFiClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelScFiClusterLinkDigi",
#endif
       "EcalBarrelScFiClusterAssociationDigi"},
      {.longitudinalShowerInfoAvailable = true, .energyWeight = "log", .logWeightBase = 6.2}, app));

  // Make sure digi and reco use the same value
  decltype(SimCalorimeterHitProcessorConfig::timeWindow) EcalBarrelImaging_timeWindow = {
      100 * edm4eic::unit::ns};

  decltype(CalorimeterHitDigiConfig::capADC) EcalBarrelImaging_capADC = 8192; //8192,  13bit ADC
  decltype(CalorimeterHitDigiConfig::dyRangeADC) EcalBarrelImaging_dyRangeADC = 3 * dd4hep::MeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) EcalBarrelImaging_pedMeanADC =
      14; // Noise floor at 5 keV: 8192 / 3 * 0.005
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) EcalBarrelImaging_pedSigmaADC =
      5; // Upper limit for sigma for AstroPix
  decltype(CalorimeterHitDigiConfig::resolutionTDC) EcalBarrelImaging_resolutionTDC =
      3.25 * dd4hep::nanosecond;
  ADD_TS(app, new JOmniFactoryGeneratorT<SimCalorimeterHitProcessor_factory>(
      "EcalBarrelImagingProcessedHitDigi", {"EcalBarrelImagingHits"},
      {"EcalBarrelImagingProcessedHitDigi", "EcalBarrelImagingProcessedHitContributionDigi"},
      {
          .readout    = "EcalBarrelImagingHits",
          .timeWindow = EcalBarrelImaging_timeWindow,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "EcalBarrelImagingRawHitDigi", {"EventHeader", "EcalBarrelImagingProcessedHitDigi"},
      {"EcalBarrelImagingRawHitDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelImagingRawHitLinkDigi",
#endif
       "EcalBarrelImagingRawHitAssociationDigi"},
      {
          .eRes          = {0.0 * sqrt(dd4hep::GeV), 0.02, 0.0 * dd4hep::GeV},
          .tRes          = 0.0 * dd4hep::ns,
          .capADC        = EcalBarrelImaging_capADC,
          .dyRangeADC    = EcalBarrelImaging_dyRangeADC,
          .pedMeanADC    = EcalBarrelImaging_pedMeanADC,
          .pedSigmaADC   = EcalBarrelImaging_pedSigmaADC,
          .resolutionTDC = EcalBarrelImaging_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "EcalBarrelImagingHits",
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "EcalBarrelImagingRecHitDigi", {"EcalBarrelImagingRawHitDigi"}, {"EcalBarrelImagingRecHitDigi"},
      {
          .capADC          = EcalBarrelImaging_capADC,
          .dyRangeADC      = EcalBarrelImaging_dyRangeADC,
          .pedMeanADC      = EcalBarrelImaging_pedMeanADC,
          .pedSigmaADC     = EcalBarrelImaging_pedSigmaADC, // not needed; use only thresholdValue
          .resolutionTDC   = EcalBarrelImaging_resolutionTDC,
          .thresholdFactor = 0.0, // use only thresholdValue
          .thresholdValue  = 41,  // 8192 ADC counts/3 MeV * 0.015 MeV (desired threshold) = 41
          .sampFrac        = "0.00429453",
          .readout         = "EcalBarrelImagingHits",
          .layerField      = "layer",
          .sectorField     = "sector",
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<ImagingTopoCluster_factory>(
      "EcalBarrelImagingProtoClusterDigi", {"EcalBarrelImagingRecHitDigi"},
      {"EcalBarrelImagingProtoClusterDigi"},
      {
          .neighbourLayersRange = 2, //  # id diff for adjacent layer
          .sameLayerDistTZ      = {2.0 * dd4hep::mm, 2 * dd4hep::mm},     //  # same layer
          .diffLayerDistEtaPhi  = {10 * dd4hep::mrad, 10 * dd4hep::mrad}, //  # adjacent layer
          .sameLayerMode        = eicrecon::ImagingTopoClusterConfig::ELayerMode::tz,
          .diffLayerMode        = eicrecon::ImagingTopoClusterConfig::ELayerMode::etaphi,
          .sectorDist           = 3.0 * dd4hep::cm,
          .minClusterHitEdep    = 0,
          .minClusterCenterEdep = 0,
          .minClusterEdep       = 100 * dd4hep::MeV,
          .minClusterNhits      = 10,
      },
      app // TODO: Remove me once fixed
      ));

  ADD_TS(app, new JOmniFactoryGeneratorT<ImagingClusterReco_factory>(
      "EcalBarrelImagingClustersWithoutShapeDigi",
      {"EcalBarrelImagingProtoClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelImagingRawHitLinkDigi",
#endif
       "EcalBarrelImagingRawHitAssociationDigi"},
      {"EcalBarrelImagingClustersWithoutShapeDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelImagingClusterLinksWithoutShapeDigi",
#endif
       "EcalBarrelImagingClusterAssociationsWithoutShapeDigi", "EcalBarrelImagingLayerDigi"},
      {
          .trackStopLayer = 6,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "EcalBarrelImagingClusterDigi",
      {"EcalBarrelImagingClustersWithoutShapeDigi",
       "EcalBarrelImagingClusterAssociationsWithoutShapeDigi"},
      {"EcalBarrelImagingClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelImagingClusterLinkDigi",
#endif
       "EcalBarrelImagingClusterAssociationDigi"},
      {.longitudinalShowerInfoAvailable = false, .energyWeight = "log", .logWeightBase = 6.2},
      app));
  ADD_TS(app, new JOmniFactoryGeneratorT<EnergyPositionClusterMerger_factory>(
      "EcalBarrelClusterDigi",
      {"EcalBarrelScFiClusterDigi", "EcalBarrelScFiClusterAssociationDigi", "EcalBarrelImagingClusterDigi",
       "EcalBarrelImagingClusterAssociationDigi"},
      {"EcalBarrelClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelClusterLinkDigi",
#endif
       "EcalBarrelClusterAssociationDigi"},
      {
          .energyRelTolerance = 0.5,
          .phiTolerance       = 0.1,
          .etaTolerance       = 0.2,
      },
      app // TODO: Remove me once fixed
      ));
  ADD_TS(app, new JOmniFactoryGeneratorT<TruthEnergyPositionClusterMerger_factory>(
      "EcalBarrelTruthClusterDigi",
      {"MCParticles", "EcalBarrelScFiClusterDigi", "EcalBarrelScFiClusterAssociationDigi",
       "EcalBarrelImagingClusterDigi", "EcalBarrelImagingClusterAssociationDigi"},
      {"EcalBarrelTruthClusterDigi",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
       "EcalBarrelTruthClusterLinkDigi",
#endif
       "EcalBarrelTruthClusterAssociationDigi"},
      app // TODO: Remove me once fixed
      ));
}
// }
