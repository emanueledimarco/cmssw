import FWCore.ParameterSet.Config as cms
process = cms.Process("testEcalRecoLocal")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoEgamma.EgammaElectronProducers.gsfElectronSequence_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi")
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")

#### CONFIGURE IT HERE
isMC = True
runOnRAW = False
useTrivial = True
#####################

if isMC:
    process.ecalEBunpacker.InputLabel = cms.InputTag('rawDataCollector')
else:
    process.ecalEBunpacker.InputLabel = cms.InputTag('source')

# get timing service up for profiling
process.TimerService = cms.Service("TimerService")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# get uncalib rechits from global method (weights for ampli, time from ratio, etc)
import RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi
process.ecalUncalibHitGlobal = RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi.ecalGlobalUncalibRecHit.clone()

# get uncalibrechits with weights method
import RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi
process.ecalUncalibHitWeights = RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi.ecalWeightUncalibRecHit.clone()
    
# get uncalibrechits with alphabeta fit method
import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi
process.ecalUncalibHitFixedAlphaBetaFit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi.ecalFixedAlphaBetaFitUncalibRecHit.clone()

# get uncalibrechits with analytic fit method  
import RecoLocalCalo.EcalRecProducers.ecalAnalFitUncalibRecHit_cfi
process.ecalUncalibHitAnalFit = RecoLocalCalo.EcalRecProducers.ecalAnalFitUncalibRecHit_cfi.ecalAnalFitUncalibRecHit.clone()

# get uncalibrechits with ratio method
import RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi
process.ecalUncalibHitGlobal = RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi.ecalGlobalUncalibRecHit.clone()

# get rechits e.g. from the weights
process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")

process.ecalRecHit.triggerPrimitiveDigiCollection = 'ecalEBunpacker:EcalTriggerPrimitives'

# get the recovered digis
if isMC:
    process.ecalDetIdToBeRecovered.ebSrFlagCollection = 'simEcalDigis:ebSrFlags'
    process.ecalDetIdToBeRecovered.eeSrFlagCollection = 'simEcalDigis:eeSrFlags'
    process.ecalRecHit.recoverEBFE = False
    process.ecalRecHit.recoverEEFE = False
    process.ecalRecHit.killDeadChannels = False
if runOnRAW:
    process.ecalDetIdToBeRecovered.ebIntegrityGainErrors = 'ecalEBunpacker:EcalIntegrityGainErrors'
    process.ecalDetIdToBeRecovered.ebIntegrityGainSwitchErrors = 'ecalEBunpacker:EcalIntegrityGainSwitchErrors'
    process.ecalDetIdToBeRecovered.ebIntegrityChIdErrors = 'ecalEBunpacker:EcalIntegrityChIdErrors'
    process.ecalDetIdToBeRecovered.eeIntegrityGainErrors = 'ecalEBunpacker:EcalIntegrityGainErrors'
    process.ecalDetIdToBeRecovered.eeIntegrityGainSwitchErrors = 'ecalEBunpacker:EcalIntegrityGainSwitchErrors'
    process.ecalDetIdToBeRecovered.eeIntegrityChIdErrors = 'ecalEBunpacker:EcalIntegrityChIdErrors'
    process.ecalDetIdToBeRecovered.integrityTTIdErrors = 'ecalEBunpacker:EcalIntegrityTTIdErrors'
    process.ecalDetIdToBeRecovered.integrityBlockSizeErrors = 'ecalEBunpacker:EcalIntegrityBlockSizeErrors'
else:
    process.ecalRecHit.ebDetIdToBeRecovered = ''
    process.ecalRecHit.eeDetIdToBeRecovered = ''
    process.ecalRecHit.ebFEToBeRecovered = ''
    process.ecalRecHit.eeFEToBeRecovered = ''


process.ecalRecHitGlobal = process.ecalRecHit.clone()
process.ecalRecHitGlobal.EBuncalibRecHitCollection = 'ecalUncalibHitGlobal:EcalUncalibRecHitsEB'
process.ecalRecHitGlobal.EEuncalibRecHitCollection = 'ecalUncalibHitGlobal:EcalUncalibRecHitsEE'
process.ecalRecHitGlobal.EBrechitCollection = 'EcalRecHitsGlobalEB'
process.ecalRecHitGlobal.EErechitCollection = 'EcalRecHitsGlobalEE'

process.ecalRecHitWeights = process.ecalRecHit.clone()
process.ecalRecHitWeights.EBuncalibRecHitCollection = 'ecalUncalibHitWeights:EcalUncalibRecHitsEB'
process.ecalRecHitWeights.EEuncalibRecHitCollection = 'ecalUncalibHitWeights:EcalUncalibRecHitsEE'
process.ecalRecHitWeights.EBrechitCollection = 'EcalRecHitsWeightsEB'
process.ecalRecHitWeights.EErechitCollection = 'EcalRecHitsWeightsEE'

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
              fileNames = cms.untracked.vstring('/store/cmst3/user/bendavid/photongun_pu25/photongun_pu25_988_1_AV6.root')
                ) 


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_ecalUncalibHit*_*_testEcalRecoLocal',
                                                                      'keep *_ecalRecHit*_*_testEcalRecoLocal',
                                                                      'keep *_offlineBeamSpot_*_*',
                                                                      'keep *_addPileupInfo_*_*'
                                                                      ),
                               fileName = cms.untracked.string('testEcalLocalRecoA.root')
                               )


process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitGlobal*
                                         process.ecalUncalibHitWeights)

process.ecalRecHitsReco = cms.Sequence(process.ecalRecHitGlobal
                                       *process.ecalRecHitWeights)

process.ecalTestRecoLocal = cms.Sequence(#process.ecalEBunpacker
                                         process.ecalAmplitudeReco
                                         *process.ecalRecHitsReco
                                         )

from PhysicsTools.PatAlgos.tools.helpers import *
#if isMC:
#     massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:ebDigis"), cms.InputTag("simEcalDigis:ebDigis"),True)
#     massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:eeDigis"), cms.InputTag("simEcalDigis:eeDigis"),True)
if runOnRAW:
    massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:ebDigis"), cms.InputTag("ecalEBunpacker:ebDigis"),True)
    massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:eeDigis"), cms.InputTag("ecalEBunpacker:eeDigis"),True)

process.p = cms.Path(process.ecalTestRecoLocal)
process.outpath = cms.EndPath(process.out)


process.GlobalTag.globaltag = 'START71_V1::All'
# if useTrivial:
#     process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")
#     process.es_prefer_EcalTrivialConditionRetriever = cms.ESPrefer("EcalTrivialConditionRetriever")

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("EcalTBWeightsRcd"),
             tag = cms.string("EcalTBWeights_3p5_time_mc"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
             )
    )
