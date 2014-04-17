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
# methods: Global, Weights, AlphaBetaFit, AnalyticFit
AmplitudeRecoMethod = "Global"
runOnRAW = False
runProfiler = False
#####################

if runProfiler:
    process.IgProfService = cms.Service("IgProfService",
                                        reportFirstEvent            = cms.untracked.int32(0),
                                        reportEventInterval         = cms.untracked.int32(1),
                                        reportToFileAtPostEvent     = cms.untracked.string("| gzip -c > XXXX.%I.gz")
                                        )

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

# get uncalibrechits with ratio method (does not run in 5_3_12)
#import RecoLocalCalo.EcalRecProducers.ecalRatioUncalibRecHit_cfi
#process.ecalUncalibHitRatio = RecoLocalCalo.EcalRecProducers.ecalRatioUncalibRecHit_cfi.ecalRatioUncalibRecHit.clone()
#process.ecalUncalibHitRatio.EBdigiCollection = 'ecalEBunpacker:ebDigis'
#process.ecalUncalibHitRatio.EEdigiCollection = 'ecalEBunpacker:eeDigis'


# get uncalibrechits with ratio method
import RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi
process.ecalUncalibHitGlobal = RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi.ecalGlobalUncalibRecHit.clone()

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

wrongAlgoMessage = "The AmplitudeRecoMethod: "+AmplitudeRecoMethod+" that you chose is not foreseen. Check the customization of the cfg "

# get rechits e.g. from the weights
process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")
process.ecalRecHit.triggerPrimitiveDigiCollection = 'ecalEBunpacker:EcalTriggerPrimitives'
if AmplitudeRecoMethod == "Global":
    print "==> Using Global as amplitude reconstruction"
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitGlobal:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitGlobal:EcalUncalibRecHitsEE'
elif AmplitudeRecoMethod == "Weights":
    print "==> Using Weights as amplitude reconstruction"
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitWeights:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitWeights:EcalUncalibRecHitsEE'
elif AmplitudeRecoMethod == "AlphaBetaFit":
    print "==> Using AlphaBetaFit as amplitude reconstruction" 
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitFixedAlphaBetaFit:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitFixedAlphaBetaFit:EcalUncalibRecHitsEE'
elif AmplitudeRecoMethod == "AnalyticFit":
    print "==> Using AnalyticFit as amplitude reconstruction" 
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitAnalFit:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitAnalFit:EcalUncalibRecHitsEE'        
else: print wrongAlgoMessage

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
#             fileNames = cms.untracked.vstring('/store/cmst3/user/bendavid/photongun_nopu/photongun_nopu_1000_1_POT.root')
              fileNames = cms.untracked.vstring('/store/cmst3/user/bendavid/photongun_pu25/photongun_pu25_988_1_AV6.root')
                ) 


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_ecalUncalibHit*_*_testEcalRecoLocal',
                                                                      'keep *_ecalRecHit_*_testEcalRecoLocal',
                                                                      'keep *_offlineBeamSpot_*_*',
                                                                      'keep *_addPileupInfo_*_*'
                                                                      ),
                               fileName = cms.untracked.string('testEcalLocalRecoA.root')
                               )



if AmplitudeRecoMethod == "Global":              process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitGlobal)
elif AmplitudeRecoMethod == "Weights":           process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitWeights)
elif AmplitudeRecoMethod == "AlphaBetaFit":      process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitFixedAlphaBetaFit)
elif AmplitudeRecoMethod == "AnalyticFit":       process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitAnalFit)
else: print wrongAlgoMessage

process.ecalTestRecoLocal = cms.Sequence(#process.ecalEBunpacker
                                         process.ecalAmplitudeReco
                                         *process.ecalUncalibHitGlobal
#                                         *process.ecalDetIdToBeRecovered
                                         *process.ecalRecHit
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
