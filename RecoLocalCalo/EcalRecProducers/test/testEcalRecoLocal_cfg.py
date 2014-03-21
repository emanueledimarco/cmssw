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
# methods: Weights, AlphaBetaFit, AlphaBetaGammaFit, AnalyticFit
AmplitudeRecoMethod = "AlphaBetaGammaFit"
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

# get uncalibrechits with weights method
import RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi
process.ecalUncalibHitWeights = RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi.ecalWeightUncalibRecHit.clone()
    
# get uncalibrechits with alphabeta fit method
import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi
process.ecalUncalibHitFixedAlphaBetaFit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi.ecalFixedAlphaBetaFitUncalibRecHit.clone()

# get uncalibrechits with alphabetagamma fit method
import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaGammaFitUncalibRecHit_cfi
process.ecalUncalibHitFixedAlphaBetaGammaFit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaGammaFitUncalibRecHit_cfi.ecalFixedAlphaBetaGammaFitUncalibRecHit.clone()

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
if AmplitudeRecoMethod == "Weights":
    print "==> Using Weights as amplitude reconstruction"
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitWeights:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitWeights:EcalUncalibRecHitsEE'
elif AmplitudeRecoMethod == "AlphaBetaFit":
    print "==> Using AlphaBetaFit as amplitude reconstruction" 
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitFixedAlphaBetaFit:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitFixedAlphaBetaFit:EcalUncalibRecHitsEE'
elif AmplitudeRecoMethod == "AlphaBetaGammaFit":
    print "==> Using AlphaBetaGammaFit as amplitude reconstruction" 
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitFixedAlphaBetaGammaFit:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitFixedAlphaBetaGammaFit:EcalUncalibRecHitsEE'    
elif AmplitudeRecoMethod == "AnalyticFit":
    print "==> Using AnalyticFit as amplitude reconstruction" 
    process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHitAnalFit:EcalUncalibRecHitsEB'
    process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHitAnalFit:EcalUncalibRecHitsEE'        
else: print wrongAlgoMessage

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",
#            fileNames = cms.untracked.vstring('/store/group/phys_egamma/emanuele/ecal/reconstruction/DYToEE_M_20_TuneZ2star_8TeV_pythia6_GEN-SIM-RAW+PU25bx25_START53_V19D-v1.root')
             fileNames = cms.untracked.vstring('file:/cmsrm/pc24_2/emanuele/data/ZEE_14TeV_DR61SLHCx_GENSIMRECO.root')
                )


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_ecalUncalibHit*_*_*',
                                                                      'keep *_ecalRecHit_*_*'
                                                                      ),
                               fileName = cms.untracked.string('testEcalLocalRecoA.root')
                               )



if AmplitudeRecoMethod == "Weights":             process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitWeights)
elif AmplitudeRecoMethod == "AlphaBetaFit":      process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitFixedAlphaBetaFit)
elif AmplitudeRecoMethod == "AlphaBetaGammaFit": process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitFixedAlphaBetaGammaFit)
elif AmplitudeRecoMethod == "AnalyticFit":       process.ecalAmplitudeReco = cms.Sequence(process.ecalUncalibHitAnalFit)
else: print wrongAlgoMessage

process.ecalTestRecoLocal = cms.Sequence(process.ecalEBunpacker
                                         *process.ecalAmplitudeReco
                                         *process.ecalUncalibHitGlobal
#                                         *process.ecalDetIdToBeRecovered
                                         *process.ecalRecHit
                                        )

from PhysicsTools.PatAlgos.tools.helpers import *
if isMC:
     massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:ebDigis"), cms.InputTag("simEcalDigis:ebDigis"),True)
     massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:eeDigis"), cms.InputTag("simEcalDigis:eeDigis"),True)
if runOnRAW:
    massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:ebDigis"), cms.InputTag("ecalEBunpacker:ebDigis"),True)
    massSearchReplaceAnyInputTag(process.ecalTestRecoLocal,cms.InputTag("ecalDigis:eeDigis"), cms.InputTag("ecalEBunpacker:eeDigis"),True)

process.p = cms.Path(process.ecalTestRecoLocal)
process.outpath = cms.EndPath(process.out)
process.GlobalTag.globaltag = 'START17_61_V5::All'
