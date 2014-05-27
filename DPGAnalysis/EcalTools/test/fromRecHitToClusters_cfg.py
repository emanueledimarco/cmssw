import FWCore.ParameterSet.Config as cms
process = cms.Process("HitsToClusters")

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
#####################

# get timing service up for profiling
process.TimerService = cms.Service("TimerService")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# get uncalib rechits from global method (weights for ampli, time from ratio, etc)
import RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi
process.ecalGlobalUncalibRecHit.subtractPU = cms.bool(True)

# get calibrated rechits
process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")

# get the recovered digis
if isMC:
    process.ecalDetIdToBeRecovered.ebSrFlagCollection = 'simEcalDigis:ebSrFlags'
    process.ecalDetIdToBeRecovered.eeSrFlagCollection = 'simEcalDigis:eeSrFlags'
    process.ecalRecHit.recoverEBFE = False
    process.ecalRecHit.recoverEEFE = False
    process.ecalRecHit.killDeadChannels = False
    process.ecalRecHit.ebDetIdToBeRecovered = ''
    process.ecalRecHit.eeDetIdToBeRecovered = ''
    process.ecalRecHit.ebFEToBeRecovered = ''
    process.ecalRecHit.eeFEToBeRecovered = ''


# and now the clusters and superclusters
process.load("RecoEcal.Configuration.RecoEcal_cff")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff")
process.particleFlowSuperClusterECAL.useRegression = False

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(100) )
process.source = cms.Source("PoolSource",
              fileNames = cms.untracked.vstring('/store/cmst3/user/bendavid/photongun_pu25/photongun_pu25_988_1_AV6.root')
                ) 


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*_*_HitsToClusters',
                                                                      'keep *_offlineBeamSpot_*_*',
                                                                      'keep *_addPileupInfo_*_*',
                                                                      'keep *_genParticles_*_*'
                                                                      ),
                               fileName = cms.untracked.string('testEcalFullReco.root')
                               )


process.ecalLocalReco = cms.Sequence( process.ecalGlobalUncalibRecHit
                                      * process.ecalRecHit )

process.clusteringSequence = cms.Sequence( process.pfClusteringECAL
                                           * process.ecalClusters )

process.p = cms.Path( process.ecalLocalReco * process.clusteringSequence )
process.outpath = cms.EndPath(process.out)

process.GlobalTag.globaltag = 'START71_V1::All'

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("EcalTBWeightsRcd"),
             tag = cms.string("EcalTBWeights_3p5_time_mc"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
             )
    )

