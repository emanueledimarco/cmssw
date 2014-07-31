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


# get uncalib rechits from global method (weights for ampli, time from ratio, etc)
import RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi
# modify the rechits sequence activating the PU sub
process.ecalGlobalUncalibRecHitNoPU = process.ecalGlobalUncalibRecHit.clone()
process.ecalGlobalUncalibRecHitNoPU.subtractPU = cms.bool(True)

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
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.source = cms.Source("PoolSource",
              fileNames = cms.untracked.vstring('/store/cmst3/user/bendavid/photongun_pu25/photongun_pu25_988_1_AV6.root')
                ) 


process.ecalRecHitNoPU = process.ecalRecHit.clone()
process.ecalRecHitNoPU.EBuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHitNoPU","EcalUncalibRecHitsEB")
process.ecalRecHitNoPU.EEuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHitNoPU","EcalUncalibRecHitsEE")

process.ecalLocalReco = cms.Sequence( process.ecalGlobalUncalibRecHit * process.ecalRecHit
                                      * process.ecalGlobalUncalibRecHitNoPU * process.ecalRecHitNoPU )


# modify the clustering sequence to use the PU-subtracted rechits

# barrel:
process.load("RecoEcal.EgammaClusterProducers.correctedHybridSuperClusters_cfi")
process.correctedHybridSuperClustersNoPU = process.correctedHybridSuperClusters.clone()
process.correctedHybridSuperClustersNoPU.recHitProducer = cms.InputTag("ecalRecHitNoPU","EcalRecHitsEB")

# endcap:
process.load("RecoEcal.EgammaClusterProducers.multi5x5BasicClusters_cfi")
process.multi5x5BasicClustersCleanedNoPU = process.multi5x5BasicClustersCleaned.clone()
process.multi5x5BasicClustersCleanedNoPU.endcapHitTag = cms.InputTag('ecalRecHitNoPU','EcalRecHitsEE')
process.load("RecoEcal.EgammaClusterProducers.multi5x5SuperClusters_cfi")
process.multi5x5SuperClustersCleanedNoPU = process.multi5x5SuperClustersCleaned.clone()
process.multi5x5SuperClustersCleanedNoPU.endcapClusterTag = cms.InputTag('multi5x5BasicClustersCleanedNoPU', 'multi5x5EndcapBasicClusters')
process.load("RecoEcal.EgammaClusterProducers.multi5x5SuperClustersWithPreshower_cfi")
process.multi5x5SuperClustersWithPreshowerNoPU = process.multi5x5SuperClustersWithPreshower.clone()
endcapSClusterProducer = cms.InputTag("multi5x5SuperClustersNoPU","multi5x5EndcapSuperClusters")
process.correctedMulti5x5SuperClustersWithPreshowerNoPU = process.correctedMulti5x5SuperClustersWithPreshower.clone()
process.correctedMulti5x5SuperClustersWithPreshowerNoPU.rawSuperClusterProducer = cms.InputTag("multi5x5SuperClustersWithPreshowerNoPU")
process.correctedMulti5x5SuperClustersWithPreshowerNoPU.recHitProducer = cms.InputTag("ecalRecHitNoPU","EcalRecHitsEE")

process.ClusterECALSequence = cms.Sequence(process.correctedHybridSuperClusters
                                           * process.multi5x5BasicClustersCleaned
                                           * process.multi5x5SuperClustersCleaned
                                           * process.multi5x5SuperClustersWithPreshower
                                           * process.correctedMulti5x5SuperClustersWithPreshower )

process.ClusterECALSequenceNoPU = cms.Sequence(process.correctedHybridSuperClustersNoPU
                                               * process.multi5x5BasicClustersCleanedNoPU
                                               * process.multi5x5SuperClustersCleanedNoPU
                                               * process.multi5x5SuperClustersWithPreshowerNoPU
                                               * process.correctedMulti5x5SuperClustersWithPreshowerNoPU )

process.ClusteringSequence = cms.Sequence(process.ClusterECALSequence # the one with the new hits, no PU-sub
                                          * process.ClusterECALSequenceNoPU # the one with the new hits, yes PU-sub
                                          )


# PFlow clustering (barrel and endcap in one go):
process.load("RecoParticleFlow/PFClusterProducer/particleFlowRecHitECALNoPU_cfi")
process.particleFlowClusterECALUncorrectedNoPU = process.particleFlowClusterECALUncorrected.clone()
process.particleFlowClusterECALUncorrectedNoPU.recHitsSource = cms.InputTag("particleFlowRecHitECALNoPU")
process.particleFlowClusterECALNoPU = process.particleFlowClusterECAL.clone()
process.particleFlowClusterECALNoPU.inputECAL = cms.InputTag("particleFlowClusterECALUncorrectedNoPU")
process.load("RecoEcal.EgammaClusterProducers.particleFlowSuperClusterECAL_cfi")
process.particleFlowSuperClusterECALNoPU = process.particleFlowSuperClusterECAL.clone()
process.particleFlowSuperClusterECALNoPU.PFClusters = cms.InputTag("particleFlowClusterECALNoPU")
process.particleFlowSuperClusterECALNoPU.ESAssociation = cms.InputTag("particleFlowClusterECALNoPU")


process.particleFlowClusterECALSequence = cms.Sequence(process.particleFlowRecHitECAL
                                                       * process.particleFlowClusterECALUncorrected * process.particleFlowClusterECAL
                                                       * process.particleFlowSuperClusterECAL )

process.particleFlowClusterECALSequenceNoPU = cms.Sequence(process.particleFlowRecHitECALNoPU
                                                           * process.particleFlowClusterECALUncorrectedNoPU * process.particleFlowClusterECALNoPU
                                                           * process.particleFlowSuperClusterECALNoPU )

process.particleFlowClusteringSequence = cms.Sequence(process.particleFlowClusterECALSequence # the one with the new hits, no PU-sub
                                                      * process.particleFlowClusterECALSequenceNoPU # the one with the new hits, yes PU-sub
                                                      )

process.p = cms.Path( process.ecalLocalReco
                      * process.ClusteringSequence
                      * process.particleFlowClusteringSequence )


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*SuperClusters*_*_*',
                                                                      'keep *_offlineBeamSpot_*_*',
                                                                      'keep *_addPileupInfo_*_*',
                                                                      'keep *_genParticles_*_*'
                                                                      ),
                               fileName = cms.untracked.string('testEcalFullReco.root')
                               )

process.outpath = cms.EndPath(process.out)



process.GlobalTag.globaltag = 'START71_V1::All'

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("EcalTBWeightsRcd"),
             tag = cms.string("EcalTBWeights_3p5_time_mc"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
             )
    )

