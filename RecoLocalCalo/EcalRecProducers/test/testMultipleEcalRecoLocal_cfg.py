import FWCore.ParameterSet.Config as cms
process = cms.Process("RECO2")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')

#### CONFIGURE IT HERE
isMC = True
#####################
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# start from RAW format for more flexibility
process.raw2digi_step = cms.Sequence(process.RawToDigi)

# get uncalibrechits with global method / time from ratio
import RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi
process.ecalGlobalUncalibRecHit = RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi.ecalGlobalUncalibRecHit.clone()
# get uncalib rechits from multifit method / time from ratio
import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi
process.ecalMultiFitUncalibRecHit =  RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHit.algoPSet.activeBXs = cms.vint32(-5,-4,-3,-2,-1,0,1,2,3,4)
process.ecalMultiFitUncalibRecHit.algoPSet.useLumiInfoRunHeader = cms.bool( False )
# get uncalib rechits from multifit method / time from weights
process.ecalMultiFit2UncalibRecHit =  RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFit2UncalibRecHit.algoPSet.timealgo = cms.string("WeightsMethod")
process.ecalMultiFit2UncalibRecHit.algoPSet.activeBXs = cms.vint32(-5,-4,-3,-2,-1,0,1,2,3,4)
process.ecalMultiFit2UncalibRecHit.algoPSet.useLumiInfoRunHeader = cms.bool ( False )

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


process.ecalRecHitGlobal = process.ecalRecHit.clone()
process.ecalRecHitGlobal.EBuncalibRecHitCollection = 'ecalGlobalUncalibRecHit:EcalUncalibRecHitsEB'
process.ecalRecHitGlobal.EEuncalibRecHitCollection = 'ecalGlobalUncalibRecHit:EcalUncalibRecHitsEE'
process.ecalRecHitGlobal.EBrechitCollection = 'EcalRecHitsGlobalEB'
process.ecalRecHitGlobal.EErechitCollection = 'EcalRecHitsGlobalEE'

process.ecalRecHitMultiFit = process.ecalRecHit.clone()
process.ecalRecHitMultiFit.EBuncalibRecHitCollection = 'ecalMultiFitUncalibRecHit:EcalUncalibRecHitsEB'
process.ecalRecHitMultiFit.EEuncalibRecHitCollection = 'ecalMultiFitUncalibRecHit:EcalUncalibRecHitsEE'
process.ecalRecHitMultiFit.EBrechitCollection = 'EcalRecHitsMultiFitEB'
process.ecalRecHitMultiFit.EErechitCollection = 'EcalRecHitsMultiFitEE'

process.ecalRecHitMultiFit2 = process.ecalRecHit.clone()
process.ecalRecHitMultiFit2.EBuncalibRecHitCollection = 'ecalMultiFit2UncalibRecHit:EcalUncalibRecHitsEB'
process.ecalRecHitMultiFit2.EEuncalibRecHitCollection = 'ecalMultiFit2UncalibRecHit:EcalUncalibRecHitsEE'
process.ecalRecHitMultiFit2.EBrechitCollection = 'EcalRecHitsMultiFit2EB'
process.ecalRecHitMultiFit2.EErechitCollection = 'EcalRecHitsMultiFit2EE'

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(1000) )
path = '/store/data/Run2012D/DoubleElectron/RAW-RECO/ZElectron-22Jan2013-v1/10000/'
process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(path+'0008202C-E78F-E211-AADB-0026189437FD.root'
                                                              ))


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_ecalUncalib*_*_RECO2',
                                                                      'keep *_ecalRecHit*_*_RECO2',
                                                                      'keep *_offlineBeamSpot_*_*',
                                                                      'keep *_addPileupInfo_*_*'
                                                                      ),
                               fileName = cms.untracked.string('reco2_pu40.root')
                               )


process.ecalAmplitudeReco = cms.Sequence( process.ecalGlobalUncalibRecHit *
                                          process.ecalMultiFitUncalibRecHit *
                                          process.ecalMultiFit2UncalibRecHit
                                          )

process.ecalRecHitsReco = cms.Sequence( process.ecalRecHitGlobal *
                                        process.ecalRecHitMultiFit *
                                        process.ecalRecHitMultiFit2 )

process.ecalTestRecoLocal = cms.Sequence( process.raw2digi_step *
                                          process.ecalAmplitudeReco *
                                          process.ecalRecHitsReco )

from PhysicsTools.PatAlgos.tools.helpers import *

process.p = cms.Path(process.ecalTestRecoLocal)
process.outpath = cms.EndPath(process.out)


