import FWCore.ParameterSet.Config as cms
process = cms.Process("piZeroAnalysis")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("DQMOffline.EGamma.piZeroAnalyzer_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cff")
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

DQMStore = cms.Service("DQMStore")

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        '/store/data/Run2015B/AlCaP0/RAW/v1/000/251/244/00000/86E14A79-AE25-E511-936E-02163E01207C.root'
#        '/store/relval/CMSSW_3_6_0_pre2/RelValSingleGammaPt35/GEN-SIM-RECO/MC_3XY_V24-v1/0001/364E7B38-6F27-DF11-91A9-0026189438D4.root',
#       '/store/relval/CMSSW_3_6_0_pre2/RelValSingleGammaPt35/GEN-SIM-RECO/MC_3XY_V24-v1/0000/48AE643B-0727-DF11-99FB-001731AF66F5.root'


       

##    '/store/relval/CMSSW_3_0_0_pre6/RelValSingleGammaPt35/GEN-SIM-RECO/IDEAL_30X_v1/0005/98C45436-41DE-DD11-9B91-000423D95220.root'


##         '/store/relval/CMSSW_3_0_0_pre6/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP_30X_v1/0005/60FEA65F-2DDE-DD11-A12E-001617E30F56.root',
##         '/store/relval/CMSSW_3_0_0_pre6/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP_30X_v1/0005/C0E99B59-34DE-DD11-9194-000423D174FE.root',
##         '/store/relval/CMSSW_3_0_0_pre6/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP_30X_v1/0005/F41CD6A5-41DE-DD11-8049-000423D99AAA.root',
##         '/store/relval/CMSSW_3_0_0_pre6/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP_30X_v1/0005/FCFD0B48-2BDE-DD11-97C8-000423D99B3E.root'



))



process.FEVT = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring("keep *_MEtoEDMConverter_*_*"),
    fileName = cms.untracked.string('photonsMEtoEDMConverter.root')
)

# start from RAW format for more flexibility
#process.raw2digi_step = cms.Sequence(process.RawToDigi)

#DUMMY RECHIT
process.dummyHits = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(True),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB' ,'HLT'),
                                     endcapHitProducer      = cms.InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE' ,'HLT'),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis","HLT"),
                                     endcapDigis            = cms.InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis","HLT"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))

# ECAL reco
import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi
process.ecalMultiFitUncalibRecHit =  RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag('dummyHits','dummyBarrelDigis','piZeroAnalysis')
process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag('dummyHits','dummyEndcapDigis','piZeroAnalysis')

#UNCALIB to CALIB
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *
process.ecalDetIdToBeRecovered =  RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi.ecalDetIdToBeRecovered.clone()
process.ecalRecHit.killDeadChannels = cms.bool( False )
process.ecalRecHit.recoverEBVFE = cms.bool( False )
process.ecalRecHit.recoverEEVFE = cms.bool( False )
process.ecalRecHit.recoverEBFE = cms.bool( False )
process.ecalRecHit.recoverEEFE = cms.bool( False )
process.ecalRecHit.recoverEEIsolatedChannels = cms.bool( False )
process.ecalRecHit.recoverEBIsolatedChannels = cms.bool( False )

process.ecalReco = cms.Sequence(process.dummyHits
                                * process.ecalMultiFitUncalibRecHit 
                                * process.ecalRecHit )

#process.options = cms.untracked.PSet(
#   wantSummary = cms.untracked.bool(True),
#   SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError','CrystalIDError2')
#)

from DQMOffline.EGamma.piZeroAnalyzer_cfi import *
piZeroAnalysis.OutputMEsInRootFile = cms.bool(True)
piZeroAnalysis.OutputFileName = 'DQMPiZeros.root'
piZeroAnalysis.Verbosity = cms.untracked.int32(0)
piZeroAnalysis.standAlone = cms.bool(True)


process.load('DQMServices.Components.DQMFileSaver_cfi')
process.dqmSaver.workflow = '/Egamma/PiZeroAnalyzer/All'

#process.p1 = cms.Path(process.MEtoEDMConverter)
#process.p1 = cms.Path(process.piZeroAnalysis*process.MEtoEDMConverter)
process.p1 = cms.Path(process.ecalReco * process.piZeroAnalysis * process.dqmSaver)
process.schedule = cms.Schedule(process.p1)

