import FWCore.ParameterSet.Config as cms

process = cms.Process("HWWJetCleanerTest")

process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWJetCleaner_cfi")
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWJetCleanerRef_cfi")

process.options = cms.untracked.PSet(
    fileMode =  cms.untracked.string('NOMERGE')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/e/emanuele/RECO/relval_Zee_CMSSW_2_1_7.root')
                            )

process.out = cms.OutputModule("PoolOutputModule",
                               verbose = cms.untracked.bool(False),
                               fileName = cms.untracked.string('JetCleaningTestOutput.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*_*_HWWJetCleanerTest',
                                                                      'keep recoCaloJets_*_*_*'
                                                                      )
                               )

process.p = cms.Path ( process.cleanedJets * process.cleanedJetsRef )
process.o = cms.EndPath( process.out )
