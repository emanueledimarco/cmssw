import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from glob import glob
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'file:/tmp/gpetrucc/CMSSW_5_2_4_patch1/src/WWAnalysis/SkimStep/test/latinosYieldSkim.SingleMu.root',
        'file:/tmp/gpetrucc/CMSSW_5_2_4_patch1/src/WWAnalysis/SkimStep/test/latinosYieldSkim.DoubleEl.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("WWAnalysis.AnalysisStep.tp_lepton_effs_cff")

process.runMuMu = cms.Path(process.tnpSimpleSequenceMuMu)
process.runElEl = cms.Path(process.tnpSimpleSequenceElEl) 

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpLeptonID.root"))

