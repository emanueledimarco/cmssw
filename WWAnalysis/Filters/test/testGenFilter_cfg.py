import FWCore.ParameterSet.Config as cms

process = cms.Process("Peaking")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.GlobalTag.globaltag = 'START38_V12::All'

process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring())
process.source.fileNames.extend(['file:/nfs/bluearc/group/trees/hww/R42X_S1_V04_S2_V01_S3_V00/001.WZtoAny/WZtoAny_1_1_X2C.root' ])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 


process.genFilter = cms.EDFilter("GenFilterDiBosons")
process.p = cms.Path(process.genFilter)
process.f = cms.Path(~process.genFilter)

# ---- endPath ----
process.out1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("peaking.root"),
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
)
process.out2 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("nonPeaking.root"),
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
)
process.end = cms.EndPath(process.out1+process.out2)


