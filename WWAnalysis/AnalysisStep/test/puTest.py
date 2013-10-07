import FWCore.ParameterSet.Config as cms

# CMSSW Regular Stuff
process = cms.Process("PU")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

#Options
process.GlobalTag.globaltag = 'START42_V13::All'

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('RMMEFN'))
process.source.fileNames = ['file:temp.root']

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

process.nPU     = cms.EDProducer("PileUpMultiplicityCounter",
    puLabel    = cms.InputTag("addPileupInfo"),
#     src        = cms.InputTag("REPLACE_ME"),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('nPU.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_*_*_PU',
    ),
)

process.p = cms.Path(process.nPU)
process.e = cms.EndPath(process.out)

