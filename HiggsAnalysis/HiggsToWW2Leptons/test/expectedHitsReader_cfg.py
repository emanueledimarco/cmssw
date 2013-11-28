import FWCore.ParameterSet.Config as cms

process = cms.Process('READER')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:./expectedHitsExample.root')
)

process.expectedHitsMapReader = cms.EDAnalyzer('ExpectedHitsMapReader'
)



process.run_step = cms.Path(process.expectedHitsMapReader)


