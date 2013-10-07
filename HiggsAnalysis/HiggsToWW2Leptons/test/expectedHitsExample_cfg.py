import FWCore.ParameterSet.Config as cms

process = cms.Process('EXAMPLE')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/m/mangano/HWW/mc/35X/H160_2W_2lnu_gluonfusion_7TeV.root')
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands =  cms.untracked.vstring(
        'keep *_*_*_RECO',
        'keep *_*_*_EXAMPLE',
        ),
    fileName = cms.untracked.string('expectedHitsExample.root'),
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START38_V12::All'

# Path and EndPath definitions
import HiggsAnalysis.HiggsToWW2Leptons.expectedHitsComputer_cfi
process.expectedHitsEle = HiggsAnalysis.HiggsToWW2Leptons.expectedHitsComputer_cfi.expectedHitsComputer.clone()
process.expectedHitsEle.inputColl   = cms.InputTag("gsfElectrons")
process.expectedHitsEle.useGsfTrack = cms.bool(True)
process.expectedHitsMu  = HiggsAnalysis.HiggsToWW2Leptons.expectedHitsComputer_cfi.expectedHitsComputer.clone()
process.expectedHitsMu.inputColl   = cms.InputTag("muons")
process.expectedHitsMu.useGsfTrack = cms.bool(False)


process.run_step = cms.Path(process.expectedHitsEle * process.expectedHitsMu)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.run_step,process.out_step)


