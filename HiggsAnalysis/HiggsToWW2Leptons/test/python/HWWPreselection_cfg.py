import FWCore.ParameterSet.Config as cms

process = cms.Process("HWWPreselection")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')

process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWPreselectionSequence_cff")

process.options = cms.untracked.PSet(
    fileMode =  cms.untracked.string('NOMERGE')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring('file:/server/01a/ceballos/samples/HWW160_31x_01.root')
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histos.root')
                                   )

process.out = cms.OutputModule("PoolOutputModule",
                               verbose = cms.untracked.bool(False),
                               fileName = cms.untracked.string('preselection.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*_*_HWWPreselection',
                                                                      'keep *_muons_*_*',
                                                                      'keep *_gsfElectrons_*_*'
                                                                      )
                               )

process.p = cms.Path ( process.higgsToWW2LeptonsPreselectionSequence )

process.o = cms.EndPath ( process.out )
