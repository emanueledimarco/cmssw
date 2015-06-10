import FWCore.ParameterSet.Config as cms


process = cms.Process("test")

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('ZEE central skim'),
    name = cms.untracked.string('$Source:  $')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_7_5_0_pre5/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_75_V4-v1/00000/00C7B7E5-7F0C-E511-9BF7-0025905A608C.root'
    ),
)

process.load("Configuration.Skimming.PDWG_ZeeSkim_cff")
process.zdiElectronFilter = cms.Path(process.zdiElectronSequence)

process.outputCsZee = cms.OutputModule("PoolOutputModule",
                                       dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('CS_ZEE')),
                                        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('zdiElectronFilter')),
                                        outputCommands = process.FEVTEventContent.outputCommands,
                                        fileName = cms.untracked.string('CS_ZEE.root')
                                        )


process.this_is_the_end = cms.EndPath(
process.outputCsZee
)
