import FWCore.ParameterSet.Config as cms

process = cms.Process("EX")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
##process.GlobalTag.globaltag = 'START42_V12::All'   # CMSSW_4XY
process.GlobalTag.globaltag = 'START52_V4::All'    # CMSSW_52Y

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'rfio:/castor/cern.ch/user/b/benedet/Fall11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_0.root'  #CMSSW_4XY
    'root://eoscms//eos/cms//store/relval/CMSSW_5_2_2/RelValProdTTbar/GEN-SIM-RECO/START52_V4-v1/0004/DCA3A2E4-C973-E111-9839-003048F1183E.root'   # CMSSW_52Y
    # 'file:/data/benedet/Fall11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_0.root'
    ),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )


# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )

#my analyzer
process.demo = cms.EDAnalyzer("MuonAnalyzer")

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("myhisto.root")
    )

process.pAna = cms.Path(process.demo)

process.schedule = cms.Schedule(process.pAna)





