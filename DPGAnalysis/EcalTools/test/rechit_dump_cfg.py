import FWCore.ParameterSet.Config as cms

# configure here #
isMC = True
##################

process = cms.Process("ECALRecHitAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

if isMC:
    process.GlobalTag.globaltag = 'START17_61_V5::All'
    inputfile =  cms.untracked.vstring('file:testEcalLocalRecoA.root')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = inputfile
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load("DPGAnalysis.EcalTools.rechit_dumper_cfi")

process.mypath = cms.Path(process.rechitDumper)
