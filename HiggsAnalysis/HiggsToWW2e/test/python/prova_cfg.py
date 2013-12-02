import FWCore.ParameterSet.Config as cms

# configure here #
isMC = False
is42X = False
is52X = False
runOnAOD = 1
doDiLeptonSkim = False
##################

process = cms.Process("VecBosAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

if(isMC):
    if(is42X):
        process.GlobalTag.globaltag = 'START42_V17::All'
    elif(is52X):
        process.GlobalTag.globaltag = 'START52_V11C::All'
    else:
        process.GlobalTag.globaltag = 'START53_V27::All'
else:
    if(is42X):
        process.GlobalTag.globaltag = 'GR_R_42_V25::All'
    elif(is52X):
        process.GlobalTag.globaltag = 'GR_R_52_V9D::All'
    else:
        process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'

if(isMC):
    if(is42X):
        inputfile = cms.untracked.vstring('file:/cmsrm/pc24_2/emanuele/data/DYeeSummer11.root')
    elif(is52X):
        inputfile = cms.untracked.vstring('file:/cmsrm/pc21_2/emanuele/data/AOD_WWSummer12.root')
    else:
        inputfile = cms.untracked.vstring('/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/GEN-SIM-RECO/HLT8E33_PU_S10_START53_V7I-v1/30000/FCC5A42A-1A7A-E211-B82D-0025901D4B08.root')
else:
    if(is42X):
        inputfile = cms.untracked.vstring('file:/cmsrm/pc24_2/emanuele/data/reRecoMay10File.root')
    elif(is52X):
        inputfile = cms.untracked.vstring('/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/191/700/00327508-AF8B-E111-8151-BCAEC53296F4.root')
    else:
        inputfile = cms.untracked.vstring('/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10027/7A1610EB-0990-E211-BC9F-002354EF3BE6.root')


process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = inputfile
                            )

process.MessageLogger.destinations = ['cout', 'cerr']
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 10

# --- electron regression and corrections ---
process.load("HiggsAnalysis.HiggsToWW2e.calibratedElectronsSequence_cff")
if(isMC):
    process.calibratedElectrons.isMC = cms.bool(True)                       
    process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")
else:
    process.calibratedElectrons.isMC = cms.bool(False)                       
    process.calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")

process.p = cms.Path(process.eCalibSequence)

process.schedule = cms.Schedule(process.p)

