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

process.treeDumper.dumpPFCandidates = True
process.treeDumper.PFPUCandidateCollection = "pfPileUp"
if(isMC):
    process.treeDumper.dumpGenMet = True
    process.treeDumper.dumpGenJets = True
    process.treeDumper.dumpGenInfo = True
    process.treeDumper.dumpMCTruth = True
    process.treeDumper.dumpLHE = False
    process.treeDumper.dumpPdfWeight = False
if (runOnAOD == 1) :
    process.treeDumper.saveFatTrk = False
    process.treeDumper.saveTrackDeDx = False
    process.treeDumper.dumpPFlowElectrons = False
    process.treeDumper.dumpHcalNoiseFlags = True
    process.treeDumper.AODHcalNoiseFlags = True
else :
    process.treeDumper.saveFatTrk = True
    process.treeDumper.saveTrackDeDx = True
    process.treeDumper.dumpPFlowElectrons = True
    process.treeDumper.dumpHcalNoiseFlags = True
    process.treeDumper.AODHcalNoiseFlags = False
process.treeDumper.dumpTree = True


process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

if(doDiLeptonSkim==True):
    print "===> 2-LEPTON SKIM IS TUNRED ON! <==="
    process.skim = cms.Sequence( process.vecbosPreFilterDiLepSkimSeq )
else:
    process.skim = cms.Sequence( )

process.prejets = cms.Sequence( process.eCalibSequence
                                * process.leptonLinkedTracks 
                                * process.electronAndPhotonSuperClustersSequence
                                * process.chargedMetProducer
                                * process.pfIsolationAllSequence )
if(isMC==False):
    process.prejets *= process.metOptionalFilterSequence

if(process.treeDumper.dumpBCs==True):
    process.prejets *= process.seedBasicClustersSequence 

if(isMC):
    process.jets = cms.Sequence( process.ourJetSequenceMCReduced )
else:
    process.jets = cms.Sequence( process.ourJetSequenceDataReduced )

process.jets *= ( process.newBtaggingSequence 
                  * process.newPFPUcorrJetBtaggingSequence
                  * process.newPFNoPUJetBtaggingSequence
                  * process.metSequence )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedElectrons = cms.PSet(initialSeed = cms.untracked.uint32(1),
                                                                                  engineName = cms.untracked.string('TRandom3')
                                                                                  )
                                                   )

process.load("HiggsAnalysis.HiggsToWW2e.collectionsFilterSequence_cff")

process.postjets = cms.Sequence( process.eIdSequence
                                 * process.FastjetForIsolation
                                 * process.pfPileUp )

if(isMC):
    process.postjets *= ( process.collectionsFilterSequenceMC
                          * process.treeDumper )
else:
    process.postjets *= ( process.collectionsFilterSequenceData
                          * process.logErrorAnalysis
                          * process.treeDumper )

# In order to use the good primary vertices everywhere (It would be nicer to set the proper inputTags in the first place)
from PhysicsTools.PatAlgos.tools.helpers import *
massSearchReplaceAnyInputTag(process.prejets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.jets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.postjets,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)

if(process.treeDumper.dumpPdfWeight == False) :
    process.p = cms.Path ( process.skim
                           * process.goodPrimaryVertices
                           * process.prejets
                           * process.jets
                           * process.postjets
                           )
else :
    process.p = cms.Path ( process.skim
                           * process.pdfWeights
                           * process.goodPrimaryVertices
                           * process.prejets
                           * process.jets
                           * process.postjets
                           )

process.schedule = cms.Schedule(process.p)

