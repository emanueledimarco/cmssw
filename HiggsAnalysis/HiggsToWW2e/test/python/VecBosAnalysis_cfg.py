import FWCore.ParameterSet.Config as cms

# configure here #
isMC = True
is42X = False
is52X = False
runOnAOD = 1
doDiLeptonSkim = True
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# --- for PF candidates isolation
process.load("CommonTools.ParticleFlow.pfPileUp_cfi")

# --- skim sequences ---
process.load("HiggsAnalysis.Tools.multiLepSkimFilter_cff")

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFPUcorrJetsProducerSequence_cff")
process.load("WWAnalysis.Tools.chargedMetProducer_cfi")
process.chargedMetProducer.collectionTag = "particleFlow"
process.chargedMetProducer.vertexTag = "offlinePrimaryVertices"

# --- noise filters ---
process.load("HiggsAnalysis.HiggsToWW2e.METOptionalFilterFlags_cff")
if(is42X):
    process.hcalLaserEventFilter.vetoByRunEventNumber = True
    process.hcalLaserEventFilter.vetoByHBHEOccupancy= False

# --- tracker failures ---
process.load("MyAnalysis.METFlags.logErrorAnalysisProducer_cff")

# --- jet corrections ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")

# process.IgProfService = cms.Service("IgProfService",
#                                     reportFirstEvent            = cms.untracked.int32(0),
#                                     reportEventInterval         = cms.untracked.int32(50),
#                                     reportToFileAtPostEvent     = cms.untracked.string("| gzip -c > YYYY.%I.gz")
#                                     )

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonLinkedTracks_cfi")

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- electron regression and corrections ---
process.load("HiggsAnalysis.HiggsToWW2e.calibratedElectronsSequence_cff")
if(isMC):
    process.calibratedElectrons.isMC = cms.bool(True)                       
    process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")
else:
    process.calibratedElectrons.isMC = cms.bool(False)                       
    process.calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")

# --- pf isolation sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonPFIsoSequence_cff")

# --- ECAL clusters merging in a unique collection and seeds filtering ---
process.load("HiggsAnalysis.HiggsToWW2e.electronAndPhotonSuperClusters_cff")
process.load("HiggsAnalysis.HiggsToWW2e.seedBasicClusters_cff")

# --- good vertex filter ---
process.load("HiggsAnalysis.HiggsToWW2e.vertexFiltering_cff")

# --- PDF systematics ---
# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
                                    # Fix POWHEG if buggy (this PDF set will also appear on output,
                                    # so only two more PDF sets can be added in PdfSetNames if not "")
                                    #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
                                    GenTag = cms.untracked.InputTag("genParticles"),
                                    PdfInfoTag = cms.untracked.InputTag("generator"),
                                    PdfSetNames = cms.untracked.vstring("cteq66.LHgrid", "MRST2006nnlo.LHgrid", "NNPDF10_100.LHgrid")
                                    )

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
if(isMC):
    process.treeDumper.nameFile = 'default_MC.root'
else:
    process.treeDumper.nameFile = 'default_data.root'
if(isMC):
    process.treeDumper.PFJetCorrectionService = 'ak5PFL1FastL2L3'
    process.treeDumper.JetCorrectionService = 'ak5CaloL1FastL2L3'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpHLTObjects = True
process.treeDumper.dumpElectrons = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = True
process.treeDumper.dumpConversions = False
process.treeDumper.dumpVertices = True
process.treeDumper.dumpParticleFlowObjects = True
# for indirect lepton veto
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

process.jets *= ( process.newPFPUcorrJetBtaggingSequence
                  * process.metSequence )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedElectrons = cms.PSet(initialSeed = cms.untracked.uint32(1),
                                                                                  engineName = cms.untracked.string('TRandom3')
                                                                                  )
                                                   )

process.load("HiggsAnalysis.HiggsToWW2e.collectionsFilterSequence_cff")

process.postjets = cms.Sequence( process.eIdSequence
                                 * process.pfPileUp )

# this is needed only for 42X, since after that the isolation is corrected with the kt6PFJets rho
if(is42X):
    process.postjet *=  process.FastjetForIsolation


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

process.MessageLogger.destinations = ['cout', 'cerr']
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

if(process.treeDumper.dumpPdfWeight == False) :
    process.vecbosPath = cms.Path ( process.skim
                                    * process.goodPrimaryVertices
                                    * process.prejets
                                    * process.jets
                                    * process.postjets
                                    )
else:
    process.vecbosPath = cms.Path ( process.skim
                                    * process.pdfWeights
                                    * process.goodPrimaryVertices
                                    * process.prejets
                                    * process.jets
                                    * process.postjets
                                    )

process.schedule = cms.Schedule(process.vecbosPath)

