import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/53X/ggH125/hzz4lSkim_4_1_am6.root'
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V10::All'   #for 53X MC  

process.load("WWAnalysis.AnalysisStep.hzz4l_selection_cff")
from WWAnalysis.AnalysisStep.zz4l.hzz4l_selection_2012_fsr_cff import *  
isMC = True

import os
cmsswVer=os.environ["CMSSW_VERSION"]
releaseVer="53X" #default
if "CMSSW_5_3_" in cmsswVer:
    releaseVer="53X"
if "CMSSW_4_2_" in cmsswVer:
    releaseVer="42X"

process.load("WWAnalysis.AnalysisStep.regressionEnergyPatElectrons_cfi")
process.boostedRegressionElectrons = process.regressionEnergyPatElectrons.clone()
process.boostedRegressionElectrons.energyRegressionType = cms.uint32(1)
process.boostedRegressionElectrons.regressionInputFile = cms.string("EGamma/EGammaAnalysisTools/data/eleEnergyRegWeights_V1.root")
process.boostedRegressionElectrons.inputPatElectronsTag = "boostedElectronsID"
process.boostedRegressionElectrons.rhoCollection = cms.InputTag("kt6PFJets","rho","RECO")
process.boostedRegressionElectrons.debug = cms.bool(False)

process.load("WWAnalysis.AnalysisStep.calibratedPatElectrons_cfi")
process.calibratedPatElectrons.energyMeasurementType = cms.uint32(1)
if not hasattr(process, 'RandomNumberGeneratorService'):
    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
process.RandomNumberGeneratorService.boostedElectrons2 = cms.PSet(
    initialSeed = cms.untracked.uint32(1),
    engineName = cms.untracked.string('TRandom3')
)

process.boostedElectrons2 = process.calibratedPatElectrons.clone()
process.boostedElectrons2.isMC = isMC
if isMC : 
    if releaseVer == "42X" : process.boostedElectrons2.inputDataset = 'Fall11'
    else                   : process.boostedElectrons2.inputDataset = 'Summer12_HCP2012'
else    : 
    if releaseVer == "42X" : process.boostedElectrons2.inputDataset = 'Jan16ReReco'
    else                   : process.boostedElectrons2.inputDataset = 'HCP2012'
process.boostedElectrons2.updateEnergyError = cms.bool(True)
process.boostedElectrons2.isAOD = cms.bool(True)
process.postreboosting = cms.Sequence(
    process.boostedRegressionElectrons * process.boostedElectrons2 * boostedElectronsEAPFIso * boostedElectrons 
)
process.boostedRegressionElectrons.inputPatElectronsTag = "boostedElectronsID"
process.boostedElectrons2.inputPatElectronsTag = "boostedRegressionElectrons"
process.boostedElectronsEAPFIso.src = "boostedElectrons2"   

process.goodEl = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string("(%s) && (%s)" % (ELID_LOOSE,ELID_GOOD)),
)


process.load("PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi")
process.electronMatch.src     = 'goodEl'
process.electronMatch.matched = 'prunedGen'
process.electronMatch3 = process.electronMatch.clone(mcStatus = [3])
process.pull1 = cms.EDProducer("PullsComputer", probes = cms.InputTag("goodEl"), matches = cms.InputTag("electronMatch"))
process.pull3 = cms.EDProducer("PullsComputer", probes = cms.InputTag("goodEl"), matches = cms.InputTag("electronMatch3"))

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140'),
    decay = cms.string('goodEl@+ goodEl@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

VARS = cms.PSet(
    pt  = cms.string("pt"),
    p   = cms.string("p"),
    eta = cms.string("eta"),
    phi = cms.string("phi"),
    charge = cms.string("charge"),
    gen_p   = cms.InputTag("pull1","p"),
    gen_pt  = cms.InputTag("pull1","pt"),
    gen_eta = cms.InputTag("pull1","eta"),
    gen_phi = cms.InputTag("pull1","phi"),
    gen3_p   = cms.InputTag("pull3","p"),
    gen3_pt  = cms.InputTag("pull3","pt"),
    gen3_eta = cms.InputTag("pull3","eta"),
    gen3_phi = cms.InputTag("pull3","phi"),
    ptErr   = cms.string("? ecalDriven ? (pt/p) * p4Error('P4_COMBINATION') : gsfTrack.ptError"),
)
FLAGS = cms.PSet(
    PF = cms.string("1"),
)

process.superPair = cms.EDProducer("SkimEvent2LProducer",
    src = cms.InputTag("tpPairs"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    doMassRes = cms.bool(True),
    doExtendedMassRes = cms.bool(True),
    isMC = cms.bool(False),
)
process.bestZ = cms.EDFilter("SkimEvent2LSorter",
    src = cms.InputTag("superPair"),
    sortBy = cms.string("abs(mass - 91.188)"),
    sortAscending = cms.bool(True), 
    maxNumber = cms.uint32(1),
)

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("bestZ"),
    arbitration   = cms.string("OnePair"),
    # probe variables: all useful ones
    variables = cms.PSet(VARS),
    flags = cms.PSet(FLAGS),
    tagVariables = cms.PSet(VARS),
    tagFlags = cms.PSet(FLAGS),
    pairVariables = cms.PSet( 
        massErr = cms.string("userFloat('massErr')"), 
        massErr_tag   = cms.string("userFloat('massErr[0]')"),
        massErr_probe = cms.string("userFloat('massErr[1]')"),
    ),
    pairFlags = cms.PSet(),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

process.tnpSimpleSequence = cms.Sequence(
    process.boostedElectronsID +
    process.postreboosting + 
    process.goodEl +
    process.electronMatch3 + process.electronMatch + process.pull1 + process.pull3 + 
    process.tpPairs    +
    process.onePair    +
    process.superPair +
    process.bestZ +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("scaleZEE_MC.root"))
