import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0022/84618C45-33C0-DF11-9C9F-00304867C29C.root',
	'/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0022/4A02C098-E1BF-DF11-9CFA-00304867C136.root',
	'/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0022/305BA9E4-F9BF-DF11-A7C5-0018F3D09652.root',
	'/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0021/EC894065-98BF-DF11-9A15-0018F3D096F0.root',
	'/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0021/7CF81D64-9FBF-DF11-AA43-0018F3D09710.root',
	'/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0021/44771660-96BF-DF11-A1FC-00304867C0C4.root',
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('START38_V9::All')

from WWAnalysis.Filters.electronSelections_cff import *

#process.idElectrons = cms.EDFilter("GsfElectronSelector",
#    src = cms.InputTag("gsfElectrons"),
#    cut = cms.string(ELE_ID_CUT + " && " +ELE_ISO_CUT)
#)
from WWAnalysis.Filters.convRejectionSelection_cfi import convRejectionSelection
#process.goodElectrons = convRejectionSelection.clone(src = "idElectrons", filter = False)
process.convRejectionElectrons = convRejectionSelection.clone(src = "gsfElectrons", filter = False)

process.expectedHitsComputer = cms.EDProducer("ExpectedHitsComputer",
    #inputColl       = cms.InputTag("goodElectrons"),
    inputColl       = cms.InputTag("gsfElectrons"),
    useGsfTrack     = cms.bool(True),
    objectSelection = cms.string(""),
    propagator         = cms.string('PropagatorWithMaterialOpposite'),
    navigationSchool   = cms.string('SimpleNavigationSchool'),
    measurementTracker = cms.string(''),
)

process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi")
process.patElectrons.addElectronID = False
process.patElectrons.addGenMatch   = False
process.patElectrons.userData.userInts.src = cms.VInputTag(cms.InputTag("expectedHitsComputer"))
#process.patElectrons.electronSource = "goodElectrons"

process.patElectronsWithConvR = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("patElectrons"), 
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        convR = cms.PSet(
           src       = cms.InputTag("convRejectionElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""), # no preselection 
           deltaR              = cms.double(0.0001), # they're the same!
           checkRecoComponents = cms.bool(False), # doesn't work
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False), # overlaps don't cause the electron to be discared
        )
    ),
    finalCut = cms.string(''),
)

## Trigger matching (here comes the evil)
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.eleTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patElectronsWithConvR" ),
    matched = cms.InputTag( "patTrigger" ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltL1IsoRecoEcalCandidate', 'hltL1NonIsoRecoEcalCandidate' ),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)
process.eleIdTriggerMatchHLT = process.eleTriggerMatchHLT.clone(collectionTags = cms.vstring('hltPixelMatchElectronsL1Iso', 'hltPixelMatchElectronsL1NonIso') )
process.patElectronsWithTrigger = cms.EDProducer( "PATTriggerMatchElectronEmbedder",
    src     = cms.InputTag(  "patElectronsWithConvR" ),
    matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT'), cms.InputTag('eleIdTriggerMatchHLT'))
)
IDTriggerPSet = cms.PSet(
    passID  = cms.string(ELE_ID_CUT),
    passIso = cms.string(ELE_ISO_CUT),
    passConvR = cms.string("hasOverlaps('convR')"),
    hasHLTSC = cms.string("!triggerObjectMatchesByCollection('hltL1IsoRecoEcalCandidate').empty || !triggerObjectMatchesByCollection('hltL1NonIsoRecoEcalCandidate').empty"),
    hasHLTID = cms.string("!triggerObjectMatchesByCollection('hltPixelMatchElectronsL1Iso').empty || !triggerObjectMatchesByCollection('hltPixelMatchElectronsL1NonIso').empty"),
)
HLTs = [
    "HLT_Ele10_SW_EleId_L1R",             # MC
    "HLT_Ele15_SW_L1R",                   # MC
    "HLT_Ele12_SW_EleIdIsol_L1R",         # 2E31
    "HLT_Ele17_SW_LooseEleId_L1R",        # 2E31
    "HLT_Ele17_SW_CaloEleId_L1R",         # 2E31
    "HLT_DoubleEle10_SW_L1R",             # 2E31, MC
    "HLT_Ele12_SW_TighterEleIdIsol_L1R",  # 6E31
    "HLT_Ele17_SW_TightEleId_L1R",        # 6E31
    "HLT_DoubleEle15_SW_L1R",             # 6E31
]
for X in HLTs:
    setattr(IDTriggerPSet, X, cms.string("!triggerObjectMatchesByPath('%s').empty || !triggerObjectMatchesByPath('%s_v1').empty" %(X,X)))
IDTriggerPSet.HLT_Any   = cms.string("||".join(["!triggerObjectMatchesByPath('%s').empty || !triggerObjectMatchesByPath('%s_v1').empty" % (X,X) for X in HLTs]))
IDTriggerPSet.HLT_Any1E = cms.string("||".join(["!triggerObjectMatchesByPath('%s').empty || !triggerObjectMatchesByPath('%s_v1').empty" % (X,X) for X in HLTs if X.find("Double") == -1]))

process.tagElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string("pt > 20 && %s && %s" % (ELE_ID_CUT, ELE_ISO_CUT)),
)

process.probeElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string("!(1.4442 < abs(eta) < 1.5666)"), 
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('40 < mass < 140'),
    decay = cms.string('tagElectrons@+ probeElectrons@-')
)

process.tagElectronsMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    src = cms.InputTag("tagElectrons"),
    matched = cms.InputTag("genParticles"),
    pdgId = cms.vint32(11),
    distMin = cms.double(0.3),
)
process.probeElectronsMCMatch = process.tagElectronsMCMatch.clone(src = "probeElectrons")

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
        fbrem = cms.string("fbrem"),
        nEihOld = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfLostHits"),
        nEihNew = cms.string("userInt('expectedHitsComputer')"),
    ),
    flags = cms.PSet(
        IDTriggerPSet,
        passEihOld = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0"),
        passEihNew = cms.string("userInt('expectedHitsComputer') == 0"),
    ),
    tagVariables =  cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
    ),
    tagFlags = cms.PSet(IDTriggerPSet),
    isMC = cms.bool(True),
    tagMatches   = cms.InputTag("tagElectronsMCMatch"),
    probeMatches = cms.InputTag("probeElectronsMCMatch"),
    motherPdgId = cms.vint32(22, 23),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff  = cms.bool(True),
    allProbes = cms.InputTag("probeElectrons"),
)
process.tnpSimpleSequence = cms.Sequence(
    process.tagElectrons   * process.tagElectronsMCMatch   +
    process.probeElectrons * process.probeElectronsMCMatch +
    process.tpPairs    +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.convRejectionElectrons *
    process.expectedHitsComputer *
    process.patElectrons *
    process.patElectronsWithConvR *
    (process.patTrigger * (process.eleTriggerMatchHLT + process.eleIdTriggerMatchHLT)) *
    process.patElectronsWithTrigger *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZee_MC.root"))
