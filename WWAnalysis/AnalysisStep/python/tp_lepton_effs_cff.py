import FWCore.ParameterSet.Config as cms


#### FIXME: ID cuts are not up to date!
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT_OLD as MUON_ID_2010_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_ISO as MUON_ISO_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_IP as MUON_IP_CUT
#
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ID as ELE_ID_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ID2 as ELE_ID2_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ISO as ELE_ISO_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_IP as ELE_IP_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_CONV as ELE_NOCONV


HLT1Es = [ 'HLT_Ele27_WP80' ]
HLT1Ms = [ 'HLT_IsoMu24_eta2p1', 'HLT_IsoMu30_eta2p1', 'HLT_Mu40_eta2p1', 'HLT_Mu50_eta2p1' ]
HLT_Any1E = "||".join(["!triggerObjectMatchesByPath('%s_v*').empty" % (X,) for X in HLT1Es])
HLT_Any1M = "||".join(["!triggerObjectMatchesByPath('%s_v*').empty" % (X,) for X in HLT1Ms])

tagElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.5 " +
                      # "&& ("+ HLT_Any1E   +")" +   # NO, because we also use TnP triggers for DoubleEl
                      "&& ("+ ELE_ID_CUT  +")" +
                      "&& ("+ ELE_IP_CUT      +")" +
                      "&& ("+ ELE_NOCONV  +")" +
                      "&& ("+ ELE_ISO_CUT +")")),

)

tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.1 "+
                      "&& (" + HLT_Any1M + ")" +
                      "&& ("+MUON_ID_2010_CUT + ")" +
                      "&& ("+MUON_IP_CUT + ")" +
                      "&& ("+MUON_ISO_CUT+ ")")),

)

tp_lepton_effs_reco = cms.Sequence( 
    tagMuons     +
    tagElectrons
)

nVerticesElectrons = cms.EDProducer("VertexMultiplicityCounter",
    probes  = cms.InputTag("boostedElectrons"),
    objects = cms.InputTag("goodPrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"), 
)
nVerticesMuons = nVerticesElectrons.clone(probes = "boostedMuons")
nJetsElectrons = cms.EDProducer("CandCleanedMultiplicityCounter", # this is a property of the *pair*
    pairs   = cms.InputTag("tpElEl"),
    objects = cms.InputTag("slimPatJetsTriggerMatch"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 30"), 
    minTagObjDR   = cms.double(0.3),
    minProbeObjDR = cms.double(0.3),
)
nJetsMuons     = nJetsElectrons.clone(pairs = "tpMuMu")
varsElectrons = cms.Sequence(nVerticesElectrons + nJetsElectrons)
varsMuons = cms.Sequence(nVerticesMuons + nJetsMuons)

tpElEl = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('tagElectrons@+ boostedElectrons@-'),
    cut = cms.string('70 < mass < 120 && daughter(1).pt > 10'),
)
tpMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('tagMuons@+ boostedMuons@-'),
    cut = cms.string('70 < mass < 120 && daughter(1).pt > 10'),
)

TPCommonVariables = cms.PSet(
    pt    = cms.string("pt"),
    eta   = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    phi    = cms.string("phi"),
)

EleTriggers = cms.PSet(HLT_Any1E = cms.string(HLT_Any1E)); EleTriggersMC = cms.PSet()
MuTriggers  = cms.PSet(HLT_Any1M = cms.string(HLT_Any1M)); MuTriggersMC = cms.PSet()

# double-el phys triggers
EleTriggers.HLT_Ele17_Ele8_Ele8Leg = cms.string("!triggerObjectMatchesByPath('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*').empty")
EleTriggers.HLT_Ele17_Ele8_Ele17Leg = cms.string("!triggerObjectMatchesByFilter('hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter').empty")
# tag-and-probe triggers (two examples)
EleTriggers.HLT_Ele17_Ele8_TnP_Ele8Leg  = cms.string("!triggerObjectMatchesByPath('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*').empty")
EleTriggers.HLT_Ele17_Ele8_TnP_Ele17Leg = cms.string("!triggerObjectMatchesByFilter('hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter').empty")
EleTriggers.HLT_Ele32_SC17_TnP_SC17Leg  = cms.string("!triggerObjectMatchesByPath('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*').empty")
EleTriggers.HLT_Ele32_SC17_TnP_Ele32Leg = cms.string("!triggerObjectMatchesByFilter('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter').empty")
## FIXME: the Ele17Leg and Ele32Leg seem to be not working, should check if the filter names are correct and if we have the matching for those collections!

# Now the di-mu triggers also contain a per-pair filter on dz which probably has to be studied separately
MuTriggers.HLT_Mu17_Mu8_Mu17Leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17').empty")
MuTriggers.HLT_Mu17_Mu8_Mu8Leg  = cms.string("!triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8').empty")
MuTriggers.HLT_Mu17_Mu8_DZCut   = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*').empty")
MuTriggers.HLT_Mu17_TkMu8_Mu17Leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17').empty")
MuTriggers.HLT_Mu17_TkMu8_Mu8Leg  = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty")
MuTriggers.HLT_Mu17_TkMu8_DZCut   = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*').empty")

# FIXME should be updated
MuIDFlags = cms.PSet(
    passID  = cms.string(MUON_ID_CUT),
    passID_2010  = cms.string(MUON_ID_2010_CUT),
    passIP  = cms.string(MUON_IP_CUT),
    passIso = cms.string(MUON_ISO_CUT),
    passGlb = cms.string("isGlobalMuon"),
    passTM  = cms.string("isTrackerMuon"),
)
# FIXME note that ID must include also the loose pre-id!
EleIDFlags = cms.PSet(
    passID      = cms.string(ELE_ID_CUT),
    passID2      = cms.string(ELE_ID2_CUT),
    passIso     = cms.string(ELE_ISO_CUT),
    passIP      = cms.string(ELE_IP_CUT),
    passConvR   = cms.string(ELE_NOCONV),
)

tpTreeElEl = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpElEl"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        TPCommonVariables,
        nVtx     = cms.InputTag("nVerticesElectrons"),
    ),
    flags = cms.PSet(
        EleIDFlags,
        EleTriggers
    ),
    tagVariables = cms.PSet(),
    tagFlags = cms.PSet(
        EleTriggers
    ),
    pairVariables = cms.PSet(
        nJet     = cms.InputTag("nJetsElectrons"),
    ),
    pairFlags = cms.PSet(),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
tpTreeMuMu = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpMuMu"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        TPCommonVariables,
        nVtx     = cms.InputTag("nVerticesMuons"),
    ),
    flags = cms.PSet(
        MuIDFlags,
        MuTriggers,
    ),
    tagVariables = cms.PSet(),
    tagFlags = cms.PSet(
        MuTriggers
    ),
    pairVariables = cms.PSet(
        nJet     = cms.InputTag("nJetsMuons"),
    ),
    pairFlags = cms.PSet(),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

tnpSimpleSequenceElEl = cms.Sequence(
    tp_lepton_effs_reco *
    tpElEl              *
    varsElectrons       *
    tpTreeElEl
)
tnpSimpleSequenceMuMu = cms.Sequence(
    tp_lepton_effs_reco *
    tpMuMu              *
    varsMuons           *
    tpTreeMuMu
)
