
import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from glob import glob
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'file:/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Electron/selected_15_1_1bu.root',
        'file:/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Electron/selected_16_1_pPn.root',
        'file:/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Electron/selected_17_1_VDf.root',
        'file:/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Mu/selected_14_1_MWu.root',
        'file:/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Mu/selected_15_1_Ozq.root',
        'file:/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Mu/selected_16_1_KDC.root',
        'file:/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Mu/selected_17_1_35q.root',
    )
)
process.source.fileNames =  []
process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Mu/*.root')]
#process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/zwwjets/Run2010A_Mu/*.root')]
process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_Electron/*.root')]
#process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/zwwjets/Run2010B_EG/*.root')]

#process.source.fileNames = [ 'file:WZEG_CentralSkim.root' ]
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V13::All')

from WWAnalysis.Filters.electronSelections_cff import *
process.goldenElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.1 && "+ELE_ID_CUT+" && "+ELE_ISO_CUT+
                      " && abs(dB) < 0.02 && " + 
                      #"!(1.4442 < abs(eta) < 1.5666) && " + 
                      "hasOverlaps('convR') && userInt('expectedHitsComputer') == 0")),
)
process.reCleanerElectrons = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("goldenElectrons"), 
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        eleSafe = cms.PSet(
           src       = cms.InputTag("patJets"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/ele.energy > 0.3 && jet.chargedHadronMultiplicity >= 4"),
           requireNoOverlaps = cms.bool(False),
        ),
        eleAggro = cms.PSet(
           src       = cms.InputTag("patJets"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/ele.energy > 0.2 && jet.chargedHadronMultiplicity >= 3"),
           requireNoOverlaps = cms.bool(False),
        ),
        eleSafeNoPU = cms.PSet(
           src       = cms.InputTag("patJetsNoPU"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/ele.energy > 0.3 && jet.chargedHadronMultiplicity >= 4"),
           requireNoOverlaps = cms.bool(False),
        ),
        eleAggroNoPU = cms.PSet(
           src       = cms.InputTag("patJetsNoPU"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/ele.energy > 0.2 && jet.chargedHadronMultiplicity >= 3"),
           requireNoOverlaps = cms.bool(False),
        ),
        any = cms.PSet(
           src       = cms.InputTag("patJets"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.pt > 10"),
           requireNoOverlaps = cms.bool(False),
        ),
        anyNoPU = cms.PSet(
           src       = cms.InputTag("patJetsNoPU"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.pt > 10"),
           requireNoOverlaps = cms.bool(False),
        ),
    ),
    finalCut = cms.string(''),
)

from WWAnalysis.Filters.muonSelections_cff import *
process.goldenMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.1 && "+MUON_ID_CUT+" && "+MUON_ISO_CUT+" && abs(dB) < 0.02")),
)
process.reCleanerMuons = cms.EDProducer("PATMuonCleaner",
    src = cms.InputTag("goldenMuons"), 
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        muSafe = cms.PSet(
           src       = cms.InputTag("patJets"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/mu.energy > 0.3 && jet.chargedHadronMultiplicity >= 4"),
           requireNoOverlaps = cms.bool(False),
        ),
        muAggro = cms.PSet(
           src       = cms.InputTag("patJets"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/mu.energy > 0.2 && jet.chargedHadronMultiplicity >= 3"),
           requireNoOverlaps = cms.bool(False),
        ),
        muSafeNoPU = cms.PSet(
           src       = cms.InputTag("patJetsNoPU"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/mu.energy > 0.3 && jet.chargedHadronMultiplicity >= 4"),
           requireNoOverlaps = cms.bool(False),
        ),
        muAggroNoPU = cms.PSet(
           src       = cms.InputTag("patJetsNoPU"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.chargedHadronEnergy/mu.energy > 0.2 && jet.chargedHadronMultiplicity >= 3"),
           requireNoOverlaps = cms.bool(False),
        ),
        any = cms.PSet(
           src       = cms.InputTag("patJets"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.pt > 10"),
           requireNoOverlaps = cms.bool(False),
        ),
        anyNoPU = cms.PSet(
           src       = cms.InputTag("patJetsNoPU"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string("jet.pt > 10"),
           requireNoOverlaps = cms.bool(False),
        ),
    ),
    finalCut = cms.string(''),
)

process.load("PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi")
process.patJets.jetSource = "ak5PFJets"
process.patJets.embedCaloTowers     = False
process.patJets.embedPFCandidates   = False
process.patJets.addJetCorrFactors   = False
process.patJets.addBTagInfo         = False
process.patJets.addAssociatedTracks = False
process.patJets.addJetCharge        = False
process.patJets.addJetID            = False
process.patJets.addGenPartonMatch   = False
process.patJets.addGenJetMatch      = False
process.patJets.getJetMCFlavour     = False
process.patJetsNoPU = process.patJets.clone(jetSource = "ak5PFJetsNoPU")

process.patJetsLeptonVeto = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("patJets"), 
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("goldenElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
        ),
        mu = cms.PSet(
           src       = cms.InputTag("goldenMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),  
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
        ),
    ),
    finalCut = cms.string(''),
)
process.patJetsLeptonVetoNoPU = process.patJetsLeptonVeto.clone(src = "patJetsNoPU")

process.reco = cms.Sequence( 
    ( process.goldenElectrons + 
      process.goldenMuons     + 
      process.patJets         + 
      process.patJetsNoPU     ) *
    ( process.reCleanerElectrons    +
      process.reCleanerMuons        +
      process.patJetsLeptonVeto     +
      process.patJetsLeptonVetoNoPU )
)


HLT1Es = [
    "HLT_Ele10_SW_EleId_L1R",             # MC
    "HLT_Ele15_SW_L1R",                   # MC
    "HLT_Ele12_SW_EleIdIsol_L1R",         # 2E31
    "HLT_Ele17_SW_LooseEleId_L1R",        # 2E31
    "HLT_Ele17_SW_CaloEleId_L1R",         # 2E31
    "HLT_Ele12_SW_TighterEleIdIsol_L1R",  # 6E31
    "HLT_Ele17_SW_TightEleId_L1R",        # 6E31
]
HLT1Ms = [
    "HLT_Mu9",
    "HLT_IsoMu9",
    "HLT_Mu11",
    "HLT_Mu15",
]
HLT_Any1E = "||".join(["!triggerObjectMatchesByPath('%s').empty || !triggerObjectMatchesByPath('%s_v1').empty" % (X,X) for X in HLT1Es])
HLT_Any1M = "||".join(["!triggerObjectMatchesByPath('%s').empty || !triggerObjectMatchesByPath('%s_v1').empty" % (X,X) for X in HLT1Ms])

process.tagElectrons =  cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("goldenElectrons"),
    cut = cms.string(HLT_Any1E),
)
process.tagMuons =  cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("goldenMuons"),
    cut = cms.string(HLT_Any1M),
)

process.nVerticesElectrons = cms.EDProducer("VertexMultiplicityCounter",
    probes  = cms.InputTag("reCleanerElectrons"),
    objects = cms.InputTag("offlinePrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"), 
)
process.nVerticesMuons = process.nVerticesElectrons.clone(probes = "reCleanerMuons")
process.nJetsElectrons = cms.EDProducer("CandMultiplicityCounter",
    probes  = cms.InputTag("reCleanerElectrons"),
    objects = cms.InputTag("patJetsLeptonVeto"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 25"), 
)
process.nJetsElectronsNoPU = process.nJetsElectrons.clone(objects = "patJetsLeptonVetoNoPU")
process.nJetsMuons = process.nJetsElectrons.clone(probes = "reCleanerMuons")
process.nJetsMuonsNoPU = process.nJetsMuons.clone(objects = "patJetsLeptonVetoNoPU")
process.varsElectrons = cms.Sequence(process.nVerticesElectrons + process.nJetsElectrons + process.nJetsElectronsNoPU)
process.varsMuons = cms.Sequence(process.nVerticesMuons + process.nJetsMuons + process.nJetsMuonsNoPU)

process.zElEl = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('tagElectrons@+ reCleanerElectrons@-'),
    cut = cms.string('mass > 40'),
)
process.zMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('tagMuons@+ reCleanerMuons@-'),
    cut = cms.string('mass > 40'),
)


process.tpTreeElEl = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("zElEl"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
        nVtx     = cms.InputTag("nVerticesElectrons"),
        nJet     = cms.InputTag("nJetsElectrons"),
        nJetNoPU = cms.InputTag("nJetsElectronsNoPU"),
    ),
    flags = cms.PSet(
        passSafe      = cms.string("!hasOverlaps('eleSafe')"),
        passSafeNoPU  = cms.string("!hasOverlaps('eleSafeNoPU')"),
        passAggro     = cms.string("!hasOverlaps('eleAggro')"),
        passAggroNoPU = cms.string("!hasOverlaps('eleAggroNoPU')"),
        hasJet        = cms.string("hasOverlaps('any')"),
        hasJetNoPU    = cms.string("hasOverlaps('anyNoPU')"),
    ),
    tagVariables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
    ),
    tagFlags = cms.PSet(),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
process.tpTreeMuMu = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("zMuMu"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
        nVtx     = cms.InputTag("nVerticesMuons"),
        nJet     = cms.InputTag("nJetsMuons"),
        nJetNoPU = cms.InputTag("nJetsMuonsNoPU"),
    ),
    flags = cms.PSet(
        passSafe      = cms.string("!hasOverlaps('muSafe')"),
        passSafeNoPU  = cms.string("!hasOverlaps('muSafeNoPU')"),
        passAggro     = cms.string("!hasOverlaps('muAggro')"),
        passAggroNoPU = cms.string("!hasOverlaps('muAggroNoPU')"),
        hasJet        = cms.string("hasOverlaps('any')"),
        hasJetNoPU    = cms.string("hasOverlaps('anyNoPU')"),
    ),
    tagVariables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
    ),
    tagFlags = cms.PSet(),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

process.tpTreeElElOnePair = process.tpTreeElEl.clone(arbitration = "OnePair")
process.tpTreeMuMuOnePair = process.tpTreeMuMu.clone(arbitration = "OnePair")

process.fElEl = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("zElEl"), minNumber = cms.uint32(1))
process.fMuMu = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("zMuMu"), minNumber = cms.uint32(1))

process.tnpSimpleSequenceElEl = cms.Path(
    process.reco          +
    process.tagElectrons  * 
    process.varsElectrons *
    process.zElEl         *
    process.fElEl         +
    process.tpTreeElEl    +
    process.tpTreeElElOnePair
)
process.tnpSimpleSequenceMuMu = cms.Path(
    process.reco        +
    process.tagMuons    * 
    process.varsMuons   *
    process.zMuMu       *
    process.fMuMu       +
    process.tpTreeMuMu  +
    process.tpTreeMuMuOnePair
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJetCleaning.root"))

if True:
    #proc = "Zmumu"
    #proc = "Zee"
    #process.source.fileNames = [ "file:z_ww_jets_MC_"+proc+"_powheg.root" ]
    #process.TFileService.fileName = "tnpJetCleaning_MC_"+proc+"_powheg.root"
    process.source.fileNames = [ "file:z_ww_jets_MC_Zee_powheg.root", "file:z_ww_jets_MC_Zmm_powheg.root" ]
    process.TFileService.fileName = "tnpJetCleaning_MC_powheg.root"
    #process.maxEvents.input = 100
    process.goldenElectrons.cut = ("pt > 20 && abs(eta) < 2.1 && "+ELE_ID_CUT+" && "+ELE_ISO_CUT)
    #process.goldenElectrons.filter = cms.bool(True)
    process.tagElectrons.cut = ""
    #process.goldenMuons.filter = cms.bool(True)
                                   #" && abs(dB) < 0.02 && ") + 
                                   #"!(1.4442 < abs(eta) < 1.5666) && " + 
                                   # "hasOverlaps('convR') && gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0")
