import FWCore.ParameterSet.Config as cms

process = cms.Process("ZWWJetsSkim")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_8_3/RelValTTbar/GEN-SIM-RECO/START38_V9-v1/0022/CA9763E0-EFBF-DF11-81C5-002618943845.root',
        #'/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/MC_38Y_V9-v1/0021/3065AA71-9EBF-DF11-A946-0018F3D09676.root',
        'file:WZEG_CentralSkim.root'
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.GlobalTag.globaltag = cms.string('START38_V9::All')
process.GlobalTag.globaltag = cms.string('GR_R_38X_V13::All')

########### Event
process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
    filter = cms.bool(True),
)

process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
process.eventFilter = cms.Sequence(process.goodVertexFilter + process.noScraping)

############ MUONS
from WWAnalysis.Filters.muonSelections_cff import *
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *

process.goodMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 10 && " + MUON_ID_CUT + " && " + MUON_ISO_CUT),# + " && abs(dB) < 0.02")
)
process.zMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goodMuons@+ patMuonsWithTrigger@-"),
    cut = cms.string('mass > 40 && daughter(1).pt > 10'),
)
process.muonSequence = cms.Sequence(process.patMuonsWithTriggerSequence + process.goodMuons + process.zMuMu)

########### ELECTRONS
from WWAnalysis.Filters.electronSelections_cff import *
from WWAnalysis.Filters.convRejectionSelection_cfi import convRejectionSelection
process.convRejectionElectrons = convRejectionSelection.clone(src = "gsfElectrons", filter = False)
process.expectedHitsComputer = cms.EDProducer("ExpectedHitsComputer",
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
process.patElectronsWithConvR = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("patElectrons"), 
    preselection = cms.string('pt > 10'),
    checkOverlaps = cms.PSet(
        convR = cms.PSet(
           src       = cms.InputTag("convRejectionElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("pt > 10"), 
           deltaR              = cms.double(0.0001), # they're the same!
           checkRecoComponents = cms.bool(False), # doesn't work
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False), # overlaps don't cause the electron to be discared
        )
    ),
    finalCut = cms.string(''),
)

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
process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(("pt > 10 && "+ELE_ID_CUT+" && "+ELE_ISO_CUT)),
                      #" && abs(dB) < 0.02 && " + 
                      #"!(1.4442 < abs(eta) < 1.5666) && " + 
                      #"!hasOverlaps('convR') && userInt('expectedHitsComputer') == 0")),
)
process.zElEl = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('goodElectrons@+ patElectronsWithTrigger@-'),
    cut = cms.string('mass > 40 && daughter(1).pt > 10'),
)

process.electronSequence = cms.Sequence( 
    process.convRejectionElectrons *
    process.expectedHitsComputer *
    process.patElectrons *
    process.patElectronsWithConvR *
    (process.eleTriggerMatchHLT + process.eleIdTriggerMatchHLT) *
    process.patElectronsWithTrigger *
    process.goodElectrons *
    process.zElEl
)

process.lep = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("goodMuons"), cms.InputTag("goodElectrons"))
)
process.diLep = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('lep lep'),
    cut = cms.string('mass > 10'),
    checkCharge = cms.bool(False)
)


process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")
process.pfPileUp.Enable = True # otherwise it does nothing at all
process.ak5PFJetsNoPU = process.ak5PFJets.clone(src = cms.InputTag("pfNoPileUp"))

process.jetSequence = cms.Sequence(
    process.pfNoPileUpSequence +
    process.ak5PFJetsNoPU
)

#VBTF
import WWAnalysis.SkimStep.simpleCutBasedElectronIDSpring10_cfi as newcb
process.eidVBTFRel95 = newcb.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFCom95 = newcb.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFRel90 = newcb.simpleCutBasedElectronID.clone( electronQuality = '90relIso' )
process.eidVBTFCom90 = newcb.simpleCutBasedElectronID.clone( electronQuality = '90cIso'   )
process.eidVBTFRel85 = newcb.simpleCutBasedElectronID.clone( electronQuality = '85relIso' )
process.eidVBTFCom85 = newcb.simpleCutBasedElectronID.clone( electronQuality = '85cIso'   )
process.eidVBTFRel80 = newcb.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFCom80 = newcb.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )
process.eidVBTFRel70 = newcb.simpleCutBasedElectronID.clone( electronQuality = '70relIso' )
process.eidVBTFCom70 = newcb.simpleCutBasedElectronID.clone( electronQuality = '70cIso'   )
#CiC (requires RecoEgamma/ElectronIdentification V00-03-13)
import RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentification_cfi as cic
process.eidCiCVeryLoose = cic.eidVeryLoose.clone()
process.eidCiCLoose = cic.eidLoose.clone()
process.eidCiCMedium = cic.eidMedium.clone()
process.eidCiCTight = cic.eidTight.clone()
process.eidCiCSuperTight = cic.eidSuperTight.clone()
process.eidCiCHyperTight1 = cic.eidHyperTight1.clone()
process.eidCiCHyperTight2 = cic.eidHyperTight2.clone()
process.eidCiCHyperTight3 = cic.eidHyperTight3.clone()
process.eidCiCHyperTight4 = cic.eidHyperTight4.clone()
#CiC for HWW  (requires RecoEgamma/ElectronIdentification V00-03-13)
import RecoEgamma.ElectronIdentification.cutsInCategoriesHWWElectronIdentificationV04_cfi as cicHww
process.eidCiCHwwVeryLoose = cicHww.eidVeryLoose.clone()
process.eidCiCHwwLoose = cicHww.eidLoose.clone()
process.eidCiCHwwMedium = cicHww.eidMedium.clone()
process.eidCiCHwwTight = cicHww.eidTight.clone()
process.eidCiCHwwSuperTight = cicHww.eidSuperTight.clone()
process.eidCiCHwwHyperTight1 = cicHww.eidHyperTight1.clone()
process.eidCiCHwwHyperTight2 = cicHww.eidHyperTight2.clone()
process.eidCiCHwwHyperTight3 = cicHww.eidHyperTight3.clone()
process.eidCiCHwwHyperTight4 = cicHww.eidHyperTight4.clone()
process.eidSequence = cms.Sequence(
    process.eidVBTFRel95 + 
    process.eidVBTFCom95 + 
    process.eidVBTFRel90 + 
    process.eidVBTFCom90 + 
    process.eidVBTFRel85 + 
    process.eidVBTFCom85 + 
    process.eidVBTFRel80 + 
    process.eidVBTFCom80 +
    process.eidVBTFRel70 + 
    process.eidVBTFCom70 +
    process.eidCiCVeryLoose +
    process.eidCiCLoose +
    process.eidCiCMedium +
    process.eidCiCTight +
    process.eidCiCSuperTight +
    process.eidCiCHyperTight1 +
    process.eidCiCHyperTight2 +
    process.eidCiCHyperTight3 +
    process.eidCiCHyperTight4 +
    process.eidCiCHwwVeryLoose +
    process.eidCiCHwwLoose +
    process.eidCiCHwwMedium +
    process.eidCiCHwwTight +
    process.eidCiCHwwSuperTight +
    process.eidCiCHwwHyperTight1 +
    process.eidCiCHwwHyperTight2 +
    process.eidCiCHwwHyperTight3 +
    process.eidCiCHwwHyperTight4
)
process.newExpectedHitsElectrons = cms.EDProducer("ExpectedHitsComputer",
    inputColl       = cms.InputTag("gsfElectrons"),
    useGsfTrack     = cms.bool(True),
    objectSelection = cms.string(""),
    propagator         = cms.string('PropagatorWithMaterialOpposite'),
    navigationSchool   = cms.string('SimpleNavigationSchool'),
    measurementTracker = cms.string(''),
)
process.newExpectedHitsMuons = cms.EDProducer("ExpectedHitsComputer",
    inputColl       = cms.InputTag("muons"),
    useGsfTrack     = cms.bool(False),
    objectSelection = cms.string(""),
    propagator         = cms.string('PropagatorWithMaterialOpposite'),
    navigationSchool   = cms.string('SimpleNavigationSchool'),
    measurementTracker = cms.string(''),
)
process.newExpectedHits = cms.Sequence(process.newExpectedHitsElectrons + process.newExpectedHitsMuons)

process.reco = cms.Sequence(
    process.muonSequence + process.electronSequence + 
    process.lep + process.diLep + 
    process.jetSequence + 
    process.eidSequence +
    process.newExpectedHits
)

process.fMM = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("zMuMu"), minNumber = cms.uint32(1))
process.fEE = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("zElEl"), minNumber = cms.uint32(1))
process.fLL = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("diLep"), minNumber = cms.uint32(1))

process.pMM = cms.Path(process.eventFilter + process.reco + process.fMM)
process.pEE = cms.Path(process.eventFilter + process.reco + process.fEE)
process.pLL = cms.Path(process.eventFilter + process.reco + process.fLL)

if False:
    process.genLeptonFiducial = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("(abs(pdgId) == 11 || abs(pdgId) == 13) && status == 3 && abs(eta) < 2.4 && pt > 20"),
    )
    process.diGenLeptonFiducial = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string("genLeptonFiducial genLeptonFiducial"),
        cut = cms.string("mass > 10"),
        checkCharge = cms.bool(False),
    )
    process.diGenLeptonFiducialFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("diGenLeptonFiducial"), minNumber = cms.uint32(1))
    process.pGG  = cms.Path(process.eventFilter + process.reco + process.genLeptonFiducial + process.diGenLeptonFiducial + process.diGenLeptonFiducialFilter)
    process.pGGX = cms.Path(process.genLeptonFiducial + process.diGenLeptonFiducial + process.diGenLeptonFiducialFilter)

if False:
    process.zLooseElEl = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string("gsfElectrons gsfElectrons"),
        cut = cms.string("mass > 40 && min(daughter(0).pt, daughter(1).pt) > 10"),
        checkCharge = cms.bool(False),
    )
    process.zLooseMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string("muons muons"),
        cut = cms.string("mass > 40 && min(daughter(0).pt, daughter(1).pt) > 10"),
        checkCharge = cms.bool(False),
    )
    process.flEE = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("zLooseElEl"), minNumber = cms.uint32(1))
    process.flMM = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("zLooseMuMu"), minNumber = cms.uint32(1))
    process.lEE = cms.Path(process.zLooseElEl + process.flEE)
    process.lMM = cms.Path(process.zLooseMuMu + process.flMM)


process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("selected.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("pEE","pMM","pLL")),
    outputCommands = cms.untracked.vstring(
        "drop *",
        # This stuff
        "keep *_patMuonsWithTrigger_*_*",
        "keep *_patElectronsWithTrigger_*_*",
        "keep *_goodMuons_*_*",
        "keep *_goodElectrons_*_*",
        "keep *_zMuMu_*_*",
        "keep *_zElEl_*_*",
        "keep *_diLep_*_*",
        "keep *_ak5PFJets_*_*",
        "keep *_ak5PFJetsNoPU_*_*",
        # Tracking
        "keep *_offlinePrimaryVertices_*_*",
        "keep *_offlinePrimaryVerticesWithBS_*_*",
        "keep *_offlineBeamSpot_*_*",
        "keep recoTracks_generalTracks_*_*",
        "keep *_newExpectedHits*_*_*",
        # Muon
        "keep *_muons_*_*",
        "keep recoTracks_globalMuons_*_*",
        "keep recoTracks_standAloneMuons_*_*",
        "keep recoTrackExtras_standAloneMuons_*_*",
        # E/Gamma
        "keep *_photons_*_*",
        "keep *_gsfElectrons_*_*",
        "keep *_gsfElectronCores_*_*",
        "keep *_pfElectronTranslator_*_*",
        "keep recoSuperClusters_*_*_*",
        "keep recoGsfTracks_electronGsfTracks_*_*",
        "keep *_eidVBTF*_*_*",
        "keep *_eidCiC*_*_*",
        # PF
        "keep recoPFCandidates_particleFlow__*",
        # tcMET
        "keep *_pfMet_*_*",
        "keep *_tcMet_*_*",
        # MC
        "keep *_genParticles_*_*",
        "keep *_genMetTrue_*_",
        # Trigger
        "keep *_hltTriggerSummaryAOD_*_*",
        "keep *_TriggerResults_*_*",
    )

)
process.end = cms.EndPath(process.out)
