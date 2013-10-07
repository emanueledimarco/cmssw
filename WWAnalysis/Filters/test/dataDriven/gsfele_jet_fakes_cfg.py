import FWCore.ParameterSet.Config as cms

process = cms.Process("EleFakes")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

# import of standard configurations
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START38_V9::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_8_3/RelValZEE/GEN-SIM-RECO/MC_38Y_V9-v1/0021/3065AA71-9EBF-DF11-A946-0018F3D09676.root',
        '/store/relval/CMSSW_3_8_3/RelValTTbar/GEN-SIM-RECO/START38_V9-v1/0022/CA9763E0-EFBF-DF11-81C5-002618943845.root',
        '/store/relval/CMSSW_3_8_3/RelValWE/GEN-SIM-RECO/START38_V9-v1/0022/B625617F-F2BF-DF11-BBA3-001A928116B8.root',
        '/store/relval/CMSSW_3_8_3/RelValWjet_Pt_80_120/GEN-SIM-RECO/MC_38Y_V9-v1/0021/F4D78269-9CBF-DF11-AF7B-0026189438FE.root',
        '/store/relval/CMSSW_3_8_3/RelValZTT/GEN-SIM-RECO/START38_V9-v1/0021/449CD467-97BF-DF11-A771-0030486790BA.root',
        '/store/relval/CMSSW_3_8_3/RelValPhotonJets_Pt_10/GEN-SIM-RECO/START38_V9-v1/0022/96D05423-EEBF-DF11-8FED-002618943946.root',
        '/store/relval/CMSSW_3_8_3/RelValQCD_Pt_80_120/GEN-SIM-RECO/START38_V9-v1/0021/949009E5-9BBF-DF11-BD41-002618FDA207.root',
	'/store/relval/CMSSW_3_8_3/RelValQCD_FlatPt_15_3000/GEN-SIM-RECO/MC_38Y_V9-v1/0022/BE356752-A2BF-DF11-9F11-00261894388A.root',
    ),
)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("ele_fakes.root"),
    outputCommands = cms.untracked.vstring('drop *'),
)

##    _____ _           _                       
##   | ____| | ___  ___| |_ _ __ ___  _ __  ___ 
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \/ __|
##   | |___| |  __/ (__| |_| | | (_) | | | \__ \
##   |_____|_|\___|\___|\__|_|  \___/|_| |_|___/
##                                              
##   
from WWAnalysis.Filters.electronSelections_cff import *
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
process.makePatElectrons = cms.Sequence(
    process.convRejectionElectrons *
    process.expectedHitsComputer *
    process.patElectrons *
    process.patElectronsWithConvR
)

##        _      _       
##       | | ___| |_ ___ 
##    _  | |/ _ \ __/ __|
##   | |_| |  __/ |_\__ \
##    \___/ \___|\__|___/
##                       
##   

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

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.patTrigger.onlyStandAlone = True

process.jetTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patJets" ),
    matched = cms.InputTag( "patTrigger" ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltMCJetCorJetIcone5HF07' ),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)
process.patJetsWithTrigger = cms.EDProducer( "PATTriggerMatchJetEmbedder",
    src     = cms.InputTag(  "patJets" ),
    matches = cms.VInputTag("jetTriggerMatchHLT")
)
process.patJetsPassingTrigger = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJetsWithTrigger"),
    cut = cms.string("!triggerObjectMatchesByPath('HLT_Jet15U' ).empty() ||"+
                     "!triggerObjectMatchesByPath('HLT_Jet30U' ).empty() ||"+
                     "!triggerObjectMatchesByPath('HLT_Jet50U' ).empty() ||"+
                     "!triggerObjectMatchesByPath('HLT_Jet70U' ).empty() ||"+
                     "!triggerObjectMatchesByPath('HLT_Jet100U').empty()"),
)

process.patJetsWithMatches = cms.EDProducer("PATJetCleaner",
    ## pat electron input source
    src = cms.InputTag("patJetsWithTrigger"), 
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string('pt > 10'),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        gsfEle = cms.PSet(
           src       = cms.InputTag("patElectronsWithConvR"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("pt > 10"),  
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False),
        ),
    ),
    finalCut = cms.string(''),
)

process.makePatJets = cms.Sequence(
    process.patTrigger +
    process.patJets * 
    process.jetTriggerMatchHLT * 
    process.patJetsWithTrigger *
    process.patJetsPassingTrigger  +
    process.patJetsWithMatches
)
##    _____                 _     _____ _ _ _                
##   | ____|_   _____ _ __ | |_  |  ___(_) | |_ ___ _ __ ___ 
##   |  _| \ \ / / _ \ '_ \| __| | |_  | | | __/ _ \ '__/ __|
##   | |___ \ V /  __/ | | | |_  |  _| | | | ||  __/ |  \__ \
##   |_____| \_/ \___|_| |_|\__| |_|   |_|_|\__\___|_|  |___/
##                                                           
##   
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
process.goodEl = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectronsWithConvR"),
    cut = cms.string("pt > 10 && "+ELE_ID_CUT+" && "+ELE_ISO_CUT),
)
process.diEl = cms.EDProducer("CandViewAlwaysShallowCloneCombiner",
    cut = cms.string('daughter(1).pt > 10 && mass > 2'),
    decay = cms.string('goodEl patElectronsWithConvR'),
    checkCharge = cms.bool(False),
)
process.eventFilters = cms.Sequence(
    process.goodVertexFilter + 
    process.noScraping +
    process.goodEl *  process.diEl
)


##    _____              ___     ____            _          
##   |_   _|_ _  __ _   ( _ )   |  _ \ _ __ ___ | |__   ___ 
##     | |/ _` |/ _` |  / _ \/\ | |_) | '__/ _ \| '_ \ / _ \
##     | | (_| | (_| | | (_>  < |  __/| | | (_) | |_) |  __/
##     |_|\__,_|\__, |  \___/\/ |_|   |_|  \___/|_.__/ \___|
##              |___/                                       
##   
process.tpPairsJJ = cms.EDProducer("CandViewAlwaysShallowCloneCombiner",
    cut = cms.string('deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.5'), # avoid self-coupling
    decay = cms.string('patJetsPassingTrigger patJetsWithMatches'),
    checkCharge = cms.bool(False),
)
process.tpPairsJE = cms.EDProducer("CandViewAlwaysShallowCloneCombiner",
    cut = cms.string('deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.5'), # avoid self-coupling
    decay = cms.string('patJetsPassingTrigger patElectronsWithConvR'),
    checkCharge = cms.bool(False),
)

# additional variables to possibly fit away true signal
process.metVariable = cms.EDProducer("OtherCandVariableComputer",
    probes  = cms.InputTag("patJetsPassingTrigger"), # actually tags, but it doesn't matter
    objects = cms.InputTag("pfMet"),
    expression = cms.string("pt"),
    default    = cms.double(0), 
)
process.meeMass = cms.EDProducer("OtherCandVariableComputer",
    probes  = cms.InputTag("patJetsPassingTrigger"), # actually tags, but it doesn't matter
    objects = cms.InputTag("diEl"),
    objectsSortDescendingBy = cms.string("mass"), # pick the one with highest mass
    expression = cms.string("mass"),
    default    = cms.double(-1), 
)
process.meeCharge = cms.EDProducer("OtherCandVariableComputer",
    probes  = cms.InputTag("patJetsPassingTrigger"), # actually tags, but it doesn't matter
    objects = cms.InputTag("diEl"),
    objectsSortDescendingBy = cms.string("mass"), # pick the one with highest mass
    expression = cms.string("charge"),
    default    = cms.double(0), 
)

process.tpTreeJJ = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairsJJ"),
    arbitration   = cms.string("NonDuplicate"),
    # probe variables: all useful ones
    variables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
        n60   = cms.string("nCarrying(0.6)"),
        n90   = cms.string("nCarrying(0.9)"),
    ),
    flags = cms.PSet(gsf = cms.string("hasOverlaps('gsfEle')")),
    tagVariables = cms.PSet(
        pfmet = cms.InputTag("metVariable"),
        mee   = cms.InputTag("meeMass"),
        qee   = cms.InputTag("meeCharge"),
    ),
    tagFlags     = cms.PSet(),
    isMC         = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
for X in (15, 30, 50, 70, 100):
    setattr(process.tpTreeJJ.tagFlags, "HLT_Jet%dU"%X, cms.string("!triggerObjectMatchesByPath('HLT_Jet%dU' ).empty()"%X))

process.tpTreeJE = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairsJE"),
    arbitration   = cms.string("NonDuplicate"),
    # probe variables: all useful ones
    variables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
    ),
    flags = cms.PSet(
        passID  = cms.string(ELE_ID_CUT),
        passIso = cms.string(ELE_ISO_CUT),
        passConvR = cms.string("hasOverlaps('convR')"),
        passEihOld = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0"),
        passEihNew = cms.string("userInt('expectedHitsComputer') == 0"),
    ),
    tagFlags     = process.tpTreeJJ.tagFlags.clone(),
    tagVariables = process.tpTreeJJ.tagVariables.clone(),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

if False:
    process.load("PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi")
    process.electronMatch.mcStatus    = [3]
    process.electronMatch.checkCharge = False
    process.patElectrons.addGenMatch  = True
    process.makePatElectrons.replace(process.patElectrons, process.electronMatch + process.patElectrons)
    process.tpTreeJE.flags.trueElectron = cms.string("genParticleRef.isNonnull")
    process.fiducialGenElectron = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("abs(pdgId) == 11 && status == 3 && abs(eta) < 2.5 && pt > 10"),
        filter = cms.bool(True),
    )
    process.fiducialRecElectron = cms.EDFilter("PATElectronSelector",
        src = cms.InputTag("patElectronsWithConvR"),
        cut = cms.string("genParticleRef.isNonnull"),
        filter = cms.bool(True),
    )
    process.oneRecElectron = cms.EDFilter("PATElectronSelector",
        src = cms.InputTag("patElectronsWithConvR"),
        cut = cms.string("pt > 10 && "+ELE_ID_CUT+" && "+ELE_ISO_CUT),
        filter = cms.bool(True),
    )
    process.testFiltersGen_D = cms.Path(process.makePatElectrons * process.oneRecElectron + process.fiducialGenElectron)
    process.testFiltersGen_N = cms.Path(process.testFiltersGen_D._seq + process.eventFilters)
    process.testFiltersRec_D = cms.Path(process.makePatElectrons * process.oneRecElectron + process.fiducialRecElectron)
    process.testFiltersRec_N = cms.Path(process.testFiltersRec_D._seq + process.eventFilters)

process.fakeRate = cms.Path(
    process.makePatElectrons *
    process.makePatJets      *
    process.eventFilters *
    (process.metVariable + process.meeMass + process.meeCharge) +
    (process.tpPairsJJ * process.tpTreeJJ +
     process.tpPairsJE * process.tpTreeJE )
)


process.out.outputCommands = [ 
    "drop *", 
    "keep *_TriggerResults_*_*",
    "keep *_patTrigger_*_*", 
    "keep *_patJetsWithMatches_*_*",
    "keep *_patJetsPassingTrigger_*_*",
    "keep *_gsfElectrons_*_*",
    "keep *_patElectrons_*_*",
    "keep *_selectedPatElectronsPF_*_*",
    "keep *_tpPairs_*_*",
]
#process.end = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo_ele_fakes.root"))
