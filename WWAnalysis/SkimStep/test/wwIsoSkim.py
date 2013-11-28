import FWCore.ParameterSet.Config as cms

process = cms.Process("WWIso")

isMC = True
# isMC = RMMEMC
isAOD = False
# isAOD = RMMEAOD

#  _____               _____                               _                
# |  __ \             |  __ \                             | |               
# | |__) |   _ _ __   | |__) |_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
# |  _  / | | | '_ \  |  ___/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# | | \ \ |_| | | | | | |  | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# |_|  \_\__,_|_| |_| |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
# 

from PhysicsTools.PatAlgos.tools.coreTools import *
 
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# process.GlobalTag.globaltag = 'START3X_V26::All'
# process.GlobalTag.globaltag = 'START38_V12::All'
# process.GlobalTag.globaltag = 'RMMEGlobalTag'
process.GlobalTag.globaltag = 'START39_V8::All'


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# process.MessageLogger.suppressError = cms.untracked.vstring('patTrigger')

# process.load("WWAnalysis.SkimStep.logger_cfi")
# process.MessageLogger._categoryCanTalk('L1GtUtils')
# process.MessageLogger._moduleCanTalk('patTrigger')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#         '/store/mc/Winter10/WJetsToLNu_TuneD6T_7TeV-madgraph-tauola/AODSIM/E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/0000/2E3D58AE-420F-E011-83B7-002618943918.root'
#         '/store/data/Run2010B/Electron/AOD/Dec22ReReco_v1/0014/F6B6366D-EB0E-E011-8DFD-0018F3D09658.root'
#        'RMMEFileName'
#         'file:/home/mwlebour/data/hww/ElDec22.147451.1.AOD.root',
#         'file:/home/mwlebour/data/hww/ElDec22.149291.79.AOD.root'
        'file:/home/mwlebour/data/WW.38XMC.Samples/DYToEEM20CT10Z2powheg.root'
#         'file:/home/mwlebour/data/Winter10/Hww160.root'
#         'file:/home/mwlebour/data/Winter10/WJetsMadgrap.root'
    )
)


#from glob import glob
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/oct29Skim/WWFull/*/*.root')]

process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string('hwwSkim.root'),
#     fileName = cms.untracked.string('RMMEFileName'),
#     fileName = cms.untracked.string('Hww160.skim.root'),
    fileName = cms.untracked.string('WJetsMad.skim.root'),
    outputCommands =  cms.untracked.vstring(
        'drop *',
        # This stuff
        'keep *_goodMuons_*_*',
        'keep *_goodElectrons_*_*',
        #'keep patJets_patJets_*_*',
        'keep patJets_cleanPatJetsPF_*_*',
        'keep patJets_cleanPatJetsNoPU_*_*',
        # Tracking
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_offlineBeamSpot_*_*',
        # MET
        'keep *_tcMet_*_*',
        'keep *_pfMet_*_*',
        # MC
        'keep *_prunedGen_*_*',
        'keep *_genMetTrue_*_*',
        'keep *_highPtGenJets_*_*',
        #'keep recoGenJets_ak5GenJets_*_*',
        # Trigger
        #'keep *_hltTriggerSummaryAOD_*_*',
        'keep *_TriggerResults_*_*',
        'keep recoTracks_generalTracks_*_*',
#         'keep *_ecalRecHit_EcalRecHitsEB_*',
#         'keep *_ecalRecHit_EcalRecHitsEE_*',
#         'keep *_towerMaker_*_*',
#         'keep recoSuperClusters_correctedHybridSuperClusters_*_*',
#         'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*',
#         'keep *_gsfElectronCores_*_*'
        'keep *_eleIsoDepositEcalFromHits_*_WWIso',
        'keep *_eleIsoDepositHcalFromTowers_*_WWIso',
        'keep *_muIsoDepositCalByAssociatorTowers_ecal_WWIso',
        'keep *_muIsoDepositCalByAssociatorTowers_hcal_WWIso',
        'keep recoPFCandidates_particleFlow__RECO',
        'keep *_addPileupInfo_*_*',
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 
        'wwAllPath',
    )),
)
process.outpath = cms.EndPath(process.out)



#  _______   _                         ______ _ _ _            
# |__   __| (_)                       |  ____(_) | |           
#    | |_ __ _  __ _  __ _  ___ _ __  | |__   _| | |_ ___ _ __ 
#    | | '__| |/ _` |/ _` |/ _ \ '__| |  __| | | | __/ _ \ '__|
#    | | |  | | (_| | (_| |  __/ |    | |    | | | ||  __/ |   
#    |_|_|  |_|\__, |\__, |\___|_|    |_|    |_|_|\__\___|_|   
#               __/ | __/ |                                    
#              |___/ |___/                                  
# 
#import HLTrigger.HLTfilters.hltHighLevel_cfi 
#process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone() 
#process.hltFilter.HLTPaths = [
#     "HLT_Ele15_LW_L1R",
#     "HLT_Mu9",
# ]
 
#process.primaryVertexFilter = cms.EDFilter("VertexSelector",
#    src = cms.InputTag("offlinePrimaryVertices"),
#    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), 
#    filter = cms.bool(True),   
# )
# 
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
# 
# process.load("WWAnalysis.Tools.incompleteECALReadout_cfi")
# 
# process.goodEvents = cms.Sequence(
#     process.hltFilter 
#     process.primaryVertexFilter 
#     process.noscraping 
#     process.incompleteECALReadout 
# )

#Trigger Matching (hopefully?!?)
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.onlyStandAlone = True
process.patTrigger.processName  = '*' 
process.patTriggerEvent.processName = '*'
#process.MessageLogger.suppressWarning += ['patTrigger'] # 1 warning per event on old runs!
    
#     process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

#  _____               _____ _    _           
# / ____|             / ____| |  (_)          
#| |  __  ___ _ __   | (___ | | ___ _ __ ___  
#| | |_ |/ _ \ '_ \   \___ \| |/ / | '_ ` _ \ 
#| |__| |  __/ | | |  ____) |   <| | | | | | |
# \_____|\___|_| |_| |_____/|_|\_\_|_| |_| |_|
#                                             

process.prunedGen = cms.EDProducer( "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *  ",
        "keep++ pdgId =   {Z0}",
        "++keep pdgId =   {Z0}",
        "keep++ pdgId =   {W+}",
        "++keep pdgId =   {W+}",
        "keep++ pdgId =   {W-}",
        "++keep pdgId =   {W-}",
        "keep++ pdgId =   {h0}",
        "++keep pdgId =   {h0}",
        "keep++ pdgId =   {e+}",
        "++keep pdgId =   {e+}",
        "keep++ pdgId =   {e-}",
        "++keep pdgId =   {e-}",
        "keep++ pdgId =  {mu+}",
        "++keep pdgId =  {mu+}",
        "keep++ pdgId =  {mu-}",
        "++keep pdgId =  {mu-}",
        "++keep pdgId =      5",
        "++keep pdgId =     -5",
        "++keep pdgId =      4",
        "++keep pdgId =     -4",
        "++keep pdgId =     12",
        "++keep pdgId =     14",
        "++keep pdgId =     16",
        "++keep pdgId =    -12",
        "++keep pdgId =    -14",
        "++keep pdgId =    -16",
        "++keep pdgId = {tau+}",
        "++keep pdgId = {tau-}",
    )
)

#  ______ _           _                     _____      _   _      
# |  ____| |         | |                   |  __ \    | | | |     
# | |__  | | ___  ___| |_ _ __ ___  _ __   | |__) |_ _| |_| |__   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  |  ___/ _` | __| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | | |  | (_| | |_| | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| |_|   \__,_|\__|_| |_| 
#                                                                 

process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")


import WWAnalysis.SkimStep.simpleCutBasedElectronIDSpring10_cfi as newcb
process.eidVBTFRel95 = newcb.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = newcb.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )

process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
process.liklihoodID = process.eidLikelihoodExt.clone() 

process.eidSequence = cms.Sequence(
    process.eidVBTFRel95 + 
    process.eidVBTFRel80 + 
    process.liklihoodID 
)

process.patElectrons.electronIDSources = cms.PSet()
process.patElectrons.electronIDSources.eidVBTFRel95 = cms.InputTag("eidVBTFRel95")
process.patElectrons.electronIDSources.eidVBTFRel80 = cms.InputTag("eidVBTFRel80")
process.patElectrons.electronIDSources.liklihoodID = cms.InputTag("liklihoodID")



ELE_ISO_CUT=("( isEB && (dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt)/pt < 0.5) ||"+
             "( isEE && (dr03TkSumPt +" +
             " dr03EcalRecHitSumEt +" +
             " dr03HcalTowerSumEt)/pt < 0.5) ")



process.eleTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patElectrons" ),
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
    src     = cms.InputTag(  "patElectrons" ),
    matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT'), cms.InputTag('eleIdTriggerMatchHLT'))
) 

process.goodElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string("pt > 10"),
)

import WWAnalysis.SkimStep.expectedHitsComputer_cfi
#process.expectedHitsEle = WWAnalysis.SkimStep.expectedHitsComputer_cfi.expectedHitsComputer.clone()
#process.expectedHitsEle.inputColl   = cms.InputTag("gsfElectrons")
#process.expectedHitsEle.useGsfTrack = cms.bool(True)

process.patElectrons.embedPFCandidate = False
process.patElectrons.embedSuperCluster = True
process.patElectrons.embedTrack = True
process.patElectrons.addElectronID = True
process.electronMatch.matched = "prunedGen"
process.patElectrons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("convValueMapProd","dist"),
    cms.InputTag("convValueMapProd","dcot"),
#     cms.InputTag("liklihoodID")
)
#process.patElectrons.userData.userInts.src = cms.VInputTag(
#    cms.InputTag('expectedHitsEle'),
#)


process.load("WWAnalysis.Tools.convValueMapProd_cfi")

if isMC: 
    process.preElectronSequence = cms.Sequence(
#        process.expectedHitsEle + 
        process.convValueMapProd + 
        process.eidSequence + 
        process.patTrigger + 
        process.prunedGen * 
        process.electronMatch 
    )
else:
    removeMCMatching(process, ['Electrons'])
    process.preElectronSequence = cms.Sequence(
#        process.expectedHitsEle + 
        process.convValueMapProd + 
        process.eidSequence + 
        process.patTrigger 
    )


process.elSeq = cms.Sequence( 
    process.preElectronSequence * 
    process.patElectrons *
    (process.eleTriggerMatchHLT + process.eleIdTriggerMatchHLT) *
    process.patElectronsWithTrigger *
    process.goodElectrons 
)

#  __  __                     _____      _   _
# |  \/  |                   |  __ \    | | | |
# | \  / |_   _  ___  _ __   | |__) |_ _| |_| |__
# | |\/| | | | |/ _ \| '_ \  |  ___/ _` | __| '_ \
# | |  | | |_| | (_) | | | | | |  | (_| | |_| | | |
# |_|  |_|\__,_|\___/|_| |_| |_|   \__,_|\__|_| |_|
#
MUON_ISO_CUT=("(isolationR03().emEt +" +
              " isolationR03().hadEt +" +
              " isolationR03().sumPt)/pt < 0.5")

#process.expectedHitsMu  = WWAnalysis.SkimStep.expectedHitsComputer_cfi.expectedHitsComputer.clone()

process.muonTriggerMatchHLT = cms.EDProducer( 'PATTriggerMatcherDRDPtLessByR',
    src     = cms.InputTag( 'patMuons' ),
    matched = cms.InputTag( 'patTrigger' ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltL3MuonCandidates' ),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.1 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)
process.patMuonsWithTrigger = cms.EDProducer( 'PATTriggerMatchMuonEmbedder',
    src     = cms.InputTag(  'patMuons' ),
    matches = cms.VInputTag('muonTriggerMatchHLT')
)

process.patMuons.embedPFCandidate = False
process.patMuons.embedTrack = True
#process.patMuons.userData.userInts.src = cms.VInputTag(
#    cms.InputTag('expectedHitsMu'),
#)
process.muonMatch.matched = "prunedGen"


if isMC: 
    if False: ## Turn this on to get extra info on muon MC origin, on GEN-SIM-RECO
        process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
        from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
        addClassByHits(process.patMuons, labels=['classByHitsGlb'], extraInfo = True)
        process.muonClassificationByHits = cms.Sequence(process.mix * process.trackingParticlesNoSimHits * process.classByHitsGlb )
        #process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger + process.expectedHitsMu + process.muonClassificationByHits)
        process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger + process.muonClassificationByHits)
        process.MessageLogger.suppressWarning += ['classByHitsGlb'] # kill stupid RPC hit associator warning
    else:
        #process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger + process.expectedHitsMu)
        process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger )
else:
    removeMCMatching(process, ['Muons'])
    #process.preMuonSequence = cms.Sequence(process.patTrigger + process.expectedHitsMu)
    process.preMuonSequence = cms.Sequence(process.patTrigger)


process.goodMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 10")
)

process.muSeq = cms.Sequence( 
    process.preMuonSequence *
    process.patMuons *
    process.muonTriggerMatchHLT *
    process.patMuonsWithTrigger *
    process.goodMuons 
)
#               _               _____      _ _           _   _                 
#    /\        | |             / ____|    | | |         | | (_)                
#   /  \  _   _| |_ _ __ ___  | |     ___ | | | ___  ___| |_ _  ___  _ __  ___ 
#  / /\ \| | | | __| '__/ _ \ | |    / _ \| | |/ _ \/ __| __| |/ _ \| '_ \/ __|
# / ____ \ |_| | |_| | |  __/ | |___| (_) | | |  __/ (__| |_| | (_) | | | \__ \
#/_/    \_\__,_|\__|_|  \___|  \_____\___/|_|_|\___|\___|\__|_|\___/|_| |_|___/
#                                                                              

from PhysicsTools.PatAlgos.tools.jetTools import *

process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")  
process.pfPileUp.Enable = True    
process.ak5PFJetsNoPU = process.ak5PFJets.clone(    
    src =   cms.InputTag("pfNoPileUp")  
)
process.jetSequence = cms.Sequence(   
    process.pfNoPileUpSequence  +   
    process.ak5PFJetsNoPU   
)  

if not isMC:
    removeMCMatching(process, ['Jets'])

addJetCollection(
    process,
    cms.InputTag("ak5PFJetsNoPU"),
    algoLabel    = "NoPU",
    typeLabel    = "",
    doJTA        = True,
    doBTagging   = True,
    #jetCorrLabel = ('AK5', 'PF'),
    jetCorrLabel = ('AK5PF',cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)

addJetCollection(
    process,
    cms.InputTag("ak5PFJets"),
    algoLabel    = "PF",
    typeLabel    = "",
    doJTA        = True,
    doBTagging   = True,
    #jetCorrLabel = ('AK5', 'PF'),
    jetCorrLabel = ('AK5PF',cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)

if not isMC:
    process.patJets.addGenJetMatch = False
    process.patJetsPF.addGenJetMatch = False
    process.patJetsPF.addGenPartonMatch = False
    process.patJetsPF.getJetMCFlavour = False
    process.patJetsNoPU.addGenJetMatch = False
    process.patJetsNoPU.addGenPartonMatch = False
    process.patJetsNoPU.getJetMCFlavour = False

process.patJets.addTagInfos = False
process.patJetsPF.addTagInfos = False
process.patJetsNoPU.addTagInfos = False
process.patJets.embedPFCandidates = False
process.patJetsPF.embedPFCandidates = False
process.patJetsNoPU.embedPFCandidates = False
process.patJets.addAssociatedTracks = False
process.patJetsPF.addAssociatedTracks = False
process.patJetsNoPU.addAssociatedTracks = False

process.patJetCorrFactors.levels = cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])

process.cleanPatJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("patJets"),
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("patElectronsWithTrigger"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False),
        ),
        mu = cms.PSet(
           src       = cms.InputTag("patMuonsWithTrigger"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False),
        ),
    ),
    finalCut = cms.string(''),
)

process.cleanPatJetsPF = process.cleanPatJets.clone( src = cms.InputTag("patJetsPF") )
process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("patJetsNoPU") )

process.noHypoGenParticles = cms.EDProducer( "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "keep  *  ",
        "drop++ pdgId = {W+} & status = 3",
        "drop++ pdgId = {W-} & status = 3",
    )
)

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.myGenParticlesForJetsNoNu = process.genParticlesForJetsNoNu.clone(src="noHypoGenParticles")
process.myGenJets = process.ak5GenJets.clone(src="myGenParticlesForJetsNoNu")
process.highPtGenJets = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("myGenJets"),
    cut = cms.string("pt()>10"),
    filter = cms.bool(False),
)


if isMC:
    process.autreSeq = cms.Sequence( 
        process.jetSequence * 
            process.makePatJets * (
            process.cleanPatJets + 
            process.cleanPatJetsPF + 
            process.cleanPatJetsNoPU 
        ) + (
            process.noHypoGenParticles *
            process.myGenParticlesForJetsNoNu *
            process.myGenJets *
            process.highPtGenJets
        )
    )
else:
    process.autreSeq = cms.Sequence( 
        process.jetSequence * 
            process.makePatJets * (
            process.cleanPatJets + 
            process.cleanPatJetsPF + 
            process.cleanPatJetsNoPU 
        ) 
    )

process.load("RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi")
process.muIsoDepositCalByAssociatorTowers.IOPSet.inputMuonCollection = 'goodMuons' 

process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequencePAT_cff")
process.eleIsoDepositEcalFromHits.src = 'goodElectrons'
process.eleIsoDepositHcalFromTowers.src = 'goodElectrons'
process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE")
process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB")
process.eleIsoDepositEcalFromHits.ExtractorPSet.energyMin = 0.0

# dZValues   = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00, 9999.0 ]
# d0Values   = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00, 9999.0 ]
# tkStrips   = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025]
# tkDeltaRs  = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025]
# outerDRs   = [0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]
# eneThresh  = [0.02 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
#               0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 ] 
# etThresh   = [0.02 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
#               0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 ] 


process.eleIso = cms.Sequence(
    process.eleIsoDepositEcalFromHits +
    process.eleIsoDepositHcalFromTowers +
    process.muIsoDepositCalByAssociatorTowers
)

#   _____ _    _             _____      _   _     
#  / ____| |  (_)           |  __ \    | | | |    
# | (___ | | ___ _ __ ___   | |__) |_ _| |_| |__  
#  \___ \| |/ / | '_ ` _ \  |  ___/ _` | __| '_ \ 
#  ____) |   <| | | | | | | | |  | (_| | |_| | | |
# |_____/|_|\_\_|_| |_| |_| |_|   \__,_|\__|_| |_|
#                                                 

process.goodLeps = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("goodMuons"), cms.InputTag("goodElectrons"))
)

process.diLep = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('goodLeps goodLeps'),
    cut = cms.string('daughter(0).pt > 10 && daughter(1).pt > 10'),
    checkCharge = cms.bool(False)
)

process.countWW  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("diLep"), minNumber = cms.uint32(1))
process.wwSeq = cms.Sequence( process.goodLeps * process.diLep * process.countWW )

if isAOD:
    process.wwAllPath = cms.Path( ( process.elSeq + process.muSeq ) * process.wwSeq + process.autreSeq )
else:
    process.wwAllPath = cms.Path( ( process.elSeq + process.muSeq ) * process.wwSeq + process.autreSeq + process.eleIso )

process.scrap = cms.Path( process.noscraping )

