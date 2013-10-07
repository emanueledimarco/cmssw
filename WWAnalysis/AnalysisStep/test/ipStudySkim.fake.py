import FWCore.ParameterSet.Config as cms

process = cms.Process("IPStudySkim")

isMC = True


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

## my stuff
#process.load("WWAnalysis.AnalysisStep.hits4MuonsSelection_cfi")
process.load("WWAnalysis.AnalysisStep.ipStudySelection_cfi")
process.load("SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi")
process.load('WWAnalysis.Filters.genFilterWW2L2Nu_cfi')


# process.GlobalTag.globaltag = 'START3X_V26::All'
process.GlobalTag.globaltag = 'START39_V8::All'



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/data/mangano/MC/WGamma/CMSSW_38X/WGToMuNuG_TuneD6T_7TeV-madgraph_Fall2010_RECO.root',
#        'file:/data/mangano/MC/WGamma/CMSSW_38X/WGToENuG_TuneD6T_7TeV-madgraph_Fall2010_RECO.root',
#        'file:/data/mangano/MC/WW/VVJetsTo4L_madgraph_RECO.root'
#        'file:/data/mangano/MC/WW/WJetsToLNu_madgraph_RECO.root'
#        'file:/data/mangano/MC/Zbb/Zbb2L0Jets.alpgen.Fall10.RECO.root'
#        'file:/data/mangano/MC/WW/WJetsToLnu.Winter10.root'
        'file:/data/mangano/MC/WW/VVJetsTo4L_madgraph_RECO.Winter10.root'
#        'file:/data/mangano/MC/WW/HWW160.Winter10.root'
    )
)


from glob import glob
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/oct29Skim/WWFull/*/*.root')]

process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('edmOut.el.root'),
    fileName = cms.untracked.string('edmOut.root'),
    outputCommands =  cms.untracked.vstring(
        'drop *',
#        'keep patMuons_muMatch_*_*',
        'keep patMuons_boostedMuons_*_*',
        'keep patMuons_muID_*_*',
        'keep patMuons_muISO_*_*',
        'keep patElectrons_boostedElectrons_*_*',
        'keep patElectrons_eleID_*_*',
        'keep patElectrons_eleISO_*_*',
        'keep patElectrons_eleNoConv_*_*',
        'keep *_prunedGen_*_*',
        'keep TrackingRecHitsOwned_electronGsfTracks__*',
        'keep recoGsfTrackExtras_electronGsfTracks__*',
        'keep recoTrackExtras_electronGsfTracks__*',
        'keep PileupSummaryInfo_addPileupInfo__*',
        'keep recoVertexs_*_*_IPStudySkim',
#        'keep SiPixelClusteredmNewDetSetVector_siPixelClusters__*',
#        'keep *_kinkFinderMu_*_*',
#        'keep *_genParticles_*_*',
#        'keep *_newGen_*_*',
#        'keep *_mergedtruthNoSimHits_*_*'
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 
            'pel',
            'pmu'
    )),
#    verbose = cms.untracked.bool(True)
)

process.outOverflow = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('overflow.root'),
    outputCommands =  cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = 
                                      cms.vstring('overflowSelectionEl',
                                                  'overflowSelectionMu')),
)

process.outpath = cms.EndPath(process.out)




#Trigger Matching (hopefully?!?)
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.onlyStandAlone = True
process.patTrigger.processName  = '*' 
process.patTriggerEvent.processName = '*'
process.MessageLogger.suppressWarning += ['patTrigger'] # 1 warning per event on old runs!
    
#     process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

#  _____               _____ _    _           
# / ____|             / ____| |  (_)          
#| |  __  ___ _ __   | (___ | | ___ _ __ ___  
#| | |_ |/ _ \ '_ \   \___ \| |/ / | '_ ` _ \ 
#| |__| |  __/ | | |  ____) |   <| | | | | | |
# \_____|\___|_| |_| |_____/|_|\_\_|_| |_| |_|
#                                             

process.newGen = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
#        particleTypes = cms.vstring("pi+"),       # also picks pi- (optional)
        filter        = cms.vstring("pt > 0.0"),  # just for testing
        genParticles   = cms.InputTag("genParticles") # original genParticle list
)   


process.prunedGen = cms.EDProducer( "GenParticlePruner",
    src = cms.InputTag("genParticles"),
##    src = cms.InputTag("newGen"),
    select = cms.vstring(
        "drop  *  ",
        "keep++ pdgId = {Z0}",
        "++keep pdgId = {Z0}",
        "keep++ pdgId = {W+}",
        "++keep pdgId = {W+}",
        "keep++ pdgId = {W-}",
        "++keep pdgId = {W-}",
        "keep++ pdgId = {h0}",
        "++keep pdgId = {h0}",
        "keep++ pdgId = {e+} & pt > 3",
        "++keep pdgId = {e+} & pt > 3",
        "keep++ pdgId = {e-} & pt > 3",
        "++keep pdgId = {e-} & pt > 3",
        "keep++ pdgId = {mu+}",
        "++keep pdgId = {mu+}",
        "keep++ pdgId = {mu-}",
        "++keep pdgId = {mu-}"
    )
)

#process.load("SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi")
#process.load("SimGeneral.MixingModule.mixNoPU_cfi")
#process.mergedtruthNoSimHits.pdgIdTP = cms.vint32(11,-11)

#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.genFromSim = cms.EDProducer("PATGenCandsFromSimTracksProducer",
#         src           = cms.InputTag("g4SimHits"),   ## use "famosSimHits" for FAMOS
#         setStatus     = cms.int32(-2),
#         particleTypes = cms.vstring("e+"),          ## picks also mu-, of course
#        filter        = cms.vstring("pt > 7.0"),     ## just for testing
#         makeMotherLink = cms.bool(True),
#         writeAncestors = cms.bool(True),             ## save also the intermediate GEANT ancestors of the muons
#         genParticles   = cms.InputTag("genParticles"),
# )




#  ______ _           _                     _____      _   _      
# |  ____| |         | |                   |  __ \    | | | |     
# | |__  | | ___  ___| |_ _ __ ___  _ __   | |__) |_ _| |_| |__   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  |  ___/ _` | __| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | | |  | (_| | |_| | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| |_|   \__,_|\__|_| |_| 
#                                                                 


ELE_ISO_CUT=("( isEB && (dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt)/pt < 0.5) ||"+
             "( isEE && (dr03TkSumPt +" +
             " dr03EcalRecHitSumEt +" +
             " dr03HcalTowerSumEt)/pt < 0.5) ")


process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")

#Code for the liklihood eID
process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
process.liklihoodID = process.eidLikelihoodExt.clone() 

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

process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(("pt > 10 && "+ELE_ISO_CUT)),
)

import WWAnalysis.SkimStep.expectedHitsComputer_cfi
process.expectedHitsEle = WWAnalysis.SkimStep.expectedHitsComputer_cfi.expectedHitsComputer.clone()
process.expectedHitsEle.inputColl   = cms.InputTag("gsfElectrons")
process.expectedHitsEle.useGsfTrack = cms.bool(True)

process.patElectrons.embedPFCandidate = False
process.patElectrons.embedSuperCluster = False
process.patElectrons.embedTrack = True
process.patElectrons.addElectronID = False


process.electronMatch.matched = "prunedGen"
### configuration for fake studies
#process.electronMatch.maxDeltaEta = cms.double(0.05)
#process.electronMatch.mcStatus = cms.vint32(1,8)
#
process.electronMatch.mcPdgId = cms.vint32(11)
process.electronMatch.mcStatus = cms.vint32(1)
process.electronMatch.resolveByMatchQuality = cms.bool(False)
process.electronMatch.checkCharge = cms.bool(True)
process.electronMatch.resolveAmbiguities = cms.bool(True)





process.patElectrons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("convValueMapProd","dist"),
    cms.InputTag("convValueMapProd","dcot"),
#    cms.InputTag("liklihoodID")
)


process.load("WWAnalysis.Tools.convValueMapProd_cfi")

if isMC: 
    process.preElectronSequence = cms.Sequence(
        process.convValueMapProd + 
#        process.liklihoodID +
        process.patTrigger + 
#        process.newGen*
        process.prunedGen * 
        process.electronMatch 
    )
else:
    removeMCMatching(process, ['Electrons'])
    process.preElectronSequence = cms.Sequence(
        process.convValueMapProd + 
#        process.liklihoodID +
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

process.muonMatch.matched = "prunedGen"


if isMC: 
    if False: ## Turn this on to get extra info on muon MC origin, on GEN-SIM-RECO
        process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
        from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
        addClassByHits(process.patMuons, labels=['classByHitsGlb'], extraInfo = True)
        process.muonClassificationByHits = cms.Sequence(process.mix * process.trackingParticlesNoSimHits * process.classByHitsGlb )
        process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger 
#                                               + process.muonClassificationByHits
                                               )
        process.MessageLogger.suppressWarning += ['classByHitsGlb'] # kill stupid RPC hit associator warning
    else:
        process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger)

else:
    removeMCMatching(process, ['Muons'])
    process.preMuonSequence = cms.Sequence(process.patTrigger )





process.muSeq = cms.Sequence( 
    process.preMuonSequence *
    process.patMuons *
    process.muonTriggerMatchHLT *
    process.patMuonsWithTrigger 
)

#############  VERTEX STUFF #######
#
#zGap=True
zGap=False

if zGap:
    process.goodPrimaryVertices = cms.EDProducer("PrimaryVertexProducer",
     PVSelParameters = cms.PSet(
        maxDistanceToBeam = cms.double(2) ## 0.1
     ),
     verbose = cms.untracked.bool(False),
     algorithm = cms.string('AdaptiveVertexFitter'),
     minNdof  = cms.double(0.0),
     TkFilterParameters = cms.PSet(
        algorithm=cms.string('filter'),
        maxNormalizedChi2 = cms.double(20.0),     # 
        minSiliconLayersWithHits = cms.int32(5), # >= 5
        minPixelLayersWithHits = cms.int32(2),   # >= 2 
        maxD0Significance = cms.double(100.0),     # keep most primary tracks
        minPt = cms.double(0.0),                 # better for softish events
        trackQuality = cms.string("any")
     ),
     beamSpotLabel = cms.InputTag("offlineBeamSpot"),
     # label of tracks to be used
     TrackLabel = cms.InputTag("generalTracks"),
     useBeamConstraint = cms.bool(False),


     # clustering
     TkClusParameters = cms.PSet(
       algorithm   = cms.string('gap'),
        TkGapClusParameters = cms.PSet( 
            zSeparation = cms.double(0.1) 
        )
     )
    )
else:
    process.goodPrimaryVertices = cms.EDProducer("PrimaryVertexProducer",
     verbose = cms.untracked.bool(False),
     algorithm = cms.string('AdaptiveVertexFitter'),
     TrackLabel = cms.InputTag("generalTracks"),
     useBeamConstraint = cms.bool(False),
     beamSpotLabel = cms.InputTag("offlineBeamSpot"),
     minNdof  = cms.double(0.0),
     PVSelParameters = cms.PSet(
        maxDistanceToBeam = cms.double(1.0)
     ),
     TkFilterParameters = cms.PSet(
        algorithm=cms.string('filter'),
        maxNormalizedChi2 = cms.double(5.0),
        minPixelLayersWithHits=cms.int32(2),
        minSiliconLayersWithHits = cms.int32(5),
        maxD0Significance = cms.double(5.0), 
        minPt = cms.double(0.0),
        trackQuality = cms.string("any")
     ),

     TkClusParameters = cms.PSet(
        algorithm   = cms.string("DA"),
        TkDAClusParameters = cms.PSet( 
            coolingFactor = cms.double(0.8),  #  rather slow annealing for now
#            Tmin = cms.double(9.),            #  end of annealing
            Tmin = cms.double(4.),            #  end of annealing
            vertexSize = cms.double(0.01)     #  ~ resolution
        )
     )
    )


#####
#####




#   _____      _              _       _      
#  / ____|    | |            | |     | |     
# | (___   ___| |__   ___  __| |_   _| | ___ 
#  \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \
#  ____) | (__| | | |  __/ (_| | |_| | |  __/
# |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|
#                                            


process.muFromW = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
#    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==13 && abs(mother.mother.pdgId)==24")
)

process.elFromW = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
#    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==11 && abs(mother.mother.pdgId)==24")
)

process.muFromTau = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
#    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==13 && abs(mother.pdgId)==15")
)

process.elFromTau = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
#    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==11 && abs(mother.pdgId)==15")
)

process.muFromZ = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
#    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==13 && abs(mother.mother.pdgId)==23")
)

process.elFromZ = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
#    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==11 && abs(mother.mother.pdgId)==23")
)

####
from PhysicsTools.PatAlgos.patSequences_cff import *
process.overFlowEl = selectedPatElectrons.clone()
process.overFlowEl.src = "boostedElectrons"
process.overFlowEl.filter = cms.bool(True)
process.overFlowEl.cut = ("abs(userFloat('ip2')/userFloat('ipErr2'))>9.999")
process.countOverFlowEl  = cms.EDFilter("CandViewCountFilter", 
   src = cms.InputTag("overFlowEl"), 
   minNumber = cms.uint32(1)
)

process.overFlowMu = selectedPatMuons.clone()
process.overFlowMu.src = "boostedMuons"
process.overFlowMu.filter = cms.bool(True)
process.overFlowMu.cut = ("abs(userFloat('ip2')/userFloat('ipErr2'))>9.999")
process.countOverFlowMu  = cms.EDFilter("CandViewCountFilter", 
   src = cms.InputTag("overFlowMu"), 
   minNumber = cms.uint32(1)
)
####

### for signal studies
#process.pmu = cms.Path(process.genFilterWW2L2Nu*process.muSeq*process.ipStudyMuSelectionSequence )
#process.pel = cms.Path(process.genFilterWW2L2Nu*process.elSeq*process.ipStudyEleSelectionSequence)
process.vertex = cms.Path(process.goodPrimaryVertices)
#process.pmu = cms.Path(process.muSeq*process.ipStudyMuSelectionSequence )
#process.pel = cms.Path(process.elSeq*process.ipStudyEleSelectionSequence)

process.overflowSelectionEl = cms.Path(process.elSeq*process.ipStudyEleSelectionSequence*
                                       process.overFlowEl*process.countOverFlowEl)

process.overflowSelectionMu = cms.Path(process.muSeq*process.ipStudyMuSelectionSequence*
                                       process.overFlowMu*process.countOverFlowMu)

### for fake studies
process.pmu = cms.Path(#process.newGen*
                       process.muSeq*~process.muFromW*~process.muFromTau*process.ipStudyMuSelectionSequence )
process.pel = cms.Path(#process.newGen*
                       process.elSeq*~process.elFromW*~process.elFromTau*process.ipStudyEleSelectionSequence)

### for fake from Zbb
#process.pmu = cms.Path(#process.newGen*
#                       process.muSeq*process.elFromZ*process.ipStudyMuSelectionSequence )
#process.pel = cms.Path(#process.newGen*
#                       process.elSeq*process.muFromZ*process.ipStudyEleSelectionSequence)


process.schedule = cms.Schedule(
    process.vertex,
    process.pmu,     
    process.pel,   
    process.outpath
 )



