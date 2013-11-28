import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.trigTools import *

eleTriggerColls = [ 
    [ 'hltL1IsoRecoEcalCandidate',       'hltL1NonIsoRecoEcalCandidate' ],
    [ 'hltPixelMatchElectronsL1Iso',     'hltPixelMatchElectronsL1NonIso' ],
    [ 'hltPixelMatch3HitElectronsL1Iso', 'hltPixelMatch3HitElectronsL1NonIso' ],
    [ 'hltPixelMatchElectronsActivity' ] ,
    [ 'hltPixelMatch3HitElectronsActivity' ] ,
    [ 'hltRecoEcalSuperClusterActivityCandidate' ] ,
]
eleTriggers = [
    "HLT_Ele10_LW_L1R_v*", "HLT_Ele15_LW_L1R_v*", "HLT_Ele15_SW_L1R_v*", 
    "HLT_Ele15_SW_CaloEleId_L1R_v*", "HLT_Ele17_SW_CaloEleId_L1R_v*", 
    "HLT_Ele17_SW_TightEleId_L1R_v*", "HLT_Ele17_SW_TighterEleIdIsol_L1R_v*",
    "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*", "HLT_Mu11_Ele8_v*", "HLT_Mu5_Ele17_v*",
    "HLT_Ele8_CaloIdL_CaloIsoVL_v*", "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
    "HLT_Ele17_CaloIdL_CaloIsoVL_v*", "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
    "HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*", "HLT_Ele45_CaloIdVT_TrkIdT_v*",
    "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*", "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*",
    "HLT_Mu17_Ele8_CaloIdL_v*", "HLT_Mu8_Ele17_CaloIdL_v*",
    "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*",
    "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*",
    "HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v*", "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*",
    "HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v*",
]
muTriggers = [
    "HLT_DoubleMu3_v*",
    "HLT_Mu5_v*",   "HLT_IsoMu12_v*",  "HLT_Mu11_Ele8_v*",         "HLT_Mu8_Jet40_v*",
    "HLT_Mu8_v*",   "HLT_IsoMu15_v*",  "HLT_Mu5_Ele17_v*",         "HLT_Mu8_Photon20_CaloIdVT_IsoT_v*",
    "HLT_Mu12_v*",  "HLT_IsoMu17_v*",  "HLT_Mu17_Ele8_CaloIdL_v*", "HLT_IsoMu12_LooseIsoPFTau10_v*",
    "HLT_Mu15_v*",  "HLT_IsoMu24_v*",  "HLT_Mu8_Ele17_CaloIdL_v*", "HLT_Mu15_LooseIsoPFTau20_v*",
    "HLT_Mu20_v*",  "HLT_IsoMu30_v*",                              "HLT_DoubleMu7_v*",
    "HLT_Mu21_v*",                
    "HLT_Mu24_v*",                
    "HLT_Mu30_v*",                
]

def addTriggerPaths(process):
    process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi" )
    process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi" )
    process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerMatchEmbedder_cfi" )
    process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi" )
    
    process.patL1Trigger = cms.EDProducer("L1Extra2PAT", src = cms.InputTag("l1extraParticles"))
    process.patDefaultSequence.replace(process.patTrigger, process.patTrigger + process.patL1Trigger)
    
    tempProd = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
        matchedCuts = cms.string('path( "HLT_Mu20_v*" )'),
        src = cms.InputTag("cleanPatMuons"),
        maxDPtRel = cms.double(0.5),
        resolveByMatchQuality = cms.bool(True),
        maxDeltaR = cms.double(0.5),
        resolveAmbiguities = cms.bool(True),
        matched = cms.InputTag("patTrigger")
    )
    tempProdNoPt = tempProd.clone()
    setattr(tempProdNoPt, '_TypedParameterizable__type', 'PATTriggerMatcherDRLessByR')
    delattr(tempProdNoPt, 'maxDPtRel')
    
    eleTriggerModules = dict(zip([ "cleanElectronTriggerMatch{0}".format(k.replace('v*','').replace('HLT_','').replace('_','')) for k in eleTriggers ],eleTriggers))
    for key in eleTriggerModules:
        setattr(process,key,tempProd.clone(src = "cleanPatElectrons", matchedCuts = 'path("{0}")'.format(eleTriggerModules[key])))
    
    eleTriggerCollModules = dict(zip([ "cleanElectronTriggerMatch{0}".format(k[0]) for k in eleTriggerColls ],eleTriggerColls))
    for key in eleTriggerCollModules:
        setattr(process,key,tempProd.clone(src = "cleanPatElectrons", matchedCuts = " || ".join(['coll("%s")' % c for c in eleTriggerCollModules[key]])))
    
    process.cleanElectronTriggerMatchL1 = tempProdNoPt.clone(
        src = "cleanPatElectrons",
        matched = "patL1Trigger",
        matchedCuts = "coll('l1extraParticles:Isolated') || coll('l1extraParticles:NonIsolated')",
    )
    
    muTriggerModules = dict(zip([ "cleanMuonTriggerMatch{0}".format(k.replace('v*','').replace('HLT_','').replace('_','')) for k in muTriggers ],muTriggers))
    for key in muTriggerModules:
        setattr(process,key,tempProd.clone(src = "cleanPatMuons", matchedCuts = 'path("{0}")'.format(muTriggerModules[key])))
    
    process.cleanMuonTriggerMatchByObject   = tempProd.clone(src = "cleanPatMuons", matchedCuts = 'coll("hltL3MuonCandidates")')
    process.cleanMuonTriggerMatchByTkObject   = tempProd.clone(src = "cleanPatMuons", matchedCuts = 'coll("hltGlbTrkMuonCands")')
    process.cleanMuonTriggerMatchByL2Object = tempProdNoPt.clone(src = "cleanPatMuons", matchedCuts = 'coll("hltL2MuonCandidates")', maxDeltaR=0.7)
    
    myDefaultTriggerMatchers = eleTriggerModules.keys()[:] + eleTriggerCollModules.keys()[:] + muTriggerModules.keys()[:] + [
        'cleanMuonTriggerMatchByObject',
        'cleanMuonTriggerMatchByTkObject',
        'cleanMuonTriggerMatchByL2Object',
    ]
    
    switchOnTriggerMatchEmbedding(process,triggerMatchers=myDefaultTriggerMatchers)
    
    process.patTrigger.onlyStandAlone = True
    process.patTrigger.processName  = '*' 
    
    process.noscraping = cms.EDFilter("FilterOutScraping",
        applyfilter = cms.untracked.bool(True),
        debugOn = cms.untracked.bool(False),
        numtrack = cms.untracked.uint32(10),
        thresh = cms.untracked.double(0.25)
    )
