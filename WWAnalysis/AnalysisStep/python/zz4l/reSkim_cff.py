import FWCore.ParameterSet.Config as cms

muons4skim = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string("(isTrackerMuon||isGlobalMuon) && abs(eta) < 2.4 && pt > 3"),
)
electrons4skim = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string("abs(eta) < 2.5 && pt > 5"),
)
leptons4skim = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag( cms.InputTag("muons4skim"),  cms.InputTag("electrons4skim"), )
)
dileptons4skim = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('leptons4skim leptons4skim'),
    cut = cms.string('deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.01'),
    checkCharge = cms.bool(False)
)

skim2010 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('min(daughter(0).pt,daughter(1).pt) > 10 && ' +
                     'max(daughter(0).pt,daughter(1).pt) > 20'),
    filter = cms.bool(True),
)

skim40 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('mass > 40 && daughter(0).pdgId = -daughter(1).pdgId'), ## and SF/OS
    filter = cms.bool(True),
)
skim40NoOF  = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('mass > 40 && abs(daughter(0).pdgId) == abs(daughter(1).pdgId)'), ## and SF only
    filter = cms.bool(True),
)


reskim     = cms.Sequence(muons4skim + electrons4skim + leptons4skim + dileptons4skim + skim2010)
reskimNoOS = cms.Sequence(muons4skim + electrons4skim + leptons4skim + dileptons4skim + skim2010)


from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *
triggerResultsFilter.l1tResults = ''
triggerResultsFilter.throw = True
triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
triggerFilter7TeV_MC = triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu17_Mu8_v*',  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 'HLT_TripleEle10_CaloIdL_TrkIdVL_v*', "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*", "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*"])
triggerFilter7TeV_MC_Mu = triggerFilter7TeV_MC.clone(triggerConditions = [ 'HLT_Mu17_Mu8_v*' ])
triggerFilter7TeV_MC_El = triggerFilter7TeV_MC.clone(triggerConditions = [ 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 'HLT_TripleEle10_CaloIdL_TrkIdVL_v*' ])

### Trigger emulation. Note that run ranges inside might have to be changed
triggerFilter7TeV_DATA = cms.EDFilter("GioTriggerEmulator",
    trigger = cms.InputTag("TriggerResults","","HLT"),
    muons   = cms.InputTag("boostedMuons"),
    doubleMu = cms.string("Mu17_Mu8"), # or "Mu13_Mu8" or "none" or "any" (=loosest unprescaled)
    doubleEl = cms.string("Ele17_Ele8_Tight"), # or "any" or "none"
    tripleEl = cms.string("HLT_TripleEle10_CaloIdL_TrkIdVL_v"), # a trigger, or "none"
    mueg = cms.string("HLT_Mu8_17_Ele17_8"), # a trigger, or "none"
    runForMC = cms.uint32(1), # run whose trigger table should be assumed in MC (180000 is good for Fall11 MC; use 1 or 178419 to test emulator on DoubleMu7 or Mu13/Mu7)
)
triggerFilter7TeV_DATA_Mu = triggerFilter7TeV_DATA.clone(doubleEl = "none", tripleEl = "none")
triggerFilter7TeV_DATA_El = triggerFilter7TeV_DATA.clone(doubleMu = "none")
triggerFilter7TeV_DATA_El_Only = cms.Sequence(triggerFilter7TeV_DATA_El + ~triggerFilter7TeV_DATA_Mu)

## 2012
ElTrig2012 = [ 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 'HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*' ]
MuTrig2012 = [ 'HLT_Mu17_Mu8_v*', 'HLT_Mu17_TkMu8_v*' ]
EMTrig2012 = [ 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*' ]
triggerFilter8TeV = triggerResultsFilter.clone(triggerConditions = ElTrig2012 + MuTrig2012 + EMTrig2012)
triggerFilter8TeV_Mu = triggerResultsFilter.clone(triggerConditions = MuTrig2012)
triggerFilter8TeV_El = triggerResultsFilter.clone(triggerConditions = ElTrig2012)
triggerFilter8TeV_ME = triggerResultsFilter.clone(triggerConditions = EMTrig2012)
triggerFilter8TeV_El_Only = cms.Sequence(triggerFilter8TeV_El + ~triggerFilter8TeV_Mu)
triggerFilter8TeV_ME_Only = cms.Sequence(triggerFilter8TeV_ME + ~triggerFilter8TeV_Mu + ~triggerFilter8TeV_El)
