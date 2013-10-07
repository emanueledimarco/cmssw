import FWCore.ParameterSet.Config as cms

# --- initial selection
muPlusSelection1 = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && abs(eta) < 2.4 && charge==+1"),
)

muMinusSelection1 = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && abs(eta) < 2.4 && charge==-1"),
)

muForVetoFirstSel = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muons"),
       filter = cms.bool(False),                              
       cut = cms.string("pt> 10.0 && abs(eta) < 2.4 "),
)



# --- IP selections
import WWAnalysis.Filters.muIPSelection_cfi
muPlusIPSelection = WWAnalysis.Filters.muIPSelection_cfi.muIPSelection.clone(
  src = cms.InputTag("muPlusSelection1"),
)

muMinusIPSelection = WWAnalysis.Filters.muIPSelection_cfi.muIPSelection.clone(
  src = cms.InputTag("muMinusSelection1"),
)

muForVetoIPSel = WWAnalysis.Filters.muIPSelection_cfi.muIPSelection.clone(
  src = cms.InputTag("muForVetoFirstSel"),
  filter = cms.bool(False),                              
)




# --- Iso selections
MUON_ISO_CUT=("(isolationR03().emEt +" +
            " isolationR03().hadEt +" +
            " isolationR03().sumPt)/pt < 0.15 ");

muPlusIsoSelection = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muPlusIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ISO_CUT),
)

muMinusIsoSelection = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muMinusIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ISO_CUT),
)

muForVetoIsoSel = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muForVetoIPSel"),
       filter = cms.bool(False),                              
       cut = cms.string(MUON_ISO_CUT),
)




# --- ID selection
MUON_ID_CUT=("(isGlobalMuon && isTrackerMuon &&" +
             " innerTrack.found >10 &&" +
             " globalTrack.normalizedChi2 <10 &&" +
             " standAloneMuon.hitPattern.numberOfValidMuonHits > 0)" 
             );

muPlusID = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muPlusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ID_CUT),
)

muMinusID = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muMinusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ID_CUT),
)

muForVeto = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muForVetoIsoSel"),
       filter = cms.bool(False),                              
       cut = cms.string(MUON_ID_CUT),
)
                           

# --- diLepton selection
diMuonSel1 = cms.EDProducer("CandViewShallowCloneCombiner",
    decay       = cms.string("muPlusID muMinusID"), 
    checkCharge = cms.bool(False),           
    cut         = cms.string("(mass > 12.0)"),
)

diMuonFilter1 = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("diMuonSel1"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


diMuonSel2 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("diMuonSel1"),
    filter = cms.bool(True),                              
    cut = cms.string("abs(mass-91.1876) > 15.0"),
)


# --- MET cuts
metSel1 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("tcMet"),
    filter = cms.bool(True),                              
    cut = cms.string("pt > 20.0"),
)


# ---- extra lepton veto
muForVetoSequence = cms.Sequence(muForVetoFirstSel*muForVetoIPSel*muForVetoIsoSel*muForVeto)

import WWAnalysis.Filters.extraLeptonVeto_cfi      
extraLeptonVetoForMu = WWAnalysis.Filters.extraLeptonVeto_cfi.extraLeptonVeto.clone(
  srcCompCands = cms.InputTag("diMuonSel2"),
  srcMuons     = cms.InputTag("muForVeto"),
  srcElectrons = cms.InputTag("eleForVeto"),
)
