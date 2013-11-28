import FWCore.ParameterSet.Config as cms

# --- initial selection
eleSelection1 = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("gsfElectrons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && eta > -2.5 && eta < 2.5"),
)

muSelection1 = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && abs(eta) < 2.4"),
)

selection1Pairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("eleSelection1@+ muSelection1@-"),
    cut = cms.string("")
)

selection1PairsCounter = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("selection1Pairs"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


# --- IP selections
import WWAnalysis.Filters.eleIPSelection_cfi
eleIPSelection = WWAnalysis.Filters.eleIPSelection_cfi.eleIPSelection.clone(
  src = cms.InputTag("eleSelection1"),
)

import WWAnalysis.Filters.muIPSelection_cfi
muIPSelection = WWAnalysis.Filters.muIPSelection_cfi.muIPSelection.clone(
  src = cms.InputTag("muSelection1"),
)

ipSelPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("eleIPSelection@+ muIPSelection@-"),
    cut = cms.string("")
)

ipSelPairsCounter = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("ipSelPairs"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)





# --- Iso selections
ELE_ISO_CUT=("(dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt  )/pt < 0.10 ");

eleIsoSelection = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ISO_CUT),
)

MUON_ISO_CUT=("(isolationR03().emEt +" +
            " isolationR03().hadEt +" +
            " isolationR03().sumPt)/pt < 0.15 ");

muIsoSelection = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ISO_CUT),
)

isoSelPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("eleIsoSelection@+ muIsoSelection@-"),
    cut = cms.string("")
)

isoSelPairsCounter = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("isoSelPairs"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)




# --- ID selection
ELE_ID_CUT=("( (isEB && sigmaIetaIeta < 0.01 &&" +           
            " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
            " hadronicOverEm < 0.04) ||" +
            " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
            " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 &&" +
            " hadronicOverEm < 0.025) )"
            );

eleID = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ID_CUT),
)

MUON_ID_CUT=("(isGlobalMuon && isTrackerMuon &&" +
             " innerTrack.found >10 &&" +
             " globalTrack.normalizedChi2 <10 &&" +
             " standAloneMuon.hitPattern.numberOfValidMuonHits > 0)" 
             );

muID = cms.EDFilter("MuonSelector",
       src = cms.InputTag("muIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(MUON_ID_CUT),
)

idSelPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("eleID@+ muID@-"),
    cut = cms.string("")
)

idSelPairsCounter = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("idSelPairs"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


# --- conversion rejection  ---    
import WWAnalysis.Filters.convRejectionSelection_cfi      
eleConvRejection = WWAnalysis.Filters.convRejectionSelection_cfi.convRejectionSelection.clone(
  src = cms.InputTag("eleID"),
)

convSelPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("eleConvRejection@+ muID@-"),
    cut = cms.string("")
)

convSelPairsCounter = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("convSelPairs"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


# ---- diLepton selection
#diLeptSel1 = cms.EDFilter("CandViewSelector",
#    src = cms.InputTag("convSelPairs"),
#    filter = cms.bool(True),                              
#    cut = cms.string("mass > 12.0"),
#)

diLeptSel1 = cms.EDProducer("CandViewShallowCloneCombiner",
    decay       = cms.string("eleConvRejection@+ muID@-"), 
    checkCharge = cms.bool(False),           
    cut         = cms.string("(mass > 12.0)"),
)

diLeptFilter1 = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("diLeptSel1"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


diLeptSel2 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("diLeptSel1"),
    filter = cms.bool(True),                              
#    cut = cms.string("abs(mass-91.1876) > 15.0"),
    cut = cms.string("mass >=0.0"), #dummy
)

# --- Jet veto
#import WWAnalysis.Filters.jetVetoFilter_cfi      
#jetVetoLept = WWAnalysis.Filters.jetVetoFilter_cfi.jetVetoFilter.clone(
#   src = cms.InputTag("diLeptSel2"),
#)

#diLeptAfterJetVeto = cms.EDFilter("CandViewCountFilter", 
#  src = cms.InputTag("jetVetoLept"),
#  filter = cms.bool(True),                              
#  minNumber = cms.uint32(1),
#)


# --- MET cuts
metSelLept = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("tcMet"),
    filter = cms.bool(True),                              
    cut = cms.string("pt > 20.0"),
)

# --- soft muon veto
import WWAnalysis.Filters.softMuonVeto_cfi
leptSoftMuonVeto = WWAnalysis.Filters.softMuonVeto_cfi.softMuonVeto.clone(
  srcCompCands = cms.InputTag("diLeptSel1"),
)


# ---- extra lepton veto
import WWAnalysis.Filters.extraLeptonVeto_cfi      
extraLeptonVetoForLept = WWAnalysis.Filters.extraLeptonVeto_cfi.extraLeptonVeto.clone(
  srcCompCands = cms.InputTag("diLeptSel2"),
  srcMuons     = cms.InputTag("muForVeto"),
  srcElectrons = cms.InputTag("eleForVeto"),
)
