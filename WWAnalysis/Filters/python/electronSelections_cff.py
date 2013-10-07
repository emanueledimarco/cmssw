import FWCore.ParameterSet.Config as cms

# --- initial selection
elePlusSelection1 = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("gsfElectrons"),
       filter = cms.bool(True),                              
#       cut = cms.string("pt> 20.0 && abs(eta) < 2.5 && charge==+1"),
       cut = cms.string("pt> 20.0 && eta > -2.5 && eta < 2.5 && charge==+1"),
)

eleMinusSelection1 = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("gsfElectrons"),
       filter = cms.bool(True),                              
       cut = cms.string("pt> 20.0 && eta > - 2.5 && eta < 2.5 && charge==-1"),
)

eleForVetoFirstSel = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("gsfElectrons"),
       filter = cms.bool(False),                              
       cut = cms.string("pt> 10.0 && eta > - 2.5 && eta < 2.5"),
)



# --- IP selections
import WWAnalysis.Filters.eleIPSelection_cfi
elePlusIPSelection = WWAnalysis.Filters.eleIPSelection_cfi.eleIPSelection.clone(
  src = cms.InputTag("elePlusSelection1"),
)

eleMinusIPSelection = WWAnalysis.Filters.eleIPSelection_cfi.eleIPSelection.clone(
  src = cms.InputTag("eleMinusSelection1"),
)

eleForVetoIPSel = WWAnalysis.Filters.eleIPSelection_cfi.eleIPSelection.clone(
  src = cms.InputTag("eleForVetoFirstSel"),
  filter = cms.bool(False),                              
)




# --- Iso selections
ELE_ISO_CUT=("(dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt  )/pt < 0.10 ");

elePlusIsoSelection = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("elePlusIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ISO_CUT),
)

eleMinusIsoSelection = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleMinusIPSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ISO_CUT),
)

eleForVetoIsoSel = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleForVetoIPSel"),
       filter = cms.bool(False),                              
       cut = cms.string(ELE_ISO_CUT),
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

elePlusID = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("elePlusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ID_CUT),
)

eleMinusID = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleMinusIsoSelection"),
       filter = cms.bool(True),                              
       cut = cms.string(ELE_ID_CUT),
)

eleForVetoIDSel = cms.EDFilter("GsfElectronSelector",
       src = cms.InputTag("eleForVetoIsoSel"),
       filter = cms.bool(False),                              
       cut = cms.string(ELE_ID_CUT),
)


# --- conversion rejection  ---    
import WWAnalysis.Filters.convRejectionSelection_cfi      
convRejectionPlus = WWAnalysis.Filters.convRejectionSelection_cfi.convRejectionSelection.clone(
  src = cms.InputTag("elePlusID"),
)

convRejectionMinus = WWAnalysis.Filters.convRejectionSelection_cfi.convRejectionSelection.clone(
  src = cms.InputTag("eleMinusID"),
)

eleForVeto = WWAnalysis.Filters.convRejectionSelection_cfi.convRejectionSelection.clone(
  src = cms.InputTag("eleForVetoIDSel"),
  filter = cms.bool(False),                              
)




# --- diLepton selection
diEleSel1 = cms.EDProducer("CandViewShallowCloneCombiner",
    decay       = cms.string("convRejectionPlus convRejectionMinus"), 
    checkCharge = cms.bool(False),           
    cut         = cms.string("(mass > 12.0)"),
)

diEleFilter1 = cms.EDFilter("CandViewCountFilter", 
  src = cms.InputTag("diEleSel1"),
  filter = cms.bool(True),                              
  minNumber = cms.uint32(1),
)


diEleSel2 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("diEleSel1"),
    filter = cms.bool(True),                              
    cut = cms.string("abs(mass-91.1876) > 15.0"),
)

# --- Jet veto
#import WWAnalysis.Filters.jetVetoFilter_cfi      
#jetVetoEle = WWAnalysis.Filters.jetVetoFilter_cfi.jetVetoFilter.clone(
#   src = cms.InputTag("diEleSel2"),
#)

#diEleAfterJetVeto = cms.EDFilter("CandViewCountFilter", 
#  src = cms.InputTag("jetVetoEle"),
#  filter = cms.bool(True),                              
#  minNumber = cms.uint32(1),
#)


# --- MET cuts
metSel1 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("tcMet"),
    filter = cms.bool(True),                              
    cut = cms.string("pt > 20.0"),
)

# --- soft muon veto
import WWAnalysis.Filters.softMuonVeto_cfi
eleSoftMuonVeto = WWAnalysis.Filters.softMuonVeto_cfi.softMuonVeto.clone(
  srcCompCands = cms.InputTag("diEleSel1"),
#  srcCompCands = cms.InputTag("diEleSel2"),
)


# ---- extra lepton veto
eleForVetoSequence = cms.Sequence(eleForVetoFirstSel*eleForVetoIPSel*eleForVetoIsoSel*eleForVetoIDSel*eleForVeto)

import WWAnalysis.Filters.extraLeptonVeto_cfi      
extraLeptonVetoForEle = WWAnalysis.Filters.extraLeptonVeto_cfi.extraLeptonVeto.clone(
  srcCompCands = cms.InputTag("diEleSel2"),
  srcMuons     = cms.InputTag("muForVeto"),
  srcElectrons = cms.InputTag("eleForVeto"),
)
