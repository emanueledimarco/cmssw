import FWCore.ParameterSet.Config as cms

extraLeptonVeto = cms.EDFilter("ExtraLeptonVeto",
  filter = cms.bool(True),
  srcCompCands  = cms.InputTag("diMuonSel2"),
  srcMuons      = cms.InputTag("muons"),
  srcElectrons  = cms.InputTag("eles"),
  dR            = cms.double(0.1),                         
)
