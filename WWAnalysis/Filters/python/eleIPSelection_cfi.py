import FWCore.ParameterSet.Config as cms

eleIPSelection = cms.EDFilter("GsfEleIPSelection",
  src = cms.InputTag("gsfElectrons"),
  filter = cms.bool(True),
  d0Cut = cms.double(0.020),                              
)
