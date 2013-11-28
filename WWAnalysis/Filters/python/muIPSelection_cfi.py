import FWCore.ParameterSet.Config as cms

muIPSelection = cms.EDFilter("MuIPSelection",
  src = cms.InputTag("muons"),
  filter = cms.bool(True),
  d0Cut = cms.double(0.020),                              
)
