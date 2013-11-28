import FWCore.ParameterSet.Config as cms

softMuonVeto = cms.EDFilter("SoftMuonVeto",
  filter = cms.bool(True),
  srcCompCands  = cms.InputTag("diMuonSel1"),
#  srcCompCands  = cms.InputTag("diMuonSel2"),
  srcMuons      = cms.InputTag("muons"),
  minTrkHits    = cms.int32(10),                         
  minPt         = cms.double(3.0),                         
  d0Cut         = cms.double(0.2),                         
)
