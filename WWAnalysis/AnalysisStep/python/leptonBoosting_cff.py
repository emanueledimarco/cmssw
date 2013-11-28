import FWCore.ParameterSet.Config as cms

boostedMuons = cms.EDProducer("PatMuonBooster",
  muonTag = cms.InputTag("cleanPatMuonsTriggerMatch"), 
  trackTag = cms.InputTag("generalTracks"),       
  vertexTag = cms.InputTag("goodPrimaryVertices"),
  deposits = cms.VPSet()
)

boostedElectrons = cms.EDProducer("PatElectronBooster",
    electronTag = cms.InputTag("cleanPatElectronsTriggerMatch"),
    trackTag = cms.InputTag("generalTracks"),       
    vertexTag = cms.InputTag("goodPrimaryVertices"),
    caloTowersTag = cms.InputTag("towerMaker"),
    electronIDSources = cms.PSet(),
    deposits = cms.VPSet(),
)
