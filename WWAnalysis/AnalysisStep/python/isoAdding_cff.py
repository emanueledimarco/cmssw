import FWCore.ParameterSet.Config as cms

isoAddedMuons = cms.EDProducer("PatMuonIsoAdder",
  muonTag = cms.untracked.InputTag("boostedMuons"), 
  deposits = cms.VPSet()
)

isoAddedElectrons = cms.EDProducer("PatElectronIsoAdder",
    electronTag = cms.untracked.InputTag("boosteElectrons"),
    deposits = cms.VPSet(),
)
