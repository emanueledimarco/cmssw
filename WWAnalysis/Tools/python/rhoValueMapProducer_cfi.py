import FWCore.ParameterSet.Config as cms

rhoValueMapProducer = cms.EDProducer('RhoValueMapProducer',
    leptonTag = cms.untracked.InputTag("muons"),
    rhoTag = cms.untracked.InputTag("kt6PFJets","rho"),
)

