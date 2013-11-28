import FWCore.ParameterSet.Config as cms

betaValueMapProducer = cms.EDProducer('BetaValueMapProducer',
    leptonTag = cms.untracked.InputTag("muons"),
    trackTag = cms.untracked.InputTag("generalTracks"),
    dEtaVeto = cms.untracked.double(0.0),
    dRVeto = cms.untracked.double(0.01),
    outerCone = cms.untracked.double(0.3),
    ptCut = cms.untracked.double(0.0),
)

