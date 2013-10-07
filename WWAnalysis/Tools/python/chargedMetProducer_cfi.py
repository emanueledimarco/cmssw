import FWCore.ParameterSet.Config as cms

chargedMetProducer = cms.EDProducer('ChargedMetProducer',
    collectionTag = cms.InputTag("chargedParticleFlow"),
    vertexTag = cms.InputTag("offlinePrimaryVertices"),
    dzCut = cms.double(0.1),
    minNeutralPt = cms.double(4.0),
    maxNeutralEta = cms.double(3.0),
)

