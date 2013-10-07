import FWCore.ParameterSet.Config as cms

interestingVertexRefProducer = cms.EDProducer('InterestingVertexRefProducer',
    leptonTags = cms.VInputTag("muons","gsfElectrons"),
    vertexTag = cms.InputTag("offlinePrimaryVertices"),
)

