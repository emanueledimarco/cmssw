import FWCore.ParameterSet.Config as cms

muonPFIsoMapProd = cms.EDProducer('MuonPFIsoMapProd',
    muonLabel = cms.untracked.InputTag("muons"),
    vtxLabel = cms.untracked.InputTag("goodPrimaryVertices"),
    pfLabel = cms.untracked.InputTag("particleFlow"),
    deltaR = cms.untracked.double(0.3)
)
