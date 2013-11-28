import FWCore.ParameterSet.Config as cms

electronPFIsoMapProd = cms.EDProducer('ElectronPFIsoMapProd',
    eleLabel = cms.untracked.InputTag("gsfElectrons"),
    vtxLabel = cms.untracked.InputTag("goodPrimaryVertices"),
    pfLabel = cms.untracked.InputTag("particleFlow"),
    deltaR = cms.untracked.double(0.4)
)
