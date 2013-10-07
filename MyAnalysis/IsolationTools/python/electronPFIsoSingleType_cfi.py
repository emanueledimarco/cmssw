import FWCore.ParameterSet.Config as cms

electronPFIsoSingleTypeMapProd = cms.EDProducer('ElectronPFIsoSingleTypeMapProd',
                                                eleLabel = cms.untracked.InputTag("gsfElectrons"),
                                                vtxLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
                                                pfLabel = cms.untracked.InputTag("pfNoPileUp"),
                                                deltaR = cms.untracked.double(0.4),
                                                neutralThreshold = cms.untracked.double(0.),
                                                pfTypes = cms.untracked.vint32(1, 2, 3, 4, 5),
                                                directional = cms.untracked.bool(False),
                                                radial = cms.untracked.bool(False)
                                                )
