import FWCore.ParameterSet.Config as cms

muonPFIsoSingleTypeMapProd = cms.EDProducer('MuonPFIsoSingleTypeMapProd',
                                            muonLabel = cms.untracked.InputTag("muons"),
                                            vtxLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
                                            pfLabel = cms.untracked.InputTag("pfNoPileUp"),
                                            deltaR = cms.untracked.double(0.4),
                                            innerDeltaRVeto = cms.untracked.double(0.0),
                                            pfTypes = cms.untracked.vint32(1, 2, 3, 4, 5),
                                            directional = cms.untracked.bool(False)
                                            )
