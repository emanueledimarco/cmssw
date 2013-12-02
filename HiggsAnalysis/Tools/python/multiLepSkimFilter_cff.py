import FWCore.ParameterSet.Config as cms

vecbosPreFilterMuons = cms.EDFilter("MuonRefSelector",
                                   src = cms.InputTag("muons"),
                                   cut = cms.string("(isGlobalMuon || isTrackerMuon) && isPFMuon && track.isNonnull && pt > 3"),
                                   )

vecbosPreFilterElectrons = cms.EDFilter("GsfElectronRefSelector",
                                       src = cms.InputTag("gsfElectrons"),
                                       cut = cms.string("pt > 0"),
                                       )

vecbosPreFilterLeps = cms.EDProducer("CandViewMerger",
                                    src = cms.VInputTag(cms.InputTag("vecbosPreFilterMuons"),
                                                        cms.InputTag("vecbosPreFilterElectrons"),
                                                        )
                                    )


vecbosPreFilterDiLepFilter = cms.EDFilter("CandViewCountFilter",
                                         src = cms.InputTag("vecbosPreFilterLeps"),
                                         minNumber = cms.uint32(2),
                                         )


vecbosPreFilterDiLepSkimSeq  = cms.Sequence(
    (vecbosPreFilterMuons + vecbosPreFilterElectrons) *
    vecbosPreFilterLeps *
    vecbosPreFilterDiLepFilter
    )


