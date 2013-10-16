import FWCore.ParameterSet.Config as cms

goodPhotons = cms.EDFilter("PhotonSelector",
                             src = cms.InputTag("photons"),
                             cut = cms.string("pt > 10.0 & abs( eta ) < 3.0")
                             )


goodPFJets = cms.EDFilter("PFJetSelector",
                            src = cms.InputTag("ak5PFJets"),
                            cut = cms.string("pt > 15.0 & abs( eta ) < 5.0")
                            )

goodPFCHSJets = cms.EDFilter("PFJetSelector",
                               src = cms.InputTag("ak5PFNoPUJets"),
                               cut = cms.string("pt > 15.0 & abs( eta ) < 5.0")
                               )

goodGenJets = cms.EDFilter("GenJetRefSelector",
                             src = cms.InputTag("ak5GenJets"),
                             cut = cms.string("pt > 15.0 & abs( eta ) < 5.0")
                             )

collectionsFilterSequence = cms.Sequence(goodPhotons *
                                         goodPFJets *
                                         goodPFCHSJets *
                                         goodGenJets)
