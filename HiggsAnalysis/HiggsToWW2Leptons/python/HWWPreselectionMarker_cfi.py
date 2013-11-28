import FWCore.ParameterSet.Config as cms

preselectionMarker = cms.EDProducer("HWWPreselectionMarker",
                                    ElectronLabel = cms.InputTag("isolatedElectrons"),
                                    MuonLabel = cms.InputTag("isolatedMuons"),
                                    JetLabel = cms.InputTag("ak5CaloJets"),
                                    MetLabel = cms.InputTag("tcMet"),
                                    LeptonPtMinMin = cms.double(10.0),
                                    LeptonPtMaxMin = cms.double(20.0),
                                    LeptonEtaMax = cms.double(2.5),
                                    LeptonChargeCombination = cms.double(-1),
                                    MetMin = cms.double(30.0),
                                    InvMassMin = cms.double(12.0)
                                    )

