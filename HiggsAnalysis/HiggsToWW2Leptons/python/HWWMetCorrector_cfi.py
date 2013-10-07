import FWCore.ParameterSet.Config as cms

muonCorrectedMET = cms.EDProducer("HWWMetCorrector",
                                  CaloMetLabel = cms.InputTag("met"),
                                  MuonLabel = cms.InputTag("muons"),
                                  MuonCorrection = cms.bool(True),
                                  MuonDepositCor = cms.bool(True),
                                  MuonPtMin = cms.double(5.0),
                                  MuonEtaMax = cms.double(2.5),
                                  MuonTrackD0Max = cms.double(999.0),
                                  MuonTrackDzMax = cms.double(999.0),
                                  MuonNHitsMin = cms.double(5),
                                  MuonDPtMax = cms.double(0.5),
                                  MuonChiSqMax = cms.double(1000),
                                  CorrMetMin = cms.double(5.0)
                                  )
