import FWCore.ParameterSet.Config as cms

boostedElectrons = cms.EDProducer("PatElectronBooster", 
                                  src = cms.InputTag("patElectrons"), 
                                  vertices = cms.InputTag("offlinePrimaryVertices"),
                                  doTrigMVA = cms.bool(True),
                                  ebGlobalRecHitCollection = cms.InputTag("ecalRecHitGlobal","EcalRecHitsGlobalEB","reRECO"),
                                  eeGlobalRecHitCollection = cms.InputTag("ecalRecHitGlobal","EcalRecHitsGlobalEE","reRECO"),
                                  ebMultifit50nsRecHitCollection = cms.InputTag("ecalRecHitMultiFit50ns","EcalRecHitsMultiFit50nsEB","reRECO"),
                                  eeMultifit50nsRecHitCollection = cms.InputTag("ecalRecHitMultiFit50ns","EcalRecHitsMultiFit50nsEE","reRECO"),
                                  ebMultifit25nsRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB","reRECO"),
                                  eeMultifit25nsRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE","reRECO"),
                                  ebMultifitNoOOTPURecHitCollection = cms.InputTag("ecalRecHitMultiFitNoOOTPU","EcalRecHitsMultiFitNoOOTPUEB","reRECO"),
                                  eeMultifitNoOOTPURecHitCollection = cms.InputTag("ecalRecHitMultiFitNoOOTPU","EcalRecHitsMultiFitNoOOTPUEE","reRECO")
                                  )
