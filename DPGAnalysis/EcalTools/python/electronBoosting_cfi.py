import FWCore.ParameterSet.Config as cms

boostedElectrons = cms.EDProducer("PatElectronBooster", 
                                  src = cms.InputTag("patElectrons"), 
                                  vertices = cms.InputTag("offlinePrimaryVertices"),
                                  doTrigMVA = cms.bool(True),
                                  ebGlobalRecHitCollection = cms.InputTag("ecalRecHitGlobal","EcalRecHitsGlobalEB","reRECO"),
                                  eeGlobalRecHitCollection = cms.InputTag("ecalRecHitGlobal","EcalRecHitsGlobalEE","reRECO"),
                                  ebMultifit50nsRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB","reRECO"),
                                  eeMultifit50nsRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE","reRECO"),
                                  ebMultifit25nsRecHitCollection = cms.InputTag("ecalRecHitMultiFit25ns","EcalRecHitsMultiFit25nsEB","reRECO"),
                                  eeMultifit25nsRecHitCollection = cms.InputTag("ecalRecHitMultiFit25ns","EcalRecHitsMultiFit25nsEE","reRECO"),
                                  ebMultifitNoOOTPURecHitCollection = cms.InputTag("ecalRecHitMultiFitNoOOTPU","EcalRecHitsMultiFitNoOOTPUEB","reRECO"),
                                  eeMultifitNoOOTPURecHitCollection = cms.InputTag("ecalRecHitMultiFitNoOOTPU","EcalRecHitsMultiFitNoOOTPUEE","reRECO")
                                  )
