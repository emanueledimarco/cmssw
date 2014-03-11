import FWCore.ParameterSet.Config as cms

# producer of rechits starting from digis
ecalFixedAlphaBetaGammaFitUncalibRecHit = cms.EDProducer("EcalUncalibRecHitProducer",
    EEdigiCollection = cms.InputTag("ecalDigis","eeDigis"),
    EBdigiCollection = cms.InputTag("ecalDigis","ebDigis"),
    EEhitCollection = cms.string("EcalUncalibRecHitsEE"),
    alphaEB = cms.double(1.138),
    alphaEE = cms.double(1.890),
    betaEB = cms.double(1.655),
    betaEE = cms.double(1.400),
    AlphaBetaFilename = cms.untracked.string("NOFILE"),
    MinAmplEndcap = cms.double(14.0),
    MinAmplBarrel = cms.double(8.0),
    UseDynamicPedestal = cms.bool(True),
    PedRMSBarrel = cms.double(1.3),
    PedRMSEndcap = cms.double(2.2),
    EBhitCollection = cms.string("EcalUncalibRecHitsEB"),
    algo = cms.string("EcalUncalibRecHitWorkerFixedAlphaBetaGammaFit")
)
