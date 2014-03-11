import FWCore.ParameterSet.Config as cms

# producer of rechits starting from digis
ecalAnalFitUncalibRecHit = cms.EDProducer("EcalUncalibRecHitProducer",
    EBdigiCollection = cms.InputTag("ecalDigis","ebDigis"),
    EEhitCollection = cms.string('EcalUncalibRecHitsEE'),
    EEdigiCollection = cms.InputTag("ecalDigis","eeDigis"),
    MinAmplEndcap = cms.double(14.0),
    MinAmplBarrel = cms.double(8.0),
    UseDynamicPedestal = cms.bool(True),
    ShapesFileName = cms.untracked.string("RecoLocalCalo/EcalRecAlgos/data/EcalShapes.root"),
    SavePlot = cms.untracked.bool(False),
    EBhitCollection = cms.string("EcalUncalibRecHitsEB"),
    algo = cms.string("EcalUncalibRecHitWorkerAnalFit")
)
