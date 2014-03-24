import FWCore.ParameterSet.Config as cms
rechitDumper = cms.EDAnalyzer("RecHitDumper",
                              EBRecHits = cms.InputTag("ecalRecHit:EcalRecHitsEB"),
                              EERecHits = cms.InputTag("ecalRecHit:EcalRecHitsEE"),
                              isMC = cms.untracked.bool(True),
                              nameFile = cms.untracked.string('tree.root'),
                              nameTree = cms.untracked.string('ntp0')
                              )


