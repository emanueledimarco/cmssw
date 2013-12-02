import FWCore.ParameterSet.Config as cms

scaledMuons = cms.EDProducer("MuScleFitMuonCorrector",
                             src = cms.InputTag("muons"),
                             debug = cms.bool(False),
                             identifier = cms.string("Data2012_53X_ReReco"),
                             applySmearing = cms.bool(False),
                             fakeSmearing = cms.bool(False)
                             )

