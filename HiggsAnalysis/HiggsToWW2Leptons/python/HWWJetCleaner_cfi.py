import FWCore.ParameterSet.Config as cms

cleanedJets = cms.EDFilter("HWWJetCleaning",
                           src = cms.InputTag("iterativeCone5CaloJets"),
                           GSFInputTag = cms.InputTag("gsfElectrons"),
                           maxDeltaR = cms.double(0.2),
                           maxHadronicOverEm = cms.double(0.2),
                           minEleOJetEratio = cms.double(0.9)
                           )
