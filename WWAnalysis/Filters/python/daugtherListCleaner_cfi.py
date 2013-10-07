import FWCore.ParameterSet.Config as cms

daugtherListCleaner = cms.EDProducer("DaugtherListCleaner",
  src = cms.InputTag("diMuonSel2"),
)
