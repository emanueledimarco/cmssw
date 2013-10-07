import FWCore.ParameterSet.Config as cms

metFilter = cms.EDFilter("METFilter",
  srcMET = cms.InputTag("tcMet"),
  srcLeptons = cms.InputTag("daugtherListCleaner"),
  filter = cms.bool(True),
  projectedMETCut = cms.double(20.0),                              
  METCut = cms.double(45.0),                              
  skipMEtCut = cms.bool(False)
)
