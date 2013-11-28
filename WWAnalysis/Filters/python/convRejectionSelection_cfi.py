import FWCore.ParameterSet.Config as cms

convRejectionSelection   = cms.EDFilter("ConvRejectionSelection",
  src            = cms.InputTag("gsfElectrons"),
  srcTracks      = cms.InputTag("generalTracks"),                               
  filter         = cms.bool(True),
  distMax        = cms.double(0.02),                              
  dCotanThetaMax = cms.double(0.02),                              
)
