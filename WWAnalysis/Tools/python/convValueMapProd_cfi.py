import FWCore.ParameterSet.Config as cms

convValueMapProd = cms.EDProducer('ConvValueMapProd',
    gsfLabel = cms.untracked.InputTag("gsfElectrons"),
    tkLabel = cms.untracked.InputTag("generalTracks"),
    beamSpotLabel = cms.untracked.InputTag("offlineBeamSpot"),
    conversionLabel = cms.untracked.InputTag("trackerOnlyConversions"), #41X
#     conversionLabel = cms.untracked.InputTag("allConversions"), #42X
)
