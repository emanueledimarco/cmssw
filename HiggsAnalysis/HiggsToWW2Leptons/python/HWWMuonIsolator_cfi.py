import FWCore.ParameterSet.Config as cms

isolatedMuons = cms.EDFilter("HWWMuonIsolation",
                             src = cms.InputTag("muons"),
                             trackIsolCut = cms.double(1.0),
                             absolute = cms.bool(False)
                             )

