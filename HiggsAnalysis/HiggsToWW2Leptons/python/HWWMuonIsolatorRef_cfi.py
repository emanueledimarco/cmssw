import FWCore.ParameterSet.Config as cms

isolatedMuonsRef = cms.EDFilter("HWWMuonIsolationRef",
                                src = cms.InputTag("muons"),
                                trackIsolCut = cms.double(1.0),
                                absolute = cms.bool(False)
                                )

