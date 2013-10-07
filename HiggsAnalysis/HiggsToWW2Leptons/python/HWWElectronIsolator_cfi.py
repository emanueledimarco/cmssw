import FWCore.ParameterSet.Config as cms

isolatedElectrons = cms.EDFilter("HWWElectronIsolation",
                                 src = cms.InputTag("selectedElectrons"),
                                 trackIsolCut = cms.double(1.0),
                                 absolute = cms.bool(False)
                                 )
