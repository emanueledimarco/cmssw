import FWCore.ParameterSet.Config as cms

isolatedElectronsRef = cms.EDFilter("HWWElectronIsolationRef",
                                    src = cms.InputTag("selectedElectrons"),
                                    trackIsolCut = cms.double(1.0),
                                    absolute = cms.bool(False)
                                    )
