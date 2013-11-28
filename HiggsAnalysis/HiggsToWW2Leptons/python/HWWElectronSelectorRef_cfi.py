import FWCore.ParameterSet.Config as cms

selectedElectronsRef = cms.EDFilter("HWWElectronSelectionRef",
                                    src                = cms.InputTag("gsfElectrons"),
				    hOverECut          = cms.double(0.150),
				    sigmaiEtaiEtaEBCut = cms.double(0.015),
				    sigmaiEtaiEtaEECut = cms.double(0.035),
				    deltaphiinCut      = cms.double(0.100),
				    deltaetainCut      = cms.double(0.010)
                                    )
