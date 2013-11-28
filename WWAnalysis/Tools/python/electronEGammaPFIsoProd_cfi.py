from WWAnalysis.Tools.electronPFIsolations_cff import *
pfEGammaIsolationSingleType = cms.Sequence( electronPFIsoChHad01 * electronPFIsoNHad01 * electronPFIsoPhoton01 *
                                            electronPFIsoChHad02 * electronPFIsoNHad02 * electronPFIsoPhoton02 *
                                            electronPFIsoChHad03 * electronPFIsoNHad03 * electronPFIsoPhoton03 *
                                            electronPFIsoChHad04 * electronPFIsoNHad04 * electronPFIsoPhoton04 *
                                            electronPFIsoChHad05 * electronPFIsoNHad05 * electronPFIsoPhoton05 *
                                            electronPFIsoChHad06 * electronPFIsoNHad06 * electronPFIsoPhoton06 *
                                            electronPFIsoChHad07 * electronPFIsoNHad07 * electronPFIsoPhoton07 )
