from WWAnalysis.Tools.muonPFIsolations_cff import *
pfMuonIsolationSingleType = cms.Sequence( muonPFIsoChHad01 * muonPFIsoNHad01 * muonPFIsoPhoton01 *
                                          muonPFIsoChHad02 * muonPFIsoNHad02 * muonPFIsoPhoton02 *
                                          muonPFIsoChHad03 * muonPFIsoNHad03 * muonPFIsoPhoton03 *
                                          muonPFIsoChHad04 * muonPFIsoNHad04 * muonPFIsoPhoton04 *
                                          muonPFIsoChHad05 * muonPFIsoNHad05 * muonPFIsoPhoton05 *
                                          muonPFIsoChHad06 * muonPFIsoNHad06 * muonPFIsoPhoton06 *
                                          muonPFIsoChHad07 * muonPFIsoNHad07 * muonPFIsoPhoton07 )
