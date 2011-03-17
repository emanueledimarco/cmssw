import FWCore.ParameterSet.Config as cms

treeDumper = cms.EDAnalyzer("HWWTreeDumper",
                            HLTObjectsInfo = cms.untracked.PSet(triggerResults = cms.InputTag("TriggerResults","","AUTO"),
                                                                processName = cms.string("AUTO"),
                                                                triggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD","","HLT")
                                                                ),
                            electronCollection = cms.InputTag("ambiguityResolvedElectrons"),
                            #pflowElectronCollection = cms.InputTag("particleFlow","electrons"),
                            pflowElectronCollection = cms.InputTag("particleFlow"),
                            photonCollection = cms.InputTag("photons"),
                            muonCollection = cms.InputTag("muons"),
                            pfTauCollection = cms.InputTag("shrinkingConePFTauProducer"),
                            ecalSCCollection = cms.InputTag("mergedSuperClusters"),
                            ecalBarrelSCCollection = cms.InputTag("correctedHybridSuperClusters"),
                            ecalEndcapSCCollection = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapSuperClusters"),
                            ecalPFClusterCollection = cms.InputTag("pfElectronTranslator","pf"),
                            ecalBCCollection = cms.InputTag("mergedBasicClusters"),
                            ecalBarrelRecHits = cms.InputTag("reducedEcalRecHitsEB"),
                            ecalEndcapRecHits = cms.InputTag("reducedEcalRecHitsEE"),
                            calotowersForIsolationProducer = cms.InputTag("towerMaker"),
                            trackCollection = cms.InputTag("generalTracks"),
                            gsfTrackCollection = cms.InputTag("electronGsfTracks"),
                            globalMuonTrackCollection = cms.InputTag("globalMuons"),
                            standAloneMuonTrackCollection = cms.InputTag("standAloneMuons"),
                            #  refittedForDeDxTrackCollection = cms.InputTag("RefitterForDeDx"),
                            refittedForDeDxTrackCollection = cms.InputTag("generalTracks"),
                            vertexCollection = cms.InputTag("offlinePrimaryVertices"),
                            K0sCollection = cms.InputTag("generalV0Candidates","Kshort"),
                            genJetCollection = cms.InputTag("ak5GenJets"),
                            jetCollection1 = cms.InputTag("ak5CaloJetsL2L3Residual"),
                            jetCollection2 = cms.InputTag("ak5CaloJets"),
                            JPTjetCollection1 = cms.InputTag("ak5JPTJetsL2L3Residual"),
                            JPTjetCollection2 = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
                            PFjetCollection1 = cms.InputTag("ak5PFJetsL2L3Residual"),
                            PFjetCollection2 = cms.InputTag("ak5PFJets"),
                            PFpuCorrJetCollection1 = cms.InputTag("ak5PFJetsL1FastL2L3Residual"),
                            PFpuCorrJetCollection2 = cms.InputTag("ak5PFJetsL1FastL2L3Residual"), # no module with L1 but no L2L3
                            metCollection = cms.InputTag("met"), # preselection
                            # corrmetCollection = cms.InputTag("metMuonJESCorAK5"), # type I and II corr applied
                            TCmetCollection = cms.InputTag("tcMet"),
                            genMetCollection = cms.InputTag("genMetCalo"),
                            PFmetCollection = cms.InputTag("pfMet"),
                            PFpreIdCollection = cms.InputTag("trackerDrivenElectronSeeds:preid"),
                            calotowerCollection = cms.InputTag("lowThrCaloTowers"),
                            hbheInput = cms.InputTag("hbhereco"),
                            hoInput = cms.InputTag("horeco"),
                            hfInput = cms.InputTag("hfreco"),
                            ecalInputs = cms.VInputTag(cms.InputTag("reducedEcalRecHitsEB"),
                                                       cms.InputTag("reducedEcalRecHitsEE")),
                            hcalNoiseSummary = cms.InputTag("hcalnoise"),
                            hepMcCollection = cms.InputTag("source"),
                            genInfoCollection = cms.InputTag("source"),
                            genWeightCollection = cms.untracked.string('CSA07WeightProducer'),
                            nameFile = cms.untracked.string('analysisTree.root'),
                            nameTree = cms.untracked.string('ntp1'),
                            # switch ON/OFF the candidate collections to dump
                            dumpRunInfo = cms.untracked.bool(True),
                            dumpElectrons = cms.untracked.bool(True),
                            dumpPFlowElectrons = cms.untracked.bool(False),
                            dumpPFpreId = cms.untracked.bool(False),
                            dumpMuons = cms.untracked.bool(True),
                            dumpPFTaus = cms.untracked.bool(True),
                            dumpTracks = cms.untracked.bool(False),
                            dumpGsfTracks = cms.untracked.bool(True),
                            dumpMuonTracks = cms.untracked.bool(True),
                            dumpVertices = cms.untracked.bool(False),
                            dumpK0s = cms.untracked.bool(False),
                            dumpCaloTowers = cms.untracked.bool(False),
                            dumpSCs = cms.untracked.bool(False),
                            dumpBCs = cms.untracked.bool(False),
                            dumpJets = cms.untracked.bool(True),
                            dumpPUcorrPFJet = cms.untracked.bool(True),
                            dumpMet = cms.untracked.bool(True),
                            dumpHcalNoiseFlags = cms.untracked.bool(False),
                            # switch ON/OFF the particle flow objects to dump
                            dumpParticleFlowObjects = cms.untracked.bool(False),
                            # switch ON/OFF the additional informations to dump
                            saveEcal = cms.untracked.bool(True),
                            saveFatEcal = cms.untracked.bool(True),
                            saveTrk = cms.untracked.bool(True),
                            saveFatTrk = cms.untracked.bool(True),
                            saveTrackDeDx = cms.untracked.bool(False),
                            saveEleID = cms.untracked.bool(True),
                            savePFEleGsfTrk = cms.untracked.bool(True),
                            savePFEleBasic = cms.untracked.bool(True),
                            saveJetBTag = cms.untracked.bool(True),
                            savePFTauBasic = cms.untracked.bool(True),
                            saveLeadPFCand = cms.untracked.bool(True),
                            savePFTauDiscriminators = cms.untracked.bool(True),
                            # MC truth
                            mcTruthCollection = cms.InputTag("genParticles"),
                            electronMatchMap = cms.InputTag("electronMatchMap"),
                            muonMatchMap = cms.InputTag("muonMatchMap"),
                            dumpGenMet = cms.untracked.bool(True),
                            dumpGenJets = cms.untracked.bool(True),
                            dumpMCTruth = cms.untracked.bool(True),
                            dumpGenInfo = cms.untracked.bool(True),
                            dumpPreselInfo = cms.untracked.bool(False),
                            dumpSignalKfactor = cms.untracked.bool(True),
                            doMCEleMatch = cms.untracked.bool(False),
                            doMCMuonMatch = cms.untracked.bool(False),
                            # trigger results
                            dumpTriggerResults = cms.untracked.bool(False),
                            dumpHLTObjects = cms.untracked.bool(False),
                            #tau discriminators
                            tauDiscrByLeadTrackFindingTag = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"),
                            tauDiscrByLeadTrackPtCutTag = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"),
                            tauDiscrByTrackIsoTag = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"),
                            tauDiscrByEcalIsoTag = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation"),
                            tauDiscrAgainstMuonsTag = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"),
                            tauDiscrAgainstElectronsTag = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"),
                            tauDiscrByTaNCTag = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
                            tauDiscrByTaNCfrHalfPercentTag = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercent"),
                            tauDiscrByTaNCfrOnePercentTag = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercent"),
                            tauDiscrByTaNCfrQuarterPercentTag = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent"),
                            tauDiscrByTaNCfrTenthPercentTag = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercent"),
                            # effectively dump the data into the tree
                            dumpTree = cms.untracked.bool(False),
                            PFJetsBTags = cms.untracked.PSet( combinedSecondaryVertexBJetTags = cms.InputTag("newCombinedSecondaryVertexBPFJetTags"),
                                                              combinedSecondaryVertexMVABJetTags = cms.InputTag("newCombinedSecondaryVertexMVABPFJetTags"),
                                                              jetBProbabilityBJetTags = cms.InputTag("newJetBProbabilityBPFJetTags"),
                                                              jetProbabilityBJetTags = cms.InputTag("newJetProbabilityBPFJetTags"),
                                                              simpleSecondaryVertexHighEffBJetTags = cms.InputTag("newSimpleSecondaryVertexHighEffBPFJetTags"),
                                                              simpleSecondaryVertexHighPurBJetTags = cms.InputTag("newSimpleSecondaryVertexHighPurBPFJetTags"),
                                                              softMuonBJetTags = cms.InputTag("newSoftMuonBPFJetTags"),
                                                              softMuonByIP3dBJetTags = cms.InputTag("newSoftMuonByIP3dBPFJetTags"),
                                                              softMuonByPtBJetTags = cms.InputTag("newSoftMuonByPtBPFJetTags"),
                                                              softElectronBJetTags = cms.InputTag("newSoftElectronBPFJetTags"),
                                                              softElectronByIP3dBJetTags = cms.InputTag("newSoftElectronByIP3dBPFJetTags"),
                                                              softElectronByPtBJetTags = cms.InputTag("newSoftElectronByPtBPFJetTags"),
                                                              trackCountingHighPurBJetTags = cms.InputTag("newTrackCountingHighPurBPFJetTags"),
                                                              trackCountingHighEffBJetTags = cms.InputTag("newTrackCountingHighEffBPFJetTags")),
                            PFPUcorrJetsBTags = cms.untracked.PSet( combinedSecondaryVertexBJetTags = cms.InputTag("newCombinedSecondaryVertexBPFPUcorrJetTags"),
                                                                    combinedSecondaryVertexMVABJetTags = cms.InputTag("newCombinedSecondaryVertexMVABPFPUcorrJetTags"),
                                                                    jetBProbabilityBJetTags = cms.InputTag("newJetBProbabilityBPFPUcorrJetTags"),
                                                                    jetProbabilityBJetTags = cms.InputTag("newJetProbabilityBPFPUcorrJetTags"),
                                                                    simpleSecondaryVertexHighEffBJetTags = cms.InputTag("newSimpleSecondaryVertexHighEffBPFPUcorrJetTags"),
                                                                    simpleSecondaryVertexHighPurBJetTags = cms.InputTag("newSimpleSecondaryVertexHighPurBPFPUcorrJetTags"),
                                                                    softMuonBJetTags = cms.InputTag("newSoftMuonBPFPUcorrJetTags"),
                                                                    softMuonByIP3dBJetTags = cms.InputTag("newSoftMuonByIP3dBPFPUcorrJetTags"),
                                                                    softMuonByPtBJetTags = cms.InputTag("newSoftMuonByPtBPFPUcorrJetTags"),
                                                                    softElectronBJetTags = cms.InputTag("newSoftElectronBPFPUcorrJetTags"),
                                                                    softElectronByIP3dBJetTags = cms.InputTag("newSoftElectronByIP3dBPFPUcorrJetTags"),
                                                                    softElectronByPtBJetTags = cms.InputTag("newSoftElectronByPtBPFPUcorrJetTags"),
                                                                    trackCountingHighPurBJetTags = cms.InputTag("newTrackCountingHighPurBPFPUcorrJetTags"),
                                                                    trackCountingHighEffBJetTags = cms.InputTag("newTrackCountingHighEffBPFPUcorrJetTags"))
                            )
