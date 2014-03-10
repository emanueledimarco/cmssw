import FWCore.ParameterSet.Config as cms
from CMGTools.External.pujetidsequence_cff import puJetMva

treeDumper = cms.EDAnalyzer("HWWTreeDumper",
                            HLTObjectsInfo = cms.untracked.PSet(triggerResults = cms.InputTag("TriggerResults","","AUTO"),
                                                                processName = cms.string("AUTO"),
                                                                triggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD","","HLT")
                                                                ),
                            electronCollection = cms.InputTag("gsfElectrons"),
                            calibElectronCollection = cms.InputTag("calibratedElectrons","calibratedGsfElectrons"),
                            #pflowElectronCollection = cms.InputTag("particleFlow","electrons"),
                            pflowElectronCollection = cms.InputTag("particleFlow"),
                            photonCollection = cms.InputTag("goodPhotons"),
                            muonCollection = cms.InputTag("muons"),
                            calibMuonCollection = cms.InputTag("scaledMuons"),
                            PFCandidateCollection = cms.InputTag("particleFlow"),
                            PFNoPUCandidateCollection = cms.InputTag("pfNoPileUp"),
                            PFPUCandidateCollection = cms.InputTag("pfPileUp"),
                            ecalSCCollection = cms.InputTag("electronAndPhotonSuperClusters"),
                            ecalBarrelSCCollection = cms.InputTag("correctedHybridSuperClusters"),
                            ecalEndcapSCCollection = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapSuperClusters"),
                            ecalElePFClusterCollection = cms.InputTag("pfElectronTranslator","pf"),
                            ecalPhoPFClusterCollection = cms.InputTag("pfPhotonTranslator","pfphot"),
                            ecalBCCollection = cms.InputTag("seedBasicClusters"),
                            ecalBarrelRecHits = cms.InputTag("reducedEcalRecHitsEB"),
                            ecalEndcapRecHits = cms.InputTag("reducedEcalRecHitsEE"),
                            esRecHits = cms.InputTag("reducedEcalRecHitsES"),
                            calotowersForIsolationProducer = cms.InputTag("towerMaker"),
                            generalTrackCollection = cms.InputTag("generalTracks"),
                            trackCollection = cms.InputTag("leptonLinkedTracks"),
                            gsfTrackCollection = cms.InputTag("electronGsfTracks"),
                            globalMuonTrackCollection = cms.InputTag("globalMuons"),
                            standAloneMuonTrackCollection = cms.InputTag("standAloneMuons"),
                            #  refittedForDeDxTrackCollection = cms.InputTag("RefitterForDeDx"),
                            refittedForDeDxTrackCollection = cms.InputTag("generalTracks"),
                            vertexCollection = cms.InputTag("offlinePrimaryVertices"),
                            K0sCollection = cms.InputTag("generalV0Candidates","Kshort"),
                            genJetCollection = cms.InputTag("goodGenJets"),
                            jetCollection1 = cms.InputTag("ak5CaloJets"),
                            jetCollection2 = cms.InputTag("ak5CaloJets","","RECO"), # the process name is needed because that collection only has the jet ID attached
                            JPTjetCollection1 = cms.InputTag("ak5JPTJetsL2L3Residual"),
                            JPTjetCollection2 = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
                            PFjetCollection1 = cms.InputTag("goodPFCHSJets"),
                            PFJetCorrectionService = cms.string("ak5PFL1FastL2L3Residual"),
                            JetCorrectionService = cms.string("ak5CaloL1FastL2L3Residual"),
                            PFpuCorrJetCollection1 = cms.InputTag("goodPFJets"),
                            # jet id MVA
                            puJetIDAlgos = puJetMva.algos,
                            metCollection = cms.InputTag("met"), # preselection
                            # corrmetCollection = cms.InputTag("metMuonJESCorAK5"), # type I and II corr applied
                            TCmetCollection = cms.InputTag("tcMet"),
                            genMetCollection = cms.InputTag("genMetCalo"),
                            PFmetCollection = cms.InputTag("pfmets"),
                            chargedMetCollection = cms.InputTag("chargedMetProducer"), # std, one per each vertex
                            PFChMetCollection = cms.InputTag("ourChPFMet"),
                            leptonLinkedPFCandidates = cms.InputTag("reducedPFCandsToSave"),
                            PFpreIdCollection = cms.InputTag("trackerDrivenElectronSeeds:preid"),
                            calotowerCollection = cms.InputTag("lowThrCaloTowers"),
                            EBRecHits = cms.InputTag("reducedEcalRecHitsEB"),
                            EERecHits = cms.InputTag("reducedEcalRecHitsEE"),
                            conversionCollection = cms.InputTag("allConversions"),
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
                            # this is off, because we now store the calibrated energy in the std ele collection
                            dumpCalibratedElectrons = cms.untracked.bool(False),
                            dumpPFlowElectrons = cms.untracked.bool(False),
                            dumpPFpreId = cms.untracked.bool(False),
                            dumpPhotons = cms.untracked.bool(True),
                            dumpConversions = cms.untracked.bool(True),
                            dumpMuons = cms.untracked.bool(True),
                            dumpPFCandidates = cms.untracked.bool(True),
                            dumpTracks = cms.untracked.bool(False),
                            dumpLinkedTracks = cms.untracked.bool(True),
                            dumpGsfTracks = cms.untracked.bool(True),
                            dumpMuonTracks = cms.untracked.bool(True),
                            dumpVertices = cms.untracked.bool(False),
                            dumpK0s = cms.untracked.bool(False),
                            dumpCaloTowers = cms.untracked.bool(False),
                            dumpEcalRecHits = cms.untracked.bool(False),
                            dumpSCs = cms.untracked.bool(False),
                            dumpBCs = cms.untracked.bool(False),
                            dumpPhoPFClusters = cms.untracked.bool(False),
                            dumpJets = cms.untracked.bool(False),
                            dumpPUcorrPFJet = cms.untracked.bool(True),
                            dumpCHSPFJet = cms.untracked.bool(False),
                            dumpMet = cms.untracked.bool(True),
                            dumpLogErrorFlags = cms.untracked.bool(True),
                            dumpHcalNoiseFlags = cms.untracked.bool(False),
                            AODHcalNoiseFlags = cms.untracked.bool(True),
                            # switch ON/OFF the particle flow objects to dump
                            dumpParticleFlowObjects = cms.untracked.bool(False),
                            # switch ON/OFF the additional informations to dump
                            saveEcal = cms.untracked.bool(True),
                            saveFatEcal = cms.untracked.bool(True),
                            saveTrk = cms.untracked.bool(True),
                            saveFatTrk = cms.untracked.bool(True),
                            saveTrackDeDx = cms.untracked.bool(False),
                            saveTrackImpactParameters = cms.untracked.bool(True),
                            saveEleID = cms.untracked.bool(True),
                            savePFEleGsfTrk = cms.untracked.bool(True),
                            savePFEleBasic = cms.untracked.bool(True),
                            saveJetBTag = cms.untracked.bool(True),
                            saveLeadPFCand = cms.untracked.bool(True),
                            # MC truth
                            mcTruthCollection = cms.InputTag("genParticles"),
                            dumpGenMet = cms.untracked.bool(False),
                            dumpGenJets = cms.untracked.bool(False),
                            dumpMCTruth = cms.untracked.bool(False),
                            dumpMCTruthExtra = cms.untracked.bool(False), 
                            is8TeV = cms.untracked.bool(True),
                            dumpGenInfo = cms.untracked.bool(False),
                            dumpLHE = cms.untracked.bool(False),
                            dumpPreselInfo = cms.untracked.bool(False),
                            # trigger results
                            dumpTriggerResults = cms.untracked.bool(False),
                            dumpHLTObjects = cms.untracked.bool(False),
                            # PDFs
                            dumpPdfWeight = cms.untracked.bool(False),
                            # pdf
                            pdfSet1 = cms.InputTag("pdfWeights:cteq66"),
                            pdfSet2 = cms.InputTag("pdfWeights:MRST2006nnlo"),
                            pdfSet3 = cms.InputTag("pdfWeights:NNPDF10"), # no _ are allowed in the names: NN10_100 => NN10
                            namepdf1 = cms.untracked.string("CTEQ66"),
                            namepdf2 = cms.untracked.string("MRST2006NNLO"),
                            namepdf3 = cms.untracked.string("NNPDF10100"),

                            dumpTree = cms.untracked.bool(False),
                            PFJetsBTags = cms.untracked.PSet( combinedSecondaryVertexBJetTags = cms.InputTag("newCombinedSecondaryVertexBPFNoPUJetTags"),
                                                              combinedSecondaryVertexMVABJetTags = cms.InputTag("newCombinedSecondaryVertexMVABPFNoPUJetTags"),
                                                              jetBProbabilityBJetTags = cms.InputTag("newJetBProbabilityBPFNoPUJetTags"),
                                                              jetProbabilityBJetTags = cms.InputTag("newJetProbabilityBPFNoPUJetTags"),
                                                              simpleSecondaryVertexHighEffBJetTags = cms.InputTag("newSimpleSecondaryVertexHighEffBPFNoPUJetTags"),
                                                              simpleSecondaryVertexHighPurBJetTags = cms.InputTag("newSimpleSecondaryVertexHighPurBPFNoPUJetTags"),
                                                              trackCountingHighPurBJetTags = cms.InputTag("newTrackCountingHighPurBPFNoPUJetTags"),
                                                              trackCountingHighEffBJetTags = cms.InputTag("newTrackCountingHighEffBPFNoPUJetTags"),
                                                              trackCountingVeryHighEffBJetTags = cms.InputTag("newTrackCountingVeryHighEffBPFNoPUJetTags")),
                            PFPUcorrJetsBTags = cms.untracked.PSet( combinedSecondaryVertexBJetTags = cms.InputTag("newCombinedSecondaryVertexBPFPUcorrJetTags"),
                                                                    combinedSecondaryVertexMVABJetTags = cms.InputTag("newCombinedSecondaryVertexMVABPFPUcorrJetTags"),
                                                                    jetBProbabilityBJetTags = cms.InputTag("newJetBProbabilityBPFPUcorrJetTags"),
                                                                    jetProbabilityBJetTags = cms.InputTag("newJetProbabilityBPFPUcorrJetTags"),
                                                                    simpleSecondaryVertexHighEffBJetTags = cms.InputTag("newSimpleSecondaryVertexHighEffBPFPUcorrJetTags"),
                                                                    simpleSecondaryVertexHighPurBJetTags = cms.InputTag("newSimpleSecondaryVertexHighPurBPFPUcorrJetTags"),
                                                                    trackCountingHighPurBJetTags = cms.InputTag("newTrackCountingHighPurBPFPUcorrJetTags"),
                                                                    trackCountingHighEffBJetTags = cms.InputTag("newTrackCountingHighEffBPFPUcorrJetTags"),
                                                                    trackCountingVeryHighEffBJetTags = cms.InputTag("newTrackCountingVeryHighEffBPFPUcorrJetTags"))
                            )
