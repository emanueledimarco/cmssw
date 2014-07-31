import FWCore.ParameterSet.Config as cms
superclusterDumper = cms.EDAnalyzer("SuperClusterDumper",
                                    EBSuperClusters1 = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel","RECO"),
                                    EESuperClusters1 = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower","RECO"),
                                    EBSuperClusters2 = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel","HitsToClusters"),
                                    EESuperClusters2 = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower","HitsToClusters"),
                                    EBSuperClusters3 = cms.InputTag("particleFlowSuperClusterECALNoPU:particleFlowSuperClusterECALBarrel","HitsToClusters"),
                                    EESuperClusters3 = cms.InputTag("particleFlowSuperClusterECALNoPU:particleFlowSuperClusterECALEndcapWithPreshower","HitsToClusters"),
                                    ecalBarrelRecHits = cms.InputTag(""),
                                    ecalEndcapRecHits = cms.InputTag(""),
                                    mcTruthCollection = cms.InputTag("genParticles"),
                                    isMC = cms.untracked.bool(True),
                                    fileName = cms.untracked.string('tree.root'),
                                    nameTree = cms.untracked.string('ntp1')
                                    )


