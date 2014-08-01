import FWCore.ParameterSet.Config as cms
caloclusterDumper = cms.EDAnalyzer("CaloClusterDumper",
                                   EBCaloClusters1 = cms.InputTag("fixed5x5ClustersReco","Barrel5x5Clusters"),
                                   EECaloClusters1 = cms.InputTag("fixed5x5ClustersReco","Endcap5x5Clusters"),

                                   EBCaloClusters2 = cms.InputTag("fixed5x5ClustersFit","Barrel5x5Clusters"),
                                   EECaloClusters2 = cms.InputTag("fixed5x5ClustersFit","Endcap5x5Clusters"),

                                   EBCaloClusters3 = cms.InputTag("fixed5x5ClustersFitNoPU","Barrel5x5Clusters"),
                                   EECaloClusters3 = cms.InputTag("fixed5x5ClustersFitNoPU","Endcap5x5Clusters"),

                                   ecalBarrelRecHits = cms.InputTag(""),
                                   ecalEndcapRecHits = cms.InputTag(""),

                                   mcTruthCollection = cms.InputTag("genParticles"),
                                   isMC = cms.untracked.bool(True),
                                   fileName = cms.untracked.string('tree.root'),
                                   nameTree = cms.untracked.string('ntp1')
                                   )


