import FWCore.ParameterSet.Config as cms
superclusterDumper = cms.EDAnalyzer("SuperClusterDumper",
                                    EBSuperClusters1 = cms.InputTag("correctedHybridSuperClusters","","RECO"),
                                    EESuperClusters1 = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","","RECO"),
                                    EBSuperClusters2 = cms.InputTag("correctedHybridSuperClusters","","HitsToClusters"),
                                    EESuperClusters2 = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","","HitsToClusters"),
                                    EBSuperClusters3 = cms.InputTag("correctedHybridSuperClustersNoPU","","HitsToClusters"),
                                    EESuperClusters3 = cms.InputTag("correctedMulti5x5SuperClustersWithPreshowerNoPU","","HitsToClusters"),
                                    ecalBarrelRecHits = cms.InputTag(""),
                                    ecalEndcapRecHits = cms.InputTag(""),
                                    mcTruthCollection = cms.InputTag("genParticles"),
                                    isMC = cms.untracked.bool(True),
                                    fileName = cms.untracked.string('tree.root'),
                                    nameTree = cms.untracked.string('ntp1')
                                    )


