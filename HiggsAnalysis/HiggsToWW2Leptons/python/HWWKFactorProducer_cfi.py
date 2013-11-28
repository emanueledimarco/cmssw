import FWCore.ParameterSet.Config as cms

KFactorProducer = cms.EDProducer("HWWKFactorProducer",
                                 genParticlesTag = cms.InputTag("genParticles"),
                                 inputFilename = cms.untracked.string('HiggsAnalysis/HiggsToWW2Leptons/data/160_10TeV.dat'),
                                 # Powheg Process to apply KFactor
                                 ProcessID = cms.untracked.int32(10010),
                                 Debug =cms.untracked.bool(False)
                                 )
