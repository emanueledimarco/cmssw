import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets

#needed for MLP tagger
goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
)

#needed for MLP tagger
kt6PFJets = kt6PFJets.clone( rParam = 0.6, doRhoFastjet = True )

#needed for Likelihood tagger
from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJetsIso = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
kt6PFJetsIso.Rho_EtaMax = cms.double(2.5)


QGTagger = cms.EDProducer('QGTagger',
  srcJets         = cms.InputTag('ak5PFJetsCorr'),
  srcRho          = cms.InputTag('kt6PFJets','rho'),
  srcRhoIso       = cms.InputTag('kt6PFJetsIso','rho'),
)

QuarkGluonTagger = cms.Sequence(goodOfflinePrimaryVertices+kt6PFJets+kt6PFJetsIso+QGTagger)
