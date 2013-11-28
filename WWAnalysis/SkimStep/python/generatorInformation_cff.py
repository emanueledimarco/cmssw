import FWCore.ParameterSet.Config as cms

def addGeneratorInfo(process,seq):
    process.prunedGen = cms.EDProducer( "GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "drop  *  ",
            "keep++ pdgId =   {Z0}",
            "++keep pdgId =   {Z0}",
            "keep++ pdgId =   {W+}",
            "++keep pdgId =   {W+}",
            "keep++ pdgId =   {W-}",
            "++keep pdgId =   {W-}",
            "keep++ pdgId =   {h0}",
            "++keep pdgId =   {h0}",
            "keep++ pdgId =   {e+}",
            "++keep pdgId =   {e+}",
            "keep++ pdgId =   {e-}",
            "++keep pdgId =   {e-}",
            "keep++ pdgId =  {mu+}",
            "++keep pdgId =  {mu+}",
            "keep++ pdgId =  {mu-}",
            "++keep pdgId =  {mu-}",
            "++keep pdgId = {tau+}",
            "keep++ pdgId = {tau+}",
            "++keep pdgId = {tau-}",
            "keep++ pdgId = {tau-}",
            "keep++ pdgId =     22",
            "++keep pdgId =     22",
            "++keep pdgId =      5",
            "++keep pdgId =     -5",
            "++keep pdgId =      4",
            "++keep pdgId =     -4",
            "++keep pdgId =     12",
            "++keep pdgId =    -12",
            "++keep pdgId =     14",
            "++keep pdgId =    -14",
            "++keep pdgId =    -16",
            "++keep pdgId =     16",
        )
    )
    
    seq += (process.prunedGen)
