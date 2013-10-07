import FWCore.ParameterSet.Config as cms

def addRhoVariables(process,seq,eleTag='gsfElectrons',muTag='muons',pfNoPUTag='pfNoPileUp'):

    process.kt6PFJets.doRhoFastjet = cms.bool(True)
    process.kt6PFJets.doAreaFastjet = cms.bool(True)
    process.kt6PFJets.Rho_EtaMax = cms.double(4.4)
    process.kt6PFJetsForIso = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(2.5) )
    process.kt6PFJetsNoPU = process.kt6PFJets.clone( src = pfNoPUTag )
    process.kt6PFJetsForIsoNoPU = process.kt6PFJetsForIso.clone( src = pfNoPUTag )

    # Re-cluster jets w/ pfNoPileUp
    process.ak5PFJetsNoPU = process.ak5PFJets.clone( src = pfNoPUTag )
    
    process.load("WWAnalysis.Tools.rhoValueMapProducer_cfi")
    
    process.rhoMu = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIso","rho",process.name_()), leptonTag = muTag)
    process.rhoEl = process.rhoMu.clone(leptonTag = eleTag)
    process.rhoMuNoPU = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIsoNoPU","rho",process.name_()), leptonTag = muTag)
    process.rhoElNoPU = process.rhoMuNoPU.clone(leptonTag = eleTag)
    
    process.valueMaps = cms.Sequence(
        process.rhoMu +
        process.rhoEl +
        process.rhoMuNoPU +
        process.rhoElNoPU 
    )

    seq += (
        process.ak5PFJets + 
        process.kt6PFJets +
        process.kt6PFJetsForIso + 
        process.pfPileUp *
        process.pfNoPileUp * (
            process.ak5PFJetsNoPU +
            process.kt6PFJetsNoPU + 
            process.kt6PFJetsForIsoNoPU ) *
        process.valueMaps 
    )

