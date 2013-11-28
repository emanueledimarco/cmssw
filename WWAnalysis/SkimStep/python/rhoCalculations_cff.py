import FWCore.ParameterSet.Config as cms

def addRhoVariables(process,seq,eleTag='gsfElectrons',muTag='muons'):
    process.load('RecoJets.Configuration.RecoPFJets_cff')
    process.kt6PFJetsForIso = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(2.5) )
    process.load("WWAnalysis.Tools.rhoValueMapProducer_cfi")  
    process.rhoMu = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIso","rho",process.name_()), leptonTag = muTag)
    process.rhoEl = process.rhoMu.clone(leptonTag = eleTag)

    process.valueMaps = cms.Sequence(
        process.rhoMu +
        process.rhoEl
    )

    seq += (
        process.ak5PFJets + 
        process.kt6PFJets +
        process.kt6PFJetsForIso +
        process.valueMaps 
    )
        
