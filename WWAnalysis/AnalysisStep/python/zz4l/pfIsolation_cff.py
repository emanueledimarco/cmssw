import FWCore.ParameterSet.Config as cms

isoLLLL = cms.EDProducer("CompositeCandidateIsoAdder",
    src = cms.InputTag("recLLLL"),
    deposits = cms.VPSet(
    cms.PSet(# Charged (excluding photon conversions)
        label = cms.string("pfiso4lCharged"),
        muonIsolation = cms.PSet(
            src    = cms.InputTag("isoDepMuFromPFCharged"),
            deltaR = cms.double(0.3),
            muonVetos     = cms.vstring("ConeVeto(0.01)"),
            electronVetos = cms.vstring("ConeVeto(0.015)"),
            photonVetos   = cms.vstring(), # note: we don't have FSR photons yet
        ),
        electronIsolation = cms.PSet(
            src    = cms.InputTag("isoDepElFromPFCharged"),
            deltaR = cms.double(0.4),
            muonVetos     = cms.vstring("ConeVeto(0.01)"),
            electronVetos = cms.vstring("ConeVeto(0.015)"),
            photonVetos   = cms.vstring(), # note: we don't have FSR photons yet
        ),
    ), cms.PSet( # Neutral hadrons
        label = cms.string("pfiso4lNeutral"),
        muonIsolation = cms.PSet(
            src    = cms.InputTag("isoDepMuFromPFNeutral"),
            deltaR = cms.double(0.3),
            muonVetos     = cms.vstring("Threshold(0.5)"),
            electronVetos = cms.vstring("Threshold(0.5)", "ConeVeto(0.07)"), # 0.07 veto cone cone
            photonVetos   = cms.vstring(), # note: we don't have FSR photons yet
        ),
        electronIsolation = cms.PSet(
            src    = cms.InputTag("isoDepElFromPFNeutral"),
            deltaR = cms.double(0.4),
            muonVetos     = cms.vstring("Threshold(0.5)"),
            electronVetos = cms.vstring("Threshold(0.5)", "ConeVeto(0.07)"), # 0.07 veto cone cone
            photonVetos   = cms.vstring(), # note: we don't have FSR photons yet
        ),
    ), cms.PSet( # Photons       
        label = cms.string("pfiso4lPhotons"),
        muonIsolation = cms.PSet(
            src    = cms.InputTag("isoDepMuFromPFPhotons"),
            deltaR = cms.double(0.3),
            muonVetos     = cms.vstring("Threshold(0.5)"),
            electronVetos = cms.vstring("Threshold(0.5)", "ConeVeto(0.08)", "RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)"),
            photonVetos   = cms.vstring(), # note: we don't have FSR photons yet
        ),
        electronIsolation = cms.PSet(
            src    = cms.InputTag("isoDepElFromPFPhotons"),
            deltaR = cms.double(0.4),
            muonVetos     = cms.vstring("Threshold(0.5)"),
            electronVetos = cms.vstring("Threshold(0.5)", "ConeVeto(0.08)", "RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)"),
            photonVetos   = cms.vstring(), # note: we don't have FSR photons yet
        ),
    ))
)


