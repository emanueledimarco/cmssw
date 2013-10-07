import FWCore.ParameterSet.Config as cms

#     __  ___                      ____  _______           
#    /  |/  /__  ______  ____     / __ \/ ____(_)_________ 
#   / /|_/ // / / / __ \/ __ \   / /_/ / /_  / // ___/ __ \
#  / /  / // /_/ / /_/ / / / /  / ____/ __/ / /(__  ) /_/ /
# /_/  /_/ \__,_/\____/_/ /_/  /_/   /_/   /_//____/\____/ 
#


updatedMuonPFIsoNHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuons"),
    pfLabel     = cms.InputTag("reducedPFCandsPfNoPU"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 130 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
updatedMuonPFIsoPhoton04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuons"),
    pfLabel     = cms.InputTag("reducedPFCandsPfNoPU"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 22 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
boostedMuonsUpdatedPFIso = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuons"),
    valueMaps = cms.PSet(
        muonPFIsoPhoton04pt05_step1 = cms.InputTag("updatedMuonPFIsoPhoton04"),
        muonPFIsoNHad04pt05_step1   = cms.InputTag("updatedMuonPFIsoNHad04"),
    ),
)

#     ________              ____  _______           
#    / ____/ /___  _____   / __ \/ ____(_)_________ 
#   / __/ / // _ \/ ___/  / /_/ / /_  / // ___/ __ \
#  / /___/ //  __/ /__   / ____/ __/ / /(__  ) /_/ /
# /_____/_/ \___/\___/  /_/   /_/   /_//____/\____/ 
#


updatedElectronPFIsoChHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedElectrons"),
    pfLabel     = cms.InputTag("reducedPFCandsPfNoPU"), 
    pfSelection = cms.string("charge != 0 && abs(pdgId) == 211"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.015), 
    directional = cms.bool(False),
)
updatedElectronPFIsoChHad01 = updatedElectronPFIsoChHad04.clone(deltaR = 0.1)
updatedElectronPFIsoChHad02 = updatedElectronPFIsoChHad04.clone(deltaR = 0.2)
updatedElectronPFIsoChHad03 = updatedElectronPFIsoChHad04.clone(deltaR = 0.3)
updatedElectronPFIsoChHad05 = updatedElectronPFIsoChHad04.clone(deltaR = 0.5)
boostedElectronsUpdatedPFIso = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectrons"),
    valueMaps = cms.PSet( 
        electronPFIsoChHad01_v2 = cms.InputTag("updatedElectronPFIsoChHad01"),
        electronPFIsoChHad02_v2 = cms.InputTag("updatedElectronPFIsoChHad02"),
        electronPFIsoChHad03_v2 = cms.InputTag("updatedElectronPFIsoChHad03"),
        electronPFIsoChHad04_v2 = cms.InputTag("updatedElectronPFIsoChHad04"),
        electronPFIsoChHad05_v2 = cms.InputTag("updatedElectronPFIsoChHad05"),
    )
)

#     _________    ____  __________          
#    / ____/   |  / __ \/ ____/  _/_________ 
#   / __/ / /| | / /_/ / /_   / / / ___/ __ \
#  / /___/ ___ |/ ____/ __/ _/ / (__  ) /_/ /
# /_____/_/  |_/_/   /_/   /___//____/\____/ 
#                                            


boostedElectronsEAPFIso = cms.EDProducer("PatElectronEffAreaIso",
    src = cms.InputTag("boostedElectronsUpdatedPFIso"),
    rho = cms.string("rhoElFullEta"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    chargedOption = cms.string("_v2"), # postfix to the userFloat values
    effectiveAreaTarget = cms.string("Data2012"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
    #separatePhotonAndHadronEAs = cms.bool(True), # use separte EA for photons and hadrons
    #truncateAtZero = cms.string("both"), # (yes|no) for total EA, (both|sum|no) for separate EA
)
boostedMuonsEAPFIso = cms.EDProducer("PatMuonEffAreaIso",
    src = cms.InputTag("boostedMuonsUpdatedPFIso"),
    rho = cms.string("rhoMuFullEta"),
    deltaR = cms.string("04"),
    neutralsOption = cms.string("pt05_step1"), # postfix to the userFloat value
    label = cms.string("pfCombIso04EACorr"),
    effectiveAreaTarget = cms.string("Data2012"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
)
boostedElectrons = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectronsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('electronPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('electronPFIsoChHad04_v2')"),
        pfChHadIso04_old = cms.string("userFloat('electronPFIsoChHad04')"),
        sip   = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
)
boostedMuons = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfChHadRelIso04 = cms.string("userFloat('muonPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('muonPFIsoChHad04')"),
        bdtIso = cms.string("userFloat('bdtisoNew')"),
        sip    = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
)


# NEW 2012 JEC on the fly

slimPatJetsTriggerMatch = cms.EDProducer("PatJetJEC2012Adder",
   jetTag = cms.InputTag("slimPatJetsTriggerMatch"),
   rhoTag = cms.InputTag("kt6PFJets","rho"),
   isData=  cms.untracked.bool(True)
)

reboosting = cms.Sequence(
    ( updatedElectronPFIsoChHad01 + updatedElectronPFIsoChHad02 + updatedElectronPFIsoChHad03 + updatedElectronPFIsoChHad04 + updatedElectronPFIsoChHad05) *
    boostedElectronsUpdatedPFIso * 
    boostedElectronsEAPFIso   * boostedElectrons +
    (updatedMuonPFIsoNHad04 + updatedMuonPFIsoPhoton04) *
    boostedMuonsUpdatedPFIso *
    boostedMuonsEAPFIso * boostedMuons * slimPatJetsTriggerMatch
)

