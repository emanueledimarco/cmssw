import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *

## Re-run POG PF isolation for neutrals with threshold at pt > 0.5 GeV and
## and vetos as per RecoMuon/MuonIsolation/python/muonPFIsolationValues_cff.py
updatedMuonPFIsoChHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuons"),
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge != 0 && abs(pdgId) == 211"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0001), # self-veto 0.0001
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
updatedMuonPFIsoNHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuons"),
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 130 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.01), # self/fsr veto 0.01
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
updatedMuonPFIsoPhoton04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuons"),
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 22 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.01), # self/fsr veto 0.01
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
boostedMuonsUpdatedPFIso = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuons"),
    valueMaps = cms.PSet(
        muonPFIsoChHad04_step1   = cms.InputTag("updatedMuonPFIsoChHad04"),
        muonPFIsoPhoton04pt05_step1 = cms.InputTag("updatedMuonPFIsoPhoton04"),
        muonPFIsoNHad04pt05_step1   = cms.InputTag("updatedMuonPFIsoNHad04"),
    ),
)

fixupMuonSequence = cms.Sequence(
    ( updatedMuonPFIsoChHad04  + updatedMuonPFIsoNHad04 + updatedMuonPFIsoPhoton04 ) * 
    boostedMuonsUpdatedPFIso
)

boostedMuonsEAPFIso.src = "boostedMuonsUpdatedPFIso"
boostedMuonsEAPFIso.neutralsOption = "pt05_step1"
boostedMuonsEAPFIso.chargedOption  = "_step1" # postfix to the userFloat values

boostedMuons.variables.pfChHadRelIso04 = "userFloat('muonPFIsoChHad04_step1')/pt"
boostedMuons.variables.pfChHadIso04    = "userFloat('muonPFIsoChHad04_step1')"
boostedMuons.variables.pfNHadIso04   = "userFloat('muonPFIsoNHad04pt05_step1')"
boostedMuons.variables.pfPhotonIso04 = "userFloat('muonPFIsoPhoton04pt05_step1')"

reboosting = cms.Sequence(
    boostedElectronsEAPFIso *  boostedElectrons +
    fixupMuonSequence   *
    boostedMuonsEAPFIso * boostedMuons 
)
