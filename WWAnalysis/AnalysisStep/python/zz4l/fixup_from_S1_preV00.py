import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *

## Re-run Muon BDT Iso (tag has changed w.r.t. skim production)
boostedMuonsUpdatedBDTIso = cms.EDProducer("PatMuonBoosterBDTIso", 
      useExistingIsoValues = cms.bool(True),
      src = cms.InputTag("boostedMuons"),
      rho = cms.string("rhoMu"), 
      outputName = cms.string("bdtisoNew"),
      effectiveAreaTarget = cms.string("Data2011"),
)

## Re-run POG PF isolation for neutrals with threshold at pt > 0.5 GeV
updatedMuonPFIsoNHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuonsUpdatedBDTIso"),
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 130 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
updatedMuonPFIsoPhoton04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedMuonsUpdatedBDTIso"),
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 22 && pt > 0.5"), # neutral hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.0), # no special veto in the endcaps
    directional = cms.bool(False),
)
boostedMuonsUpdatedPFIso = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsUpdatedBDTIso"),
    valueMaps = cms.PSet(
        muonPFIsoPhoton04pt05_step1 = cms.InputTag("updatedMuonPFIsoPhoton04"),
        muonPFIsoNHad04pt05_step1   = cms.InputTag("updatedMuonPFIsoNHad04"),
    ),
)

fixupMuonSequence = cms.Sequence(
    boostedMuonsUpdatedBDTIso *
    updatedMuonPFIsoNHad04 * updatedMuonPFIsoPhoton04 * boostedMuonsUpdatedPFIso
)

#===========================================================================

## Recompute iso without vetos that are not in POG code
updatedElectronPFIsoChHad04 = cms.EDProducer("LeptonPFIsoFromStep1",
    leptonLabel = cms.InputTag("boostedElectrons"),
    endcapDefinition = cms.string("abs(superCluster.eta) > 1.479"), # for electrons we use supercluster eta, not the default momentum eta
    pfLabel     = cms.InputTag("reducedPFNoPUCands"), 
    pfSelection = cms.string("charge != 0 && abs(pdgId) == 211"), # charged hadrons
    deltaR     = cms.double(0.4), # radius
    deltaRself = cms.double(0.0), # no self veto
    vetoConeEndcaps = cms.double(0.015), 
    directional = cms.bool(False),
)
updatedElectronPFIsoNHad04 = updatedElectronPFIsoChHad04.clone(
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 130"), # neutral hadrons
    vetoConeEndcaps = cms.double(0.0), # No veto for NH, EGamma has abandoned it
)
updatedElectronPFIsoPhoton04 = updatedElectronPFIsoChHad04.clone(
    pfSelection = cms.string("charge == 0 && abs(pdgId) == 22"), # neutral hadrons
    vetoConeEndcaps = cms.double(0.08), # Preserve veto for Photons
)
updatedElectronPFIsoChHad01 = updatedElectronPFIsoChHad04.clone(deltaR = 0.1)
updatedElectronPFIsoChHad02 = updatedElectronPFIsoChHad04.clone(deltaR = 0.2)
updatedElectronPFIsoChHad03 = updatedElectronPFIsoChHad04.clone(deltaR = 0.3)
updatedElectronPFIsoChHad05 = updatedElectronPFIsoChHad04.clone(deltaR = 0.5)
updatedElectronPFIsoNHad01 = updatedElectronPFIsoNHad04.clone(deltaR = 0.1)
updatedElectronPFIsoNHad02 = updatedElectronPFIsoNHad04.clone(deltaR = 0.2)
updatedElectronPFIsoNHad03 = updatedElectronPFIsoNHad04.clone(deltaR = 0.3)
updatedElectronPFIsoNHad05 = updatedElectronPFIsoNHad04.clone(deltaR = 0.5)
updatedElectronPFIsoPhoton01 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.1)
updatedElectronPFIsoPhoton02 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.2)
updatedElectronPFIsoPhoton03 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.3)
updatedElectronPFIsoPhoton05 = updatedElectronPFIsoPhoton04.clone(deltaR = 0.5)

boostedElectronsUpdatedPFIso = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectrons"),
    valueMaps = cms.PSet(
        electronPFIsoChHad01_v2 = cms.InputTag("updatedElectronPFIsoChHad01"),
        electronPFIsoChHad02_v2 = cms.InputTag("updatedElectronPFIsoChHad02"),
        electronPFIsoChHad03_v2 = cms.InputTag("updatedElectronPFIsoChHad03"),
        electronPFIsoChHad04_v2 = cms.InputTag("updatedElectronPFIsoChHad04"),
        electronPFIsoChHad05_v2 = cms.InputTag("updatedElectronPFIsoChHad05"),
        electronPFIsoNHad01_v2 = cms.InputTag("updatedElectronPFIsoNHad01"),
        electronPFIsoNHad02_v2 = cms.InputTag("updatedElectronPFIsoNHad02"),
        electronPFIsoNHad03_v2 = cms.InputTag("updatedElectronPFIsoNHad03"),
        electronPFIsoNHad04_v2 = cms.InputTag("updatedElectronPFIsoNHad04"),
        electronPFIsoNHad05_v2 = cms.InputTag("updatedElectronPFIsoNHad05"),
        electronPFIsoPhoton01_v2 = cms.InputTag("updatedElectronPFIsoPhoton01"),
        electronPFIsoPhoton02_v2 = cms.InputTag("updatedElectronPFIsoPhoton02"),
        electronPFIsoPhoton03_v2 = cms.InputTag("updatedElectronPFIsoPhoton03"),
        electronPFIsoPhoton04_v2 = cms.InputTag("updatedElectronPFIsoPhoton04"),
        electronPFIsoPhoton05_v2 = cms.InputTag("updatedElectronPFIsoPhoton05"),
    )
)

## Fix BDT ID bug (userInt vs userFloat)
boostedElectronsFixBDTID = cms.EDProducer("PatElectronBoosterBDTID", 
    src = cms.InputTag("boostedElectronsUpdatedPFIso"),
    postfix = cms.string("_Fix"),
)
boostedElectronsFixBDTIso = cms.EDProducer("PatElectronBoosterBDTIso", 
    src = cms.InputTag("boostedElectronsFixBDTID"),
    rho = cms.string("rhoEl"),
    effectiveAreaTarget = cms.string("Data2011"),
    chargedOption = cms.string("_v2"), # postfix to the userFloat values
    neutralsOption = cms.string("_v2"), # postfix to the userFloat value
    outputName = cms.string("bdtisoFix"),
)

fixupElectronSequence = cms.Sequence(
    ( updatedElectronPFIsoPhoton01 + updatedElectronPFIsoPhoton02 + updatedElectronPFIsoPhoton03 + updatedElectronPFIsoPhoton04 + updatedElectronPFIsoPhoton05 +
      updatedElectronPFIsoNHad01 + updatedElectronPFIsoNHad02 + updatedElectronPFIsoNHad03 + updatedElectronPFIsoNHad04 + updatedElectronPFIsoNHad05 +
      updatedElectronPFIsoChHad01 + updatedElectronPFIsoChHad02 + updatedElectronPFIsoChHad03 + updatedElectronPFIsoChHad04 + updatedElectronPFIsoChHad05) *
    boostedElectronsUpdatedPFIso * 
    boostedElectronsFixBDTID * 
    boostedElectronsFixBDTIso  
)
 
##================================================================================================0

## Compute EA-corrected isolations
boostedElectronsEAPFIso.src = "boostedElectronsFixBDTIso"
boostedElectronsEAPFIso.chargedOption  = "_v2" # postfix to the userFloat values
boostedElectronsEAPFIso.neutralsOption = "_v2"
boostedElectronsEAPFIso.rho = 'rhoEl'

boostedMuonsEAPFIso.src = "boostedMuonsUpdatedPFIso"
boostedMuonsEAPFIso.neutralsOption = "pt05_step1"
boostedMuonsEAPFIso.rho = 'rhoMu'

## Embed final values and discriminators. Beyond here, fixups should be invisible
boostedElectrons.variables.pfChHadRelIso04 = "userFloat('electronPFIsoChHad04_v2')/pt"
boostedElectrons.variables.pfChHadIso04    = "userFloat('electronPFIsoChHad04_v2')"
boostedElectrons.variables.pfNHadIso04     = "userFloat('electronPFIsoNHad04_v2')"
boostedElectrons.variables.pfPhotonIso04   = "userFloat('electronPFIsoPhoton04_v2')"
boostedElectrons.variables.bdtID      = "userFloat('%s_Fix')" % EL_BDT
boostedElectrons.variables.bdtIDTrig  = "userFloat('bdttrig_Fix')"
boostedElectrons.variables.bdtIso     = "userFloat('bdtisoFix')"

boostedMuons.variables.pfNHadIso04   = "userFloat('muonPFIsoNHad04pt05_step1')"
boostedMuons.variables.pfPhotonIso04 = "userFloat('muonPFIsoPhoton04pt05_step1')"
boostedMuons.variables.bdtIso        = "userFloat('bdtisoNew')"

reboosting = cms.Sequence(
    fixupElectronSequence     *
    boostedElectronsEAPFIso   *  boostedElectrons +
    fixupMuonSequence   *
    boostedMuonsEAPFIso * boostedMuons 
)
