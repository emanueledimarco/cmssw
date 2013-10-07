import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.Skimming.higgsToWW2Leptons_HLTPaths_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWMetCorrector_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWMuonIsolator_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWMuonIsolatorRef_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWMuonSelector_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWMuonSelectorRef_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWElectronIsolator_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWElectronIsolatorRef_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWElectronSelector_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWElectronSelectorRef_cfi import *
from HiggsAnalysis.HiggsToWW2Leptons.HWWPreselectionMarker_cfi import *

higgsToWW2LeptonsPreselectionSequence = cms.Sequence(
    higgsToWW2LeptonsHLTFilter*
    #muonCorrectedMET *
    selectedMuons *
    #selectedMuonsRef *
    isolatedMuons *
    #isolatedMuonsRef *
    selectedElectrons *
    #selectedElectronsRef *
    isolatedElectrons *
    #isolatedElectronsRef *
    preselectionMarker
    )
    
