import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import *

from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi import *
from WWAnalysis.SkimStep.simpleCutBasedElectronIDSpring11_cfi import simpleCutBasedElectronID
from WWAnalysis.SkimStep.cutsInCategoriesHWWElectronIdentificationV04_cfi import *
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi import *
import RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi as cicHZZ

def addElectronIDs(process,seq):
    process.eIdSequence = cms.Sequence()
    
    process.vbtf11WP60 = simpleCutBasedElectronID.clone( electronQuality = '60' )
    process.vbtf11WP70 = simpleCutBasedElectronID.clone( electronQuality = '70' )
    process.vbtf11WP80 = simpleCutBasedElectronID.clone( electronQuality = '80' )
    process.vbtf11WP85 = simpleCutBasedElectronID.clone( electronQuality = '85' )
    process.vbtf11WP90 = simpleCutBasedElectronID.clone( electronQuality = '90' )
    process.vbtf11WP95 = simpleCutBasedElectronID.clone( electronQuality = '95' )
    process.eIdSequence += process.vbtf11WP60
    process.eIdSequence += process.vbtf11WP70
    process.eIdSequence += process.vbtf11WP80
    process.eIdSequence += process.vbtf11WP85
    process.eIdSequence += process.vbtf11WP90
    process.eIdSequence += process.vbtf11WP95
    
    process.cicVeryLooseHWW   = eidHWWVeryLoose.clone()
    process.cicLooseHWW       = eidHWWLoose.clone()
    process.cicMediumHWW      = eidHWWMedium.clone()
    process.cicTightHWW       = eidHWWTight.clone()
    process.cicSuperTightHWW  = eidHWWSuperTight.clone()
    process.cicHyperTight1HWW = eidHWWHyperTight1.clone()
    process.cicHyperTight2HWW = eidHWWHyperTight2.clone()
    process.cicHyperTight3HWW = eidHWWHyperTight3.clone()
    process.eIdSequence += process.cicVeryLooseHWW
    process.eIdSequence += process.cicLooseHWW
    process.eIdSequence += process.cicMediumHWW
    process.eIdSequence += process.cicTightHWW
    process.eIdSequence += process.cicSuperTightHWW 
    process.eIdSequence += process.cicHyperTight1HWW 
    process.eIdSequence += process.cicHyperTight2HWW 
    process.eIdSequence += process.cicHyperTight3HWW 
    
    process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
    process.egammaIDLikelihood = process.eidLikelihoodExt.clone()
    process.eIdSequence += process.egammaIDLikelihood 
    
    process.cicVeryLoose   = eidVeryLoose.clone()
    process.cicLoose       = eidLoose.clone()
    process.cicMedium      = eidMedium.clone()
    process.cicTight       = eidTight.clone()
    process.cicSuperTight  = eidSuperTight.clone()
    process.cicHyperTight1 = eidHyperTight1.clone()
    process.cicHyperTight2 = eidHyperTight2.clone()
    process.cicHyperTight3 = eidHyperTight3.clone()
    process.eIdSequence += process.cicVeryLoose 
    process.eIdSequence += process.cicLoose 
    process.eIdSequence += process.cicMedium 
    process.eIdSequence += process.cicTight 
    process.eIdSequence += process.cicSuperTight 
    process.eIdSequence += process.cicHyperTight1 
    process.eIdSequence += process.cicHyperTight2 
    process.eIdSequence += process.cicHyperTight3 
    
    process.cicVeryLooseMC   = eidVeryLooseMC.clone()
    process.cicLooseMC       = eidLooseMC.clone()
    process.cicMediumMC      = eidMediumMC.clone()
    process.cicTightMC       = eidTightMC.clone()
    process.cicSuperTightMC  = eidSuperTightMC.clone()
    process.cicHyperTight1MC = eidHyperTight1MC.clone()
    process.cicHyperTight2MC = eidHyperTight2MC.clone()
    process.cicHyperTight3MC = eidHyperTight3MC.clone()
    process.eIdSequence += process.cicVeryLooseMC 
    process.eIdSequence += process.cicLooseMC 
    process.eIdSequence += process.cicMediumMC 
    process.eIdSequence += process.cicTightMC 
    process.eIdSequence += process.cicSuperTightMC 
    process.eIdSequence += process.cicHyperTight1MC 
    process.eIdSequence += process.cicHyperTight2MC 
    process.eIdSequence += process.cicHyperTight3MC 

    for X in "VeryLoose", "Loose", "Medium", "Tight", "SuperTight", "HyperTight1": 
        setattr(process, 'cicHZZ'+X, getattr(cicHZZ, 'eidHZZ'+X).clone())
        process.eIdSequence += getattr(process, 'cicHZZ'+X)

    seq += process.eIdSequence
    return listModules(process.eIdSequence)
    
