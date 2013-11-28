import FWCore.ParameterSet.Config as cms

        
from PhysicsTools.PatAlgos.patSequences_cff import *

eleBDTSelection = cms.EDFilter("PatEleBDTSelection",
  src = cms.InputTag("boostedElectrons"),
  filter = cms.bool(False),
)


selectedRefPatElectrons = cms.EDFilter("PATElectronRefSelector",
  src = cms.InputTag("boostedElectrons"),
#    src = cms.InputTag("eleBDTSelection"),
   filter = cms.bool(False),
   cut = cms.string(""),
)


ELE_BASE = "( pt > 10 && abs(eta)<2.5 )"
wwEleMatch = selectedRefPatElectrons.clone( cut = ELE_BASE )

from WWAnalysis.AnalysisStep.electronIDs_cff import *

#                     
# Loose
#wwEleLoose      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LOOSE_2011 )
wwEleLoose      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LOOSE_SMURFS + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP )


# LHL
wwEleIDLHL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 )
wwEleISOLHL     = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 + " && " + ELE_ISO_LH_95_2011 ) 
wwEleCONVLHL    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 + " && " + ELE_ISO_LH_95_2011  + " && " + ELE_NOCONV ) 
wwEleIPLHL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_95_2011 + " && " + ELE_ISO_LH_95_2011  + " && " + ELE_NOCONV + " && " + ELE_IP ) 

# LHT
wwEleIDLHT      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 ) 
wwEleISOLHT     = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISO_LH_90_2011 ) 
wwEleCONVLHT    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISO_LH_90_2011 + " && " + ELE_NOCONV ) 
wwEleIPLHT      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISO_LH_90_2011 + " && " + ELE_NOCONV + " && " + ELE_IP ) 

wwEleISOPFLHT   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISOPF_LH_90_2011 ) 
wwEleCONVPFLHT  = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISOPF_LH_90_2011 + " && " + ELE_NOCONV ) 
wwEleIPPFLHT    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_LH_90_2011 + " && " + ELE_ISOPF_LH_90_2011 + " && " + ELE_NOCONV + " && " + ELE_IP )

# CBL
wwEleIDCBL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 ) 
wwEleISOCBL     = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 + " && " + ELE_ISO_CB_95_2011 ) 
wwEleCONVCBL    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 + " && " + ELE_ISO_CB_95_2011 + " && " + ELE_NOCONV ) 
wwEleIPCBL      = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_ID_CB_95_2011 + " && " + ELE_ISO_CB_95_2011 + " && " + ELE_NOCONV + " && " + ELE_IP )

# CBT
wwEleIDCBT      = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 ) 
wwEleISOCBT     = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 + " && " + ELE_ISO_CB_90_2011 ) 
wwEleCONVCBT    = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 + " && " + ELE_ISO_CB_90_2011 + " && " + ELE_NOCONV ) 
wwEleIPCBT      = selectedRefPatElectrons.clone( cut = ELE_BASE+ " && " + ELE_ID_CB_90_2011 + " && " + ELE_ISO_CB_90_2011 + " && " + ELE_NOCONV + " && " + ELE_IP )

# Merge
wwEleIDMerge    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID ) 
wwEleISOMerge   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO ) 
wwEleCONVMerge  = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV ) 
wwEleIPMerge    = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP )

# Merge2
wwEleIDMerge2   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 ) 
wwEleISOMerge2  = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 + " && " + ELE_MERGE_ISO ) 
wwEleCONVMerge2 = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV ) 
wwEleIPMerge2   = selectedRefPatElectrons.clone( cut = ELE_BASE + " && " + ELE_MERGE_ID2 + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP ) 

#       _                 _                         
#      | |               | |                        
#      | | __ _ _ __ ___ | |__   ___  _ __ ___  ___ 
#  _   | |/ _` | '_ ` _ \| '_ \ / _ \| '__/ _ \/ _ \
# | |__| | (_| | | | | | | |_) | (_) | | |  __/  __/
#  \____/ \__,_|_| |_| |_|_.__/ \___/|_|  \___|\___|
#                                                   
#                                                   
#   _____                           _           
#  / ____|                         (_)          
# | (___   ___ ___ _ __   __ _ _ __ _  ___  ___ 
#  \___ \ / __/ _ \ '_ \ / _` | '__| |/ _ \/ __|
#  ____) | (_|  __/ | | | (_| | |  | | (_) \__ \
# |_____/ \___\___|_| |_|\__,_|_|  |_|\___/|___/
#                                               

# Scenario 1 is LP selection:
Scenario1_LP_ELECTRONS    =  ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP
wwEleScenario1 = selectedRefPatElectrons.clone( cut = Scenario1_LP_ELECTRONS )

# Scenario 2 is LP selection + muon kink finder, no change in the electrons
Scenario2_LP_ELECTRONS    =  ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP
wwEleScenario2 = selectedRefPatElectrons.clone( cut = Scenario2_LP_ELECTRONS )

# Scenario 3 is LP selection + muon kink finder switching to new LH ele ID
# the loose ID is applied before the LH ID 
Scenario3_LH_ELECTRONS    =  ELE_BASE + " && " + ELE_ID_LOOSE_SMURFS + " && " + ELE_JAMLH_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP
wwEleScenario3 = selectedRefPatElectrons.clone( cut = Scenario3_LH_ELECTRONS )

# Scenario 4 is LP selection + muon kink finder switching to new BDT ele ID 
Scenario4_BDT_ELECTRONS   =  ELE_BASE + " && " + ELE_ID_LOOSE_SMURFS + " && " + ELE_BDT_ID_SMURF + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP
wwEleScenario4 = selectedRefPatElectrons.clone( cut = Scenario4_BDT_ELECTRONS)

# Scenario 5 is a loose loose tree to do the Wjets studies
Scenario5_LOOSE_ELECTRONS =  ELE_BASE + " && " + ELE_ID_LOOSE_SMURFS + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP
wwEleScenario5 = selectedRefPatElectrons.clone( cut = Scenario5_LOOSE_ELECTRONS )

# 
#  ______           _ 
# |  ____|         | |
# | |__   _ __   __| |
# |  __| | '_ \ / _` |
# | |____| | | | (_| |
# |______|_| |_|\__,_|
#                     

wwElectronSequence = cms.Sequence(  
    eleBDTSelection *
    wwEleMatch *
#     wwEleIDLHL *
#     wwEleISOLHL *
#     wwEleCONVLHL *
#     wwEleIPLHL *
#    wwEleIDLHT *
#    wwEleISOLHT *
#    wwEleCONVLHT *
#    wwEleIPLHT *
    wwEleIDMerge *
    wwEleISOMerge *
    wwEleCONVMerge *
    wwEleIPMerge *
    wwEleIDMerge2 *
    wwEleISOMerge2 *
    wwEleCONVMerge2 *
    wwEleIPMerge2 
#     wwEleISOPFLHT *
#     wwEleCONVPFLHT *
#     wwEleIPPFLHT 
#     wwEleIDCBL *
#     wwEleISOCBL *
#     wwEleCONVCBL *
#     wwEleIPCBL *
#     wwEleIDCBT *
#     wwEleISOCBT *
#     wwEleCONVCBT *
#     wwEleIPCBT 
)


