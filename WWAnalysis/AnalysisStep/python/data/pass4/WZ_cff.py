import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_10_1_jUO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_11_1_I71.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_12_1_Igf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_13_1_8qo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_14_1_nla.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_15_1_Pms.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_16_1_WUu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_17_1_iCt.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_18_1_mi8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_19_1_2rE.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_1_1_oXL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_20_1_1Cz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_21_1_UkY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_22_1_qmN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_23_1_P6R.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_24_1_lkM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_25_1_4rh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_26_1_b15.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_27_1_G2b.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_28_1_YrY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_29_1_mSt.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_2_1_4Sx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_30_1_ZOC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_31_1_miN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_32_1_wUw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_33_1_GSN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_34_1_O7s.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_35_1_5JB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_36_1_Lz6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_37_1_A7h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_38_1_2T7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_39_1_iZZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_3_1_L2v.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_40_1_vZB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_4_1_qvs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_5_1_uoZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_6_1_Z71.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_7_1_8eR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_8_1_EPj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WZ/hwwSkim_9_1_ScH.root",
] );


secFiles.extend( [
               ] )

