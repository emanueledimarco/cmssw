import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_10_2_e2D.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_12_3_fQx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_13_3_0vU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_14_2_VRg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_15_2_TGs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_16_2_cgZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_17_2_Ark.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_18_2_YjW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_1_3_U0p.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_20_2_AUW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_21_3_ItT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_22_2_ijC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_28_2_CDH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_30_3_VcT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_33_2_G1U.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_34_3_47q.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_35_2_9rz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_36_2_oL9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_37_2_Y2b.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_38_2_Xyv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_39_2_17p.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_3_3_pDw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_40_2_MgV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_41_4_TDZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_42_3_tFL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_43_3_F9h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_45_2_FcW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_46_2_FPX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_47_2_Fnk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_48_2_zlc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_49_2_rcI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_4_3_f2g.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_50_2_0Wb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_51_3_Gih.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_52_2_Bd3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_53_3_9IK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_54_2_duX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_55_2_3ny.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_56_2_SKM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_57_2_KsP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_58_2_haL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_5_2_21I.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_60_3_2hP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_61_3_DQj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_62_3_JDR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_63_3_jt0.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_6_2_tdJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_7_3_FRL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_8_3_8BO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma_2l/hwwSkim_9_2_FV1.root",
] );


secFiles.extend( [
               ] )

