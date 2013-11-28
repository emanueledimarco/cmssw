import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_10_2_4ck.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_11_2_bQV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_12_2_UVD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_13_2_yU2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_14_2_t9Z.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_15_2_EAR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_16_2_49r.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_17_2_gSo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_18_2_fr6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_19_2_cQr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_1_2_qiM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_20_2_Uh7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_21_2_nvD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_22_2_vox.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_23_2_Shx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_24_2_2S1.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_25_2_jzl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_26_2_OLN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_27_2_qi2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_28_2_1ui.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_29_2_2Lg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_2_2_pNM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_30_2_SwP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_31_2_kiW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_32_2_6Sr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_33_2_4X2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_34_2_FSN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_35_2_HOg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_36_2_jkL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_37_2_bD0.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_38_2_x4h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_39_2_DmW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_3_2_pZ1.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_40_2_Vmx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_41_2_lPJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_4_2_ZK8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_5_2_MVz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_6_2_Aqf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_7_2_57F.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_8_2_nTV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/gg_WW_2l2nu/hwwSkim_9_2_0eq.root",
] );


secFiles.extend( [
               ] )

