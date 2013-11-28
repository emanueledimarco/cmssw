import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_10_1_4XI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_11_1_JqQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_12_1_pyL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_13_1_G1T.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_14_1_Hnp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_15_1_fWQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_16_1_RBG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_17_1_FyW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_18_1_Cfs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_19_1_0zb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_1_1_z6K.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_20_1_aZp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_21_1_Wqf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_22_1_X23.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_23_1_UMP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_24_1_YVf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_25_1_sPR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_26_1_Wrm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_27_1_ldR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_28_1_hnQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_29_1_uMG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_2_1_yJZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_30_1_GjA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_31_1_M2H.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_32_1_98h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_33_1_0xY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_34_1_6cB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_35_1_9qG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_36_1_Bwj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_37_1_Uq9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_38_1_uhV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_39_1_n4X.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_3_1_Tql.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_40_1_vSD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_41_1_1o9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_4_1_Y6g.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_5_1_2Sr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_6_1_fMb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_7_1_0CR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_8_1_uq5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_2l/hwwSkim_9_1_1fg.root",
] );


secFiles.extend( [
               ] )

