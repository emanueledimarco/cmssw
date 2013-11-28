import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_10_1_TPc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_11_1_wx7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_12_1_n3D.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_13_1_QY5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_14_1_B9C.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_15_1_Lg9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_16_1_Ajf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_17_1_aba.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_18_1_Bdc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_19_1_454.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_1_1_D31.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_20_1_llb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_21_1_15X.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_22_1_gZk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_23_1_E3W.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_24_1_PqO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_25_1_BOf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_26_1_gGJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_27_1_OsH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_28_1_9OJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_29_1_6Uy.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_2_1_rK7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_30_1_TSn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_31_1_Wz0.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_32_2_Nir.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_33_1_7LW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_34_1_xB2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_35_1_7XW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_36_1_O4p.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_37_1_YWU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_38_1_e1X.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_39_1_8XF.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_3_1_2uA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_40_1_sqz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_41_1_R4t.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_4_1_c9h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_5_1_lfd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_6_1_LQn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_7_1_dBL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_8_2_hSv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wtaunu/hwwSkim_9_1_05O.root",
] );


secFiles.extend( [
               ] )

