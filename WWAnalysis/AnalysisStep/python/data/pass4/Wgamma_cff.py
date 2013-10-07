import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_10_1_B9R.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_11_1_lCX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_12_1_vsv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_13_1_aRs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_14_1_qJL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_15_1_9tC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_16_1_f0s.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_17_1_UXX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_18_1_NXp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_19_1_KS9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_1_2_8Uh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_20_1_LJb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_21_1_TVD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_22_1_1fb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_23_1_UHC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_24_1_KVN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_25_1_6Zl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_26_1_yXB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_27_1_l7m.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_28_1_zIj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_29_1_kf4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_2_1_p6r.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_30_1_2jL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_31_1_IFs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_32_1_1U8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_33_1_K6y.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_34_1_yFN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_35_1_fYk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_36_1_NTz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_37_1_IEE.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_38_1_TPf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_39_1_sSV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_3_1_yb5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_40_1_sso.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_41_1_puz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_4_1_MLp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_5_1_slg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_6_1_TqK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_7_1_WEj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_8_1_xdX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wgamma/hwwSkim_9_1_5In.root",
] );


secFiles.extend( [
               ] )

