import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_10_1_lsT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_11_1_s5z.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_12_1_sSZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_13_1_YWH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_14_1_YVP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_15_1_6q8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_16_1_Lwz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_17_1_8eU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_18_1_Spl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_19_1_4gM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_1_1_qXM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_20_1_obs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_21_1_kPQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_22_1_KfJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_23_1_mol.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_24_1_tDc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_25_1_1HR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_26_1_NGL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_27_1_9up.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_28_1_gNf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_29_1_jnr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_2_1_brp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_30_1_CLY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_31_1_1kf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_32_1_y2P.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_33_1_CtX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_34_1_xNB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_35_1_oN9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_36_1_bJ7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_37_1_QVi.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_38_1_BvU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_39_1_beQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_3_1_Fex.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_40_1_HQY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_41_1_C9W.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_4_1_CXa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_5_1_qyC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_6_1_Cnw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_7_1_qW5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_8_1_Lgb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_t/hwwSkim_9_1_0Vv.root",
] );


secFiles.extend( [
               ] )

