import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_10_1_NYE.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_11_1_Sf4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_12_1_kld.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_13_1_EQH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_14_1_WGS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_15_1_YIh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_16_1_GHw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_17_1_sHi.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_18_1_qe8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_19_1_d5e.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_1_1_cit.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_20_1_HPi.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_21_1_q97.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_22_1_vhH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_23_1_ByZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_24_1_z90.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_25_1_rPr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_26_1_BvK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_27_1_m9J.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_28_1_U9s.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_29_1_rAa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_2_1_t22.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_30_1_bbk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_31_1_5HR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_32_1_uoK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_33_1_7TL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_34_1_nO3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_35_1_FFF.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_36_1_Yf8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_37_1_qBN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_38_1_fBO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_39_1_Flq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_3_1_XPf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_40_1_ft1.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_41_1_Jp8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_4_1_evs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_5_1_Zo2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_6_1_Xpk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_7_1_F2N.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_8_1_7Qc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_s/hwwSkim_9_1_2s3.root",
] );


secFiles.extend( [
               ] )

