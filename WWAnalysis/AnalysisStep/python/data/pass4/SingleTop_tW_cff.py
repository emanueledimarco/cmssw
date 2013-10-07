import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_10_1_Zj8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_11_1_zue.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_12_1_Ss6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_13_1_myu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_14_1_WeV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_15_1_Wpo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_16_1_xaI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_17_1_7Eh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_18_1_L2I.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_19_1_a9s.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_1_1_YjI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_20_1_AA4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_21_1_wmh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_22_1_Ww8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_23_1_wEx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_24_1_lho.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_25_1_6HJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_26_1_wVX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_27_1_QAz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_28_1_XbU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_29_1_faz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_2_1_mkj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_30_1_hML.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_31_1_D9o.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_32_1_xjK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_33_1_uMb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_34_1_4OW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_35_1_rI5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_36_1_K1d.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_37_1_Ume.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_38_1_Rqj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_39_1_gSV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_3_1_sBn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_40_1_8WY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_41_1_yiV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_4_1_jz4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_5_1_uPb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_6_1_iQr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_7_1_8B7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_8_1_zc6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/SingleTop_tW/hwwSkim_9_1_gcC.root",
] );


secFiles.extend( [
               ] )

