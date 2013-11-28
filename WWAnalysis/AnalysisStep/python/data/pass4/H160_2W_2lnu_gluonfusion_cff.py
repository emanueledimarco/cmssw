import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_10_1_tSa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_11_1_npF.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_12_1_Qa3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_13_1_NcJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_14_1_cxu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_15_1_Gg5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_16_1_eXa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_17_1_Qnu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_18_1_n5b.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_19_1_8Ie.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_1_1_kpS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_20_1_E8y.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_21_1_lCv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_22_1_IMM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_23_1_YH2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_24_1_vwp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_25_1_uLi.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_26_1_CGB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_27_1_ZxK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_28_1_tZS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_29_1_sQq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_2_1_9bg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_30_1_eFz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_31_1_xbz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_32_1_zU7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_33_1_F9V.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_34_1_e5Q.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_35_1_Iom.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_36_1_3Hz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_37_1_puV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_38_1_6YM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_39_1_4xG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_3_1_A7j.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_40_1_zft.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_41_1_UjY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_4_1_4an.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_5_1_37E.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_6_1_ZaL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_7_1_FFh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_8_1_lWk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/H160_2W_2lnu_gluonfusion/hwwSkim_9_1_ovN.root",
] );


secFiles.extend( [
               ] )

