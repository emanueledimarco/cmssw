import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_10_1_Wwe.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_11_1_2e6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_12_1_EF6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_13_1_bdF.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_14_1_Wi7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_15_1_Ogz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_16_1_13i.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_17_1_pRZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_18_1_4uu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_19_1_VA9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_1_1_v3M.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_20_1_q90.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_21_1_cLC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_22_1_Z2t.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_23_1_qAR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_24_1_c5G.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_25_1_hiB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_26_1_1UL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_27_1_Zf3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_28_1_sA3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_29_1_t67.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_2_1_4Tu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_30_1_4jx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_31_1_fei.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_32_1_GoX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_33_1_Upb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_34_1_9ip.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_35_1_w0D.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_36_1_01N.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_37_1_4ol.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_38_1_TqN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_39_1_6LD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_3_1_d1O.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_40_1_zZc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_41_1_ysE.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_4_1_Y1o.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_5_1_wNr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_6_1_d3T.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_7_1_QgW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_8_1_Arb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/ZZ/hwwSkim_9_1_OPL.root",
] );


secFiles.extend( [
               ] )

