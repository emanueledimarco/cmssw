import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_10_1_xpM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_11_1_0Xz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_12_1_CyG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_13_1_WhS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_14_1_yUp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_15_1_6hN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_16_1_DHC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_17_1_4gN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_18_1_DpZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_19_1_h23.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_1_1_zx3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_20_1_lTj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_21_1_AMV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_22_1_KaY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_23_1_37g.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_24_1_hSH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_25_1_lQj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_26_1_f2T.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_27_1_aSv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_28_1_Nh5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_29_1_wiS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_2_1_sVx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_30_1_o1a.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_31_1_Z47.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_32_1_207.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_33_1_YnZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_34_1_3AV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_35_1_lGN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_36_1_XD8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_37_1_o6y.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_38_1_bGQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_39_1_Ws4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_3_1_jmb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_40_1_n5A.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_41_1_HpM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_4_1_tGi.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_5_1_D9w.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_6_1_HRp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_7_1_A2M.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_8_1_YXK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Zgamma/hwwSkim_9_1_W3Q.root",
] );


secFiles.extend( [
               ] )

