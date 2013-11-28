import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_10_1_skn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_11_1_nmW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_1_1_N1U.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_12_1_QYg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_13_1_peT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_14_1_YFf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_15_1_AoO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_16_1_BwL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_17_1_uSy.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_18_1_cRD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_19_1_8kz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_20_1_ZGJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_21_1_SNM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_2_1_cnh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_22_1_31f.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_23_1_QJ9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_24_1_OGO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_25_1_eXd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_26_1_wRg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_27_1_wQV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_28_1_9b8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_29_1_K9D.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_30_1_DMx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_31_1_zyr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_3_1_ea6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_32_1_Sh7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_33_1_I8u.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_34_1_Vsk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_4_1_Uu1.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_5_1_hiT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_6_1_mHN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_7_1_TIq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_8_1_Fx4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma_2l/hwwSkim_9_1_NoA.root",
] );


secFiles.extend( [
               ] )

