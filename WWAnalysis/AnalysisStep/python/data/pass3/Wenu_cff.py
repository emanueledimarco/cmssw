import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_10_1_UOX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_11_1_Z4h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_1_1_Bb7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_12_1_dqd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_13_1_fTB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_14_1_1gA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_15_1_4uq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_17_1_XIu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_18_1_Rab.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_19_1_caQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_20_1_M5c.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_21_1_BTK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_2_1_Mn0.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_3_1_J99.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_4_1_uqu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_5_1_jkb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_6_1_Wjl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_7_1_YVc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_8_1_4ve.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wenu/hwwSkim_9_1_Vs1.root",
] );


secFiles.extend( [
               ] )

