import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_10_1_XKJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_11_2_tjb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_1_1_gte.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_12_1_2Oj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_13_1_uUT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_14_1_hUC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_15_1_XnL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_16_1_MFk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_17_1_Syr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_18_1_rVI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_19_1_xtM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_20_1_evp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_21_1_gsf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_2_1_61W.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_3_1_2c0.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_4_1_l7q.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_5_1_XZZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_6_1_TLQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_7_1_wnQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wgamma/hwwSkim_8_1_PqZ.root",
] );


secFiles.extend( [
               ] )

