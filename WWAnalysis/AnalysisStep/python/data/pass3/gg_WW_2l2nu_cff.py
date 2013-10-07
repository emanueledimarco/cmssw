import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_10_1_APf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_11_2_dbM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_12_2_x7c.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_1_2_xWb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_13_1_DnH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_14_2_jOK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_15_2_n79.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_16_1_wbb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_17_2_njP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_18_2_zJe.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_19_2_9cd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_20_2_05p.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_21_2_glI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_2_2_uWV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_3_2_K1K.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_4_2_2iR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_5_2_wmY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_6_2_RhQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_7_5_Gt6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_8_1_BCf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//gg_WW_2l2nu/hwwSkim_9_2_EuY.root",
] );


secFiles.extend( [
               ] )

