import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_10_1_MOl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_1_1_1Dc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_11_1_dnf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_2_1_oub.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_3_1_Kkn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_4_1_6Wb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_5_1_GAY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_6_1_xXl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_7_1_b4H.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_8_1_xLO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//gg_WW_2l2nu/hwwSkim_9_1_u2Y.root",
] );


secFiles.extend( [
               ] )

