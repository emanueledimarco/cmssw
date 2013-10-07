import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_10_1_xpG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_1_1_8ty.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_2_1_ZM1.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_3_1_OLr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_4_1_wRF.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_5_1_8kO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_6_1_q4q.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_7_1_KqB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_8_1_p6k.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wgamma/hwwSkim_9_1_wju.root",
] );


secFiles.extend( [
               ] )

