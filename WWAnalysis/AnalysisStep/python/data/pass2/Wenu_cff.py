import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_10_1_vfj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_11_1_KQS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_1_1_yB8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_2_1_TJG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_3_1_H4a.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_4_1_iW7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_5_1_nqO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_6_1_4jp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_7_1_qqp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_8_1_WFN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wenu/hwwSkim_9_1_cYv.root",
] );


secFiles.extend( [
               ] )

