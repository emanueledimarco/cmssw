import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_10_1_6jt.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_11_1_yXu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_1_1_5bW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_2_1_Icw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_3_1_gnv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_4_1_AwT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_5_1_COL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_6_1_jfk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_7_1_1WA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_8_1_wuP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//ZZ/hwwSkim_9_1_Qd4.root",
] );


secFiles.extend( [
               ] )

