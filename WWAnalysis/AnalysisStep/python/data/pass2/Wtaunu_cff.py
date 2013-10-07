import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_10_1_PAn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_11_1_21J.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_1_1_Pom.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_2_1_06t.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_3_1_8fX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_4_1_mcM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_5_1_8tt.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_6_1_bXO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_7_1_3Xq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_8_1_M6x.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wtaunu/hwwSkim_9_1_u7k.root",
] );


secFiles.extend( [
               ] )

