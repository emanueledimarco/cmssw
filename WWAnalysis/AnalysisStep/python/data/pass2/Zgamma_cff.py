import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_10_1_Kq6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_1_1_vpC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_2_1_HjM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_3_2_f5O.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_4_1_SqZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_5_1_ziH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_6_1_S9C.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_7_1_0Xl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_8_1_3Sw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Zgamma/hwwSkim_9_1_7RZ.root",
] );


secFiles.extend( [
               ] )

