import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_10_1_SCm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_11_1_Pis.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_1_1_4dV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_2_1_ZfR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_3_1_oiB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_4_1_ruj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_5_1_sWu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_6_1_HKr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_7_1_tTX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_8_1_DO0.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//Wmunu/hwwSkim_9_1_enT.root",
] );


secFiles.extend( [
               ] )

