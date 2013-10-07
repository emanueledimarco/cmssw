import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_10_1_503.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_1_1_usv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_2_1_cBj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_3_1_8Hg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_4_1_FsS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_5_1_zR6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_6_1_2nG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_7_2_HeQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_8_1_NeY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/WW_signal/hwwSkim_9_1_O5F.root",
] );


secFiles.extend( [
               ] )

