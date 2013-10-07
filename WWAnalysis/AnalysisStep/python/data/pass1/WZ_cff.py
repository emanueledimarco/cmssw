import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/WZ/hwwSkim_1_1_1gB.root',
    'file:/nfs/bluearc/group/skims/hww/WZ/hwwSkim_2_1_Eq4.root',
    'file:/nfs/bluearc/group/skims/hww/WZ/hwwSkim_3_1_bzW.root',
] );


secFiles.extend( [
               ] )

