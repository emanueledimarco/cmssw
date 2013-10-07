import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/gg_WW_2l2nu/hwwSkim_1_1_dLO.root',
    'file:/nfs/bluearc/group/skims/hww/gg_WW_2l2nu/hwwSkim_2_1_NZG.root',
    'file:/nfs/bluearc/group/skims/hww/gg_WW_2l2nu/hwwSkim_3_1_qBY.root',
] );


secFiles.extend( [
               ] )

