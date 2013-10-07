import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/Zgamma/hwwSkim_1_1_vgr.root',
    'file:/nfs/bluearc/group/skims/hww/Zgamma/hwwSkim_2_1_BxR.root',
    'file:/nfs/bluearc/group/skims/hww/Zgamma/hwwSkim_3_1_h88.root',
] );


secFiles.extend( [
               ] )

