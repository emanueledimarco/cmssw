import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/ZZ/hwwSkim_1_1_S5w.root',
    'file:/nfs/bluearc/group/skims/hww/ZZ/hwwSkim_2_1_CIC.root',
    'file:/nfs/bluearc/group/skims/hww/ZZ/hwwSkim_3_1_wbq.root',
] );


secFiles.extend( [
               ] )

