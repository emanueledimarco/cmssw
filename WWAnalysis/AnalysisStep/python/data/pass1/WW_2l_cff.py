import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/WW_2l/hwwSkim_1_1_rrw.root',
    'file:/nfs/bluearc/group/skims/hww/WW_2l/hwwSkim_2_1_Eju.root',
    'file:/nfs/bluearc/group/skims/hww/WW_2l/hwwSkim_3_1_EHI.root',
    'file:/nfs/bluearc/group/skims/hww/WW_2l/hwwSkim_4_1_B9X.root',
    'file:/nfs/bluearc/group/skims/hww/WW_2l/hwwSkim_5_1_Tyl.root',
    'file:/nfs/bluearc/group/skims/hww/WW_2l/hwwSkim_6_1_6lE.root',
] );


secFiles.extend( [
               ] )

