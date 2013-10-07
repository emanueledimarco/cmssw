import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_10_1_5tr.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_1_1_UTO.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_2_1_MTQ.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_3_1_bfB.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_4_1_8GT.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_5_1_wRN.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_6_1_unB.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_7_1_LON.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_8_1_K2y.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_tW/hwwSkim_9_1_DYW.root',
] );


secFiles.extend( [
               ] )

