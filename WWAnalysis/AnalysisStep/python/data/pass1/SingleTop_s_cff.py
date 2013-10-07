import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_1_1_wEE.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_2_1_QNL.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_3_1_1uT.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_4_1_Z8D.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_5_1_GBJ.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_6_1_JaD.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_7_1_0Ac.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_8_1_1mF.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_s/hwwSkim_9_1_b6Q.root',
] );


secFiles.extend( [
               ] )

