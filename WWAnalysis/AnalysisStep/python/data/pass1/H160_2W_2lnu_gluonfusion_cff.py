import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/H160_2W_2lnu_gluonfusion/hwwSkim_1_1_LWL.root',
    'file:/nfs/bluearc/group/skims/hww/H160_2W_2lnu_gluonfusion/hwwSkim_2_1_01g.root',
    'file:/nfs/bluearc/group/skims/hww/H160_2W_2lnu_gluonfusion/hwwSkim_3_1_sXd.root',
] );


secFiles.extend( [
               ] )

