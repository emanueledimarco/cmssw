import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_10_1_ZLy.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_11_1_dvb.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_1_1_up6.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_2_1_txj.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_3_1_T7q.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_4_1_8ej.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_5_1_GX0.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_6_1_5Oo.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_7_1_CHN.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_8_1_ChG.root',
    'file:/nfs/bluearc/group/skims/hww/SingleTop_t/hwwSkim_9_1_Lnj.root',
] );


secFiles.extend( [
               ] )

