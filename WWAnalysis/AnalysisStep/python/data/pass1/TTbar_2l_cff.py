import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_10_1_49c.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_11_1_4Xu.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_12_1_pgp.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_13_1_OUj.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_14_1_o85.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_15_1_iTn.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_16_1_OIU.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_17_1_YaX.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_18_1_Lnu.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_19_1_QQp.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_1_1_IIV.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_20_1_9Ou.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_21_1_KD8.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_22_1_u5E.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_2_1_1Ok.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_3_1_myI.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_4_1_AGK.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_5_1_yUO.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_6_1_LTI.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_7_1_ILM.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_8_1_pEk.root',
    'file:/nfs/bluearc/group/skims/hww/TTbar_2l/hwwSkim_9_1_bqq.root',
] );


secFiles.extend( [
               ] )

