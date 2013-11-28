import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_10_1_YkW.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_11_1_1tY.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_12_1_FqK.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_13_1_NrN.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_14_1_m16.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_15_1_Ahw.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_16_1_LWF.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_17_1_cIa.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_18_1_U3u.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_19_1_VDh.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_1_1_7zl.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_20_1_13D.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_21_1_tie.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_22_1_vYv.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_23_1_H2X.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_24_1_tqo.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_25_1_MCB.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_26_1_9yz.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_27_1_pRs.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_28_1_UIv.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_29_1_Imi.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_2_1_v7t.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_30_1_hcT.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_31_1_8gu.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_32_1_or5.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_33_1_avA.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_34_1_RuD.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_35_1_liw.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_36_1_suV.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_37_1_fnU.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_38_1_0Fv.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_39_1_4jY.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_3_1_oFZ.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_40_1_V9f.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_41_1_a2t.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_42_1_j7I.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_43_1_yv2.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_4_1_TTK.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_5_1_kG2.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_6_1_AQ3.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_7_1_qhg.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_8_1_Sy6.root',
    'file:/nfs/bluearc/group/skims/hww/Wmunu/hwwSkim_9_1_vxD.root',
] );


secFiles.extend( [
               ] )

