import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_10_1_7uf.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_11_1_kXN.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_12_1_Q0l.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_13_1_dGR.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_14_1_T6p.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_15_1_jbB.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_16_1_wHW.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_17_1_ZUK.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_18_1_nts.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_19_1_tnl.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_1_1_e8m.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_20_1_5xz.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_21_1_RHh.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_22_1_rhs.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_23_1_wTX.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_24_1_ols.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_25_1_h9F.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_26_1_s5S.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_27_1_fqx.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_28_1_jGC.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_29_1_pZB.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_2_1_j90.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_30_1_Hrb.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_31_1_ent.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_32_1_2Ze.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_33_1_YZK.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_34_1_HC2.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_35_1_zNc.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_36_1_Uiz.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_37_1_PMV.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_38_1_eNO.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_39_1_Hqa.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_3_1_tFR.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_40_1_Znt.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_41_1_KPK.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_42_1_tKW.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_4_1_kl5.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_5_1_wxb.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_6_1_pk6.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_7_1_8EI.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_8_1_8OV.root',
    'file:/nfs/bluearc/group/skims/hww/Wenu/hwwSkim_9_1_bth.root',
] );


secFiles.extend( [
               ] )

