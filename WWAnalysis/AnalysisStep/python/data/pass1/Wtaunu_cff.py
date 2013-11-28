import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_10_2_mfD.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_11_1_u4b.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_12_1_mj4.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_13_1_1Ib.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_14_1_eJ9.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_15_1_0d9.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_16_1_nW8.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_17_1_kc7.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_18_1_GK8.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_19_1_jsH.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_1_1_dt8.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_20_1_rZN.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_21_1_WTA.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_22_1_MWj.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_23_1_PUl.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_24_1_vM9.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_25_1_4GG.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_26_1_N3i.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_27_1_N27.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_28_1_Cvq.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_29_1_qIp.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_2_1_caH.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_30_1_3Hv.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_31_1_eBv.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_32_1_9bn.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_33_1_IEc.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_34_1_jt3.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_35_1_Jue.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_36_1_MAV.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_37_1_Fgr.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_38_1_5mE.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_39_1_AET.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_3_1_CeZ.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_40_1_afL.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_41_1_vL5.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_42_1_jdF.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_43_1_v0d.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_44_1_Gyj.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_4_1_wiL.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_5_1_mfz.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_6_1_h61.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_7_1_1Ko.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_8_1_Y8h.root',
    'file:/nfs/bluearc/group/skims/hww/Wtaunu/hwwSkim_9_1_DiP.root',
] );


secFiles.extend( [
               ] )

