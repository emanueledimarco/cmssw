import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_10_1_WRG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_11_1_kTQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_12_1_mLB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_13_1_UdM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_14_1_8ea.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_15_1_G6q.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_16_1_DEj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_17_1_hSL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_18_1_J22.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_19_1_itR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_1_1_1ZV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_20_1_h0l.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_21_1_cvB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_22_1_sBH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_23_1_EHf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_24_1_2Kp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_25_1_pMg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_26_1_Udu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_27_1_WeK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_28_1_Lum.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_29_1_ysz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_2_1_xYn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_30_1_0nm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_31_1_edG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_32_1_ahx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_33_1_6aQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_34_1_avi.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_35_1_Y8o.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_36_1_lbs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_37_1_3cs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_38_1_DWI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_39_1_kUw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_3_1_1zf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_40_1_gA3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_41_1_VDO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_4_1_rwp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_5_1_LqK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_6_1_HlV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_7_1_poz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_8_1_xxI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wenu/hwwSkim_9_1_j4F.root",
] );


secFiles.extend( [
               ] )

