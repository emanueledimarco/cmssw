import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_10_2_xXZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_11_2_AnE.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_1_1_glv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_12_2_cZe.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_13_1_2ks.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_14_1_MS8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_15_2_0lZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_16_1_P0x.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_17_5_P9s.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_18_5_L2u.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_19_2_DGK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_20_2_mwc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_21_2_Qob.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_2_2_kCd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_3_2_z2q.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_4_2_s7c.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_5_1_fDW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_6_2_O8A.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_7_1_p0A.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_8_1_GX7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//ZZ/hwwSkim_9_1_dZc.root",
] );


secFiles.extend( [
               ] )

