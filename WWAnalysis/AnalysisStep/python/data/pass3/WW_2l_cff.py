import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_10_2_eLv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_11_2_38X.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_1_2_1j7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_12_2_b9u.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_13_2_FZ3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_14_2_gYo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_15_2_Ear.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_16_2_Kp1.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_17_2_xAW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_18_2_RGA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_19_2_kkr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_20_2_p7T.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_21_2_ejT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_2_2_plP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_3_2_xwm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_4_2_xon.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_5_2_VNP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_6_2_inw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_7_2_IrT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_8_2_E82.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WW_2l/hwwSkim_9_2_72a.root",
] );


secFiles.extend( [
               ] )

