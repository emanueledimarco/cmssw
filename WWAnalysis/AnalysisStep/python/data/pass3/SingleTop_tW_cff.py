import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_10_1_7j9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_11_1_IXA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_1_1_usZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_12_2_0rR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_13_1_7jK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_14_1_9iO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_15_1_5Bo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_16_1_U0s.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_17_1_WVs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_18_1_Vnb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_19_2_vcs.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_20_1_sEf.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_21_1_OoB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_2_1_UtC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_3_1_FbL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_4_1_WE2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_5_2_2oI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_6_1_AQP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_7_1_Ckn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_8_1_WB3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_tW/hwwSkim_9_1_qEw.root",
] );


secFiles.extend( [
               ] )

