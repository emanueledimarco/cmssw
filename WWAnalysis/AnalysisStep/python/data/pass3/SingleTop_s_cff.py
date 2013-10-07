import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_10_1_4Ty.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_11_1_mY6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_1_1_SWn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_12_1_mD0.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_13_1_TF8.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_14_1_DuJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_15_1_Yoj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_16_1_Rhj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_17_1_suL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_18_1_bYx.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_19_1_ffN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_20_1_Nhb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_21_1_BdH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_2_1_XmY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_3_1_f6F.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_4_1_Xf1.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_5_1_dnr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_6_1_DKT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_7_1_Wdt.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_8_1_QAr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_s/hwwSkim_9_1_fLx.root",
] );


secFiles.extend( [
               ] )

