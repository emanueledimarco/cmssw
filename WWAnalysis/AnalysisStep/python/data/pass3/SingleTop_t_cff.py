import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_10_1_yMe.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_11_1_ZzK.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_1_1_X7b.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_12_1_0Me.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_13_1_7qa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_14_1_74i.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_15_1_A2r.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_16_1_nFy.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_17_1_FWV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_18_1_pjh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_19_1_7mP.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_20_1_qXI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_21_1_HUN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_2_1_sMV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_3_1_jgM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_4_1_hSY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_5_1_HKA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_6_1_jnd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_7_1_Fhz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_8_1_gWo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//SingleTop_t/hwwSkim_9_1_2mX.root",
] );


secFiles.extend( [
               ] )

