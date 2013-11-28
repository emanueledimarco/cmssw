import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_10_1_D1A.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_11_1_xsC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_1_1_vel.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_12_1_iru.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_13_1_gDJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_14_1_XlL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_15_1_oGo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_16_1_FRH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_17_1_jLR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_18_1_jig.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_19_1_q7H.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_20_1_bSm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_21_1_Nxp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_2_1_zoo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_3_1_jLE.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_4_1_FJa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_5_1_XKZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_6_1_2Qn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_7_1_alO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_8_1_HYm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wtaunu/hwwSkim_9_1_Ht2.root",
] );


secFiles.extend( [
               ] )

