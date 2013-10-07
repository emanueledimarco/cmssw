import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_10_2_dXn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_11_1_SJH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_1_1_on3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_13_1_O2I.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_14_2_ZTn.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_15_2_HwG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_16_2_b6R.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_17_1_LW7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_18_1_p7N.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_19_2_pyb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_20_2_4cb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_2_1_OLr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_3_1_S9H.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_4_1_Jes.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_5_3_rIe.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_6_2_VbR.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_7_1_2uA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_8_3_A7h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Wmunu/hwwSkim_9_1_4AP.root",
] );


secFiles.extend( [
               ] )

