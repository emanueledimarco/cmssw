import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_10_1_1FD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_11_1_5rM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_12_1_T77.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_1_2_maY.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_13_2_Snk.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_14_1_RV2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_15_2_Tc5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_16_1_VAB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_17_2_4It.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_18_1_Dpd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_19_2_uL2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_20_1_Row.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_21_2_8fA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_2_2_6kl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_3_2_RlI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_4_2_UhW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_5_1_ZRI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_6_1_HYc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_7_1_ofg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//Zgamma/hwwSkim_8_2_pwN.root",
] );


secFiles.extend( [
               ] )

