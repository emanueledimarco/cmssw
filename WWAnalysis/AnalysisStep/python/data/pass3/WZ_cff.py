import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_10_1_2WQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_1_1_1KD.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_11_1_YHU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_12_1_FKQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_13_1_eim.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_14_2_qGT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_15_2_dE6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_16_1_t7l.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_17_2_OHg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_18_1_TAu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_19_2_hoN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_20_1_7QF.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_2_1_5YS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_3_2_R92.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_4_1_GuG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_5_2_84h.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_6_2_dJp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_7_1_Uj9.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_8_2_mJq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//WZ/hwwSkim_9_1_2rp.root",
] );


secFiles.extend( [
               ] )

