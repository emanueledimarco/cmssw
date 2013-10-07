import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_10_1_3Wc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_11_2_7lb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_1_1_qek.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_2_1_F4i.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_3_2_C6i.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_4_1_wXT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_5_2_Qeu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_6_2_08i.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_7_2_GrO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_8_2_hjM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_2l/hwwSkim_9_2_LD1.root",
] );


secFiles.extend( [
               ] )

