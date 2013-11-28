import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_10_1_wI6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_1_1_W4A.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_2_1_gTu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_3_1_1KH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_4_1_yO2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_5_1_nDG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_6_1_JQW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_7_1_gTQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_8_1_khV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WW_signal/hwwSkim_9_1_qrA.root",
] );


secFiles.extend( [
               ] )

