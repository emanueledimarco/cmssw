import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_10_1_kwU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_1_1_1J2.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_11_1_L2s.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_2_1_dzy.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_3_1_pZm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_4_1_EZj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_5_1_wB6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_6_1_uQC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_7_1_hFU.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_8_1_CO5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//TTbar_2l/hwwSkim_9_1_V8c.root",
] );


secFiles.extend( [
               ] )

