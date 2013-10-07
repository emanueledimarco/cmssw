import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_10_1_eD5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_11_1_Fto.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_1_1_rF5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_2_1_4iq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_3_1_d7r.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_4_1_LSp.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_5_1_N1q.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_6_1_r03.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_7_1_LFd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_8_1_yMA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_t/hwwSkim_9_1_b6T.root",
] );


secFiles.extend( [
               ] )

