import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_10_1_x29.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_11_1_tO7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_1_1_bLa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_2_1_3JQ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_3_1_BDz.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_4_1_Eij.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_5_1_i2F.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_6_1_SQ4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_7_1_NRG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_8_1_4ka.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_s/hwwSkim_9_1_3w7.root",
] );


secFiles.extend( [
               ] )

