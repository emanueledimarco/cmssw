import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_10_1_I2R.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_11_1_WGq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_1_1_z5G.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_2_1_w4m.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_3_1_omi.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_4_1_4aG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_5_1_txt.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_6_1_yti.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_7_1_kyL.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_8_1_l9N.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//SingleTop_tW/hwwSkim_9_1_06S.root",
] );


secFiles.extend( [
               ] )

