import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_10_1_VJ3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_1_1_TrB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_2_1_5NW.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_3_1_fQl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_4_1_imc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_5_1_8Fc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_6_1_vU5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_7_1_JuC.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_8_1_W0H.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//WZ/hwwSkim_9_1_kKQ.root",
] );


secFiles.extend( [
               ] )

