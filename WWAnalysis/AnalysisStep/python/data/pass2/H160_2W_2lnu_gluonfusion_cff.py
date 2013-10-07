import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_10_1_9NM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_1_1_nVh.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_2_1_UQX.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_3_1_fJj.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_4_1_8sB.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_5_1_zFb.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_6_1_McG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_7_1_hKo.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_8_1_u7D.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass002//H160_2W_2lnu_gluonfusion/hwwSkim_9_1_lDo.root",
] );


secFiles.extend( [
               ] )

