import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_10_1_PJq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_11_1_FV7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_1_1_1jO.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_12_1_Fp5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_13_1_wmd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_14_1_Rt6.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_15_1_FgE.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_16_1_wMm.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_17_1_vKw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_18_1_Sge.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_19_1_WUr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_20_1_y6R.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_21_2_E46.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_2_1_vVJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_3_1_Cpu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_4_1_pBl.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_5_1_KDS.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_6_1_JcZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_7_1_I0f.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_8_1_wW7.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass003//H160_2W_2lnu_gluonfusion/hwwSkim_9_1_iSu.root",
] );


secFiles.extend( [
               ] )

