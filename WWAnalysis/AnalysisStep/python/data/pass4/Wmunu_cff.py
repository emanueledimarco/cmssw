import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_10_1_8yv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_11_1_yMT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_12_1_aWJ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_13_2_q0Z.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_14_1_o7d.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_15_1_KS4.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_16_1_i4J.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_17_1_sg5.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_18_1_JOM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_19_1_ozy.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_1_1_Evg.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_20_1_TZv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_21_1_sor.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_22_1_HXc.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_23_1_JCZ.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_24_1_HDr.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_25_1_a4a.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_26_1_3BA.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_27_1_F4y.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_28_1_Jyw.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_29_1_W4X.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_2_1_sWV.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_30_1_iLM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_31_1_4KH.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_32_1_euv.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_33_1_UdT.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_34_1_8ot.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_35_1_Tic.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_36_1_jmN.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_37_1_Qb3.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_38_1_Bln.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_39_1_0Ei.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_3_1_UZd.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_40_1_AMG.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_41_1_oyM.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_4_1_JHa.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_5_1_xUu.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_6_1_QIq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_7_2_SAq.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_8_1_qmI.root",
"file:/nfs/bluearc/group/skims/hww/SkimPass004/Wmunu/hwwSkim_9_1_OTi.root",
] );


secFiles.extend( [
               ] )

