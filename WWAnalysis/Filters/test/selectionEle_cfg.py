import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",
                     fileNames = readFiles, 
                     secondaryFileNames = secFiles,
#                     skipEvents = cms.untracked.uint32(324),
#                     eventsToProcess = cms.untracked.VEventRange('1:12356'),
#                     lumisToProcess = cms.untracked.VLuminosityBlockRange('1:23')

)
readFiles.extend( [  
#'file:/home/mangano/WW/TTbarJets_-START3X_V26_S09.root',

"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_10_1_u74.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_11_1_lXF.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_12_1_ozc.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_13_1_Pmh.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_1_1_u0x.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_2_1_KQI.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_3_1_fTy.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_4_1_XjO.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_5_1_j8A.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_6_1_xKg.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_7_1_IKi.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_8_1_33L.root",
"file:/data/gpetrucc/7TeV/ww/store/user/gpetrucc/WW/GenLevelSkim_Spring10-START3X_V26_S09-v1/2d32bcc4f128a3b5966c4cd81a8a2e76/test_9_1_Rxw.root"

#"rfio:/castor/cern.ch/user/m/mangano/newIP/test_1_1_iVV.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/test_2_1_YKG.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/test_3_1_pIN.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/test_4_1_epg.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/test_5_1_RT7.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/test_6_1_JG1.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/test_7_1_KZc.root",

#'file:/data/mangano/MC/WW/WW_START3X_V26_S09.root',
]);



secFiles.extend( [ ]);


process = cms.Process("WWSel2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START36_V9::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = source


#======================================
# HLT
#======================================
process.selectionHLT = cms.EDFilter("HLTHighLevel",
                        TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI"),
                        HLTPaths = cms.vstring('HLT_Mu9','HLT_Ele15_LW_L1R'), 
#                        HLTPaths = cms.vstring('HLT_Ele15_LW_L1R'), 
                        eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                        andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                        throw = cms.bool(True)    # throw exception on unknown path names
                               )


#process.TFileService = cms.Service("TFileService", 
#      fileName = cms.string("SET_OUTPUT"),
#      closeFileFast = cms.untracked.bool(True)
#)


process.load('WWAnalysis.Filters.electronSelections_cff')
process.load('WWAnalysis.Filters.muonSelections_cff')

process.load('WWAnalysis.Filters.daugtherListCleaner_cfi')
process.daugtherListCleaner.src = cms.InputTag("diEleSel2")

process.load('WWAnalysis.Filters.metFilter_cfi')
process.load('WWAnalysis.Filters.jetVetoFilter_cfi')
process.load('WWAnalysis.Filters.softMuonVeto_cfi')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

#process.Tracer = cms.Service("Tracer")
#process.contentAna = cms.EDAnalyzer("EventContentAnalyzer")

process.p2 = cms.Path(process.selectionHLT*
                      process.elePlusSelection1*process.eleMinusSelection1*     #sel1  
                      process.elePlusIPSelection*process.eleMinusIPSelection*   #sel2
                      process.elePlusIsoSelection*process.eleMinusIsoSelection* #sel3
                      process.elePlusID*process.eleMinusID*                     #sel4
                      process.convRejectionPlus*process.convRejectionMinus*     #sel5 (dummy, at the moment)
                      process.metSel1*                                          #sel6
                      process.diEleSel1*process.diEleFilter1*                   #sel7
                      process.diEleSel2*                                        #sel8
                      process.daugtherListCleaner*                             
                      process.metFilter*                                        #sel9
                      process.eleForVetoSequence*
                      process.muForVetoSequence*
                      process.jetVetoFilter*                                    #sel10
                      process.eleSoftMuonVeto*                                  #sel11
                      process.extraLeptonVetoForEle                             #sel12
                      )

# ---- endPath ----
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    outputCommands = cms.untracked.vstring('drop *',
                                           'keep *_*_*_WWSel'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p2"))
)
#process.end = cms.EndPath(process.out)


