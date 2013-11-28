import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [  
#'file:/home/mangano/WW/TTbarJets_-START3X_V26_S09.root',
#'file:/data/mangano/MC/WW/WW_START3X_V26_S09.root',
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
]);



secFiles.extend( [ ]);


process = cms.Process("WWSel")

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


process.load('WWAnalysis.Filters.muonSelections_cff')
process.load('WWAnalysis.Filters.electronSelections_cff')
process.load('WWAnalysis.Filters.daugtherListCleaner_cfi')
process.load('WWAnalysis.Filters.metFilter_cfi')
process.load('WWAnalysis.Filters.jetVetoFilter_cfi')
process.load('WWAnalysis.Filters.softMuonVeto_cfi')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

process.p2 = cms.Path(process.selectionHLT*
                      process.muPlusSelection1*process.muMinusSelection1*      #sel1  
                      process.muPlusIPSelection*process.muMinusIPSelection*    #sel2
                      process.muPlusIsoSelection*process.muMinusIsoSelection*  #sel3
                      process.muPlusID*process.muMinusID*                      #sel4
#                      process.conversionRejection*                             #sel5
                      process.metSel1*                                         #sel6
                      process.diMuonSel1*process.diMuonFilter1*                #sel7
                      process.diMuonSel2*                                      #sel8
                      process.daugtherListCleaner*                             
                      process.metFilter*                                       #sel9
                      process.eleForVetoSequence*
                      process.muForVetoSequence*
                      process.jetVetoFilter*                                   #sel10
                      process.softMuonVeto*                                     #sel11
                      process.extraLeptonVetoForMu                             #sel12
                      )

# ---- endPath ----
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    outputCommands = cms.untracked.vstring('drop *',
                                           'keep *_*_*_WWSel'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p2"))
)
process.end = cms.EndPath(process.out)


