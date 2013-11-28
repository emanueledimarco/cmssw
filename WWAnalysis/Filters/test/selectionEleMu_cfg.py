import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",
                     fileNames = readFiles, 
                     secondaryFileNames = secFiles,
#                     skipEvents = cms.untracked.uint32(324),
)
readFiles.extend( [  
#"rfio:/castor/cern.ch/user/m/mangano/newIP/testEleMu2_1_1_gP8.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/testEleMu2_2_1_Z3J.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/testEleMu2_3_1_A0q.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/testEleMu2_4_1_CIj.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/testEleMu2_5_1_RjY.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/testEleMu2_6_1_NJp.root",
#"rfio:/castor/cern.ch/user/m/mangano/newIP/testEleMu2_7_1_eHf.root",

"file:/data/mangano/MC/WW/testEleMu2_1_1_gP8.root",
"file:/data/mangano/MC/WW/testEleMu2_2_1_Z3J.root",
"file:/data/mangano/MC/WW/testEleMu2_3_1_A0q.root",
"file:/data/mangano/MC/WW/testEleMu2_4_1_CIj.root",
"file:/data/mangano/MC/WW/testEleMu2_5_1_RjY.root",
"file:/data/mangano/MC/WW/testEleMu2_6_1_NJp.root",
"file:/data/mangano/MC/WW/testEleMu2_7_1_eHf.root",

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
process.load('WWAnalysis.Filters.eleMuSelections_cff')

process.load('WWAnalysis.Filters.daugtherListCleaner_cfi')
process.daugtherListCleaner.src = cms.InputTag("diLeptSel2")

process.load('WWAnalysis.Filters.metFilter_cfi')
process.metFilter.skipMEtCut = cms.bool(True)

process.load('WWAnalysis.Filters.jetVetoFilter_cfi')
process.load('WWAnalysis.Filters.softMuonVeto_cfi')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

#process.Tracer = cms.Service("Tracer")
#process.contentAna = cms.EDAnalyzer("EventContentAnalyzer")

process.p2 = cms.Path(process.selectionHLT*
                      process.eleSelection1*process.muSelection1*
                      process.selection1Pairs*process.selection1PairsCounter* # cut1
                      process.eleIPSelection*process.muIPSelection*
                      process.ipSelPairs*process.ipSelPairsCounter*           # cut2
                      process.eleIsoSelection*process.muIsoSelection*
                      process.isoSelPairs*process.isoSelPairsCounter*         # cut3
                      process.eleID*process.muID*
                      process.idSelPairs*process.idSelPairsCounter*           # cut4
                      process.eleConvRejection*
                      process.convSelPairs*process.convSelPairsCounter*       # cut5
                      process.metSelLept*                                     # cut6
                      process.diLeptSel1*process.diLeptFilter1*               # cut7
                      process.diLeptSel2*                                     # cut8
                      process.daugtherListCleaner*                            
                      process.metFilter*                                      # cut9
                      process.eleForVetoSequence*
                      process.muForVetoSequence*
                      process.jetVetoFilter*                                  # cut10
                      process.leptSoftMuonVeto*                               # cut11
                      process.extraLeptonVetoForLept                          # cut12
                      )

# ---- endPath ----
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    outputCommands = cms.untracked.vstring('drop *',
                                           'keep *_diLeptSel1_*_WWSel2',
                                           'keep *_diLeptSel2_*_WWSel2'
                                           ),
#    outputCommands = process.AODEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p2"))
)
#process.end = cms.EndPath(process.out)


