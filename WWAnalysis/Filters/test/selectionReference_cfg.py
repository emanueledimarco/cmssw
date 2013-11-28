import FWCore.ParameterSet.Config as cms

process = cms.Process("WWSel2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# import of standard configurations
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START36_V9::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 


process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(
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
    ),
)


#======================================
# HLT
#======================================
process.selectionHLT = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI"),
        HLTPaths = cms.vstring('HLT_Mu9','HLT_Ele15_LW_L1R'), 
        eventSetupPathsKey = cms.string(''), 
        andOr = cms.bool(True),            
        throw = cms.bool(True),
)

process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 40 && position.Rho < 2 && abs(position.Z) < 24"),
)

process.load('WWAnalysis.Filters.electronSelections_cff')
process.load('WWAnalysis.Filters.muonSelections_cff')
process.load('WWAnalysis.Filters.eleMuSelections_cff')

## No Missing Hit requirement
process.elePlusMissingHits = cms.EDFilter("GsfElectronSelector",
    src = cms.InputTag("convRejectionPlus"),
    cut = cms.string('gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0'),
    filter = cms.bool(True),
)
process.eleMinusMissingHits = process.elePlusMissingHits.clone(src = 'convRejectionMinus')
process.eleMissingHits      = process.elePlusMissingHits.clone(src = 'eleConvRejection')
process.diEleSel1NoMH = process.diEleSel1.clone(decay = 'elePlusMissingHits eleMinusMissingHits')
process.diEleSel2NoMH = process.diEleSel2.clone(src = 'diEleSel1NoMH')
process.diLeptSel1NoMH = process.diLeptSel1.clone(decay = 'eleMissingHits@+ muID@-')
process.diLeptSel2NoMH = process.diLeptSel2.clone(src = 'diLeptSel1NoMH')

process.load('WWAnalysis.Filters.daugtherListCleaner_cfi')
process.daugtherListCleanerDiMuon = process.daugtherListCleaner.clone(src = cms.InputTag("diMuonSel2"))
process.daugtherListCleanerDiEle  = process.daugtherListCleaner.clone(src = cms.InputTag("diEleSel2"))
process.daugtherListCleanerEleMu  = process.daugtherListCleaner.clone(src = cms.InputTag("diLeptSel2"))
process.daugtherListCleanerDiEleNoMH  = process.daugtherListCleaner.clone(src = cms.InputTag("diEleSel2NoMH"))
process.daugtherListCleanerEleMuNoMH  = process.daugtherListCleaner.clone(src = cms.InputTag("diLeptSel2NoMH"))


process.load('WWAnalysis.Filters.metFilter_cfi')
process.metFilterDiMuon = process.metFilter.clone(srcLeptons = "daugtherListCleanerDiMuon")
process.metFilterDiEle  = process.metFilter.clone(srcLeptons = "daugtherListCleanerDiEle" )
process.metFilterEleMu  = process.metFilter.clone(srcLeptons = "daugtherListCleanerEleMu", skipMEtCut = True )
process.metFilterDiMuonV1 = process.metFilterDiMuon.clone(projectedMETCut = 35,  skipMEtCut = True )
process.metFilterDiEleV1  = process.metFilterDiMuonV1.clone(srcLeptons = "daugtherListCleanerDiEleNoMH")

process.load('WWAnalysis.Filters.jetVetoFilter_cfi')
process.eleForVetoNoMH = process.eleMissingHits.clone(src = 'eleForVeto', filter = False)
process.eleForVetoSequenceNoMH = cms.Sequence(process.eleForVetoSequence + process.eleForVetoNoMH)

process.load('WWAnalysis.Filters.softMuonVeto_cfi')

process.DiMuSequence = cms.Sequence(
    process.muPlusSelection1*process.muMinusSelection1*        
    process.muPlusIPSelection*process.muMinusIPSelection*    
    process.muPlusIsoSelection*process.muMinusIsoSelection*  
    process.muPlusID*process.muMinusID*                      
    process.diMuonSel1*                
    process.diMuonSel2*                                      
    process.daugtherListCleanerDiMuon
)
process.EleSequence = cms.Sequence(
    process.elePlusSelection1*process.eleMinusSelection1*       
    process.elePlusIPSelection*process.eleMinusIPSelection*   
    process.elePlusIsoSelection*process.eleMinusIsoSelection* 
    process.elePlusID*process.eleMinusID*                     
    process.convRejectionPlus*process.convRejectionMinus 
)
process.DiEleSequence = cms.Sequence(
    process.EleSequence*
    process.diEleSel1*                   
    process.diEleSel2*                                        
    process.daugtherListCleanerDiEle                             
)
process.DiEleSequenceNoMH = cms.Sequence(
    process.EleSequence*
    process.elePlusMissingHits*process.eleMinusMissingHits*
    process.diEleSel1NoMH*                   
    process.diEleSel2NoMH*                                        
    process.daugtherListCleanerDiEleNoMH                         
)
process.EleMuSequence = cms.Sequence(
    process.eleSelection1*process.muSelection1*
    process.selection1Pairs*process.selection1PairsCounter* 
    process.eleIPSelection*process.muIPSelection*
    process.ipSelPairs*process.ipSelPairsCounter*           
    process.eleIsoSelection*process.muIsoSelection*
    process.isoSelPairs*process.isoSelPairsCounter*         
    process.eleID*process.muID*
    process.idSelPairs*process.idSelPairsCounter*           
    process.eleConvRejection*
    process.convSelPairs*process.convSelPairsCounter
)
process.DiLeptSequence = cms.Sequence(
    process.EleMuSequence*
    process.diLeptSel1*               
    process.diLeptSel2*                                     
    process.daugtherListCleanerEleMu
)
process.DiLeptSequenceNoMH = cms.Sequence(
    process.EleMuSequence*
    process.eleMissingHits* 
    process.diLeptSel1NoMH*               
    process.diLeptSel2NoMH*                                     
    process.daugtherListCleanerEleMuNoMH
)

process.jetVetoSequence = cms.Sequence(
    process.ak5PFJets*
    process.eleForVetoSequence*
    process.muForVetoSequence*
    process.jetVetoFilter                                   
)
process.jetVetoFilterNew = process.jetVetoFilter.clone(
    srcElectrons = 'eleForVetoNoMH',
    etMax = 25., etaMax = 5.
)
process.jetVetoSequenceNew = cms.Sequence(
    process.ak5PFJets*
    process.eleForVetoSequenceNoMH*
    process.muForVetoSequence*
    process.jetVetoFilterNew    
)

process.DiMuon_V0 = cms.Path(
    process.selectionHLT*
    process.DiMuSequence*
    process.metSel1*                                          
    process.metFilterDiMuon*                                        
    process.jetVetoSequence*
    process.softMuonVeto*                                     
    process.extraLeptonVetoForMu                             
)
process.DiEle_V0 = cms.Path(
    process.selectionHLT*
    process.DiEleSequence*
    process.metSel1*                                         
    process.metFilterDiEle* 
    process.jetVetoSequence*
    process.eleSoftMuonVeto*                                  
    process.extraLeptonVetoForEle                             
)

process.EleMu_V0 = cms.Path(
    process.selectionHLT*
    process.DiLeptSequence*
    process.metFilterEleMu*                                      
    process.jetVetoSequence*
    process.leptSoftMuonVeto*                               
    process.extraLeptonVetoForLept                          
)

process.DiMuon_V1a = cms.Path(
    process.selectionHLT*
    process.DiMuSequence*
    process.metSel1*                                          
    process.metFilterDiMuonV1*                                        
    process.jetVetoSequence*
    process.softMuonVeto*                                     
    process.extraLeptonVetoForMu                             
)
process.DiEle_V1a = cms.Path(
    process.selectionHLT*
    process.DiEleSequenceNoMH*
    process.metSel1*                                         
    process.metFilterDiEleV1* 
    process.jetVetoSequence*
    process.eleSoftMuonVeto*                                  
    process.extraLeptonVetoForEle                             
)
process.EleMu_V1a = cms.Path(
    process.selectionHLT*
    process.DiLeptSequenceNoMH*
    process.metFilterEleMu*                                      
    process.jetVetoSequence*
    process.leptSoftMuonVeto*                               
    process.extraLeptonVetoForLept                          
)

process.DiMuon_V1b = cms.Path(
    process.selectionHLT*
    process.DiMuSequence*
    process.metSel1*                                          
    process.metFilterDiMuonV1*                                        
    process.jetVetoSequenceNew*
    process.softMuonVeto*                                     
    process.extraLeptonVetoForMu                             
)
process.DiEle_V1b = cms.Path(
    process.selectionHLT*
    process.DiEleSequenceNoMH*
    process.metSel1*                                         
    process.metFilterDiEleV1* 
    process.jetVetoSequenceNew*
    process.eleSoftMuonVeto*                                  
    process.extraLeptonVetoForEle                             
)
process.EleMu_V1b = cms.Path(
    process.selectionHLT*
    process.DiLeptSequenceNoMH*
    process.metFilterEleMu*                                      
    process.jetVetoSequenceNew*
    process.leptSoftMuonVeto*                               
    process.extraLeptonVetoForLept                          
)



# ---- endPath ----
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    outputCommands = cms.untracked.vstring('drop *',
        'keep *_diLeptSel1_*_WWSel2',
        'keep *_diLeptSel2_*_WWSel2'
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring(
        "DiMuon_V0", "DiEle_V0", "DiEleMu_V0",
        "DiMuon_V1a","DiEle_V1a","DiEleMu_V1a",
        "DiMuon_V1b","DiEle_V1b","DiEleMu_V1b",
    ))
)
#process.end = cms.EndPath(process.out)
process.outDebug = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("selectedReference.root"),
    outputCommands = cms.untracked.vstring(
        'keep *', 
        'drop *_*_*_WWSel2', 
        'keep *_ak5PFJets_*_WWSel2',
        'keep recoMuons_*_*_*',
        'keep recoGsfElectrons_*_*_*',
        'keep *_TriggerResults_*_*'
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring(
        "DiMuon_V0", "DiEle_V0", "DiEleMu_V0",
        "DiMuon_V1a","DiEle_V1a","DiEleMu_V1a",
        "DiMuon_V1b","DiEle_V1b","DiEleMu_V1b",
    ))
)
#process.endDebug = cms.EndPath(process.outDebug)

if False:
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    for P in process.paths_().keys():
        massSearchReplaceAnyInputTag(getattr(process,P), "tcMet", "pfMet")
    process.outDebug.fileName = "selectedReference_pfMet.root"

