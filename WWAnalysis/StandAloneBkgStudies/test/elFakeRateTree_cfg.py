import FWCore.ParameterSet.Config as cms

process = cms.Process("LepFakeRate")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("WWAnalysis.StandAloneBkgStudies.looseLeptonDefinitions_cff")

process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
#from WWAnalysis.AnalysisStep.wwElectrons_cfi import wwEleIPMerge  as tightElectrons
process.looseElectrons = process.wwEleLoose
process.tightElectrons = process.wwEleIPMerge

from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_BASE
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_ID_LOOSE_2011
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ID
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ISO
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_CONV
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_IP

#from WWAnalysis.AnalysisStep.wwMuons_cfi import PRESEL_MU
#from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT
#from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_ISO
#from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_IP




process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/home/mangano/productionFolder/reloadLP/forStep1/CMSSW_4_1_8_patch1/myTestOnFakes/QCD.root'
        'file:/home/mangano/newDevelopment/CMSSW_4_2_7_patch2/src/WWAnalysis/StandAloneBkgStudies/test/WJetFakeSkim.root'
#        'file:/home/mangano/newDevelopment/CMSSW_4_2_7_patch2/src/WWAnalysis/StandAloneBkgStudies/test/EleFakeSkim.root'
    )
)

process.realJets = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("slimPatJetsTriggerMatch"),
        cut = cms.string("pt > 15 && (chargedMultiplicity+neutralMultiplicity) > 1"),
        filter = cms.bool(True),
)

process.leadingJet = cms.EDFilter("LargestPtCandViewSelector",
    src = cms.InputTag("realJets"), 
    maxNumber = cms.uint32(1)
)

process.leadingEle = cms.EDFilter("LargestPtCandViewSelector",
    src = cms.InputTag("looseElectrons"), 
    maxNumber = cms.uint32(1)
)

#### ele pairs
process.diLeptons4Veto = cms.EDProducer("CandViewShallowCloneCombiner",
#        decay = cms.string('tightElectrons@+ tightElectrons@-'),
        decay = cms.string('boostedElectrons@+ boostedElectrons@-'),
        cut = cms.string('abs(daughter(0).pdgId) == abs(daughter(1).pdgId)')
)
#        cut = cms.string('abs(daughter(0).pdgId) == abs(daughter(1).pdgId) && (mass < 12 || abs(mass-91.1876) < 15)')

process.diLeptons4VetoFilter  = cms.EDFilter("CandViewCountFilter", 
 src = cms.InputTag("diLeptons4Veto"), minNumber = cms.uint32(1)
)





############ loose
process.looseEleBooster = cms.EDProducer("LooseElectronBooster",
  eleTag = cms.InputTag("looseElectrons"),   
  elePairTag = cms.InputTag("diLeptons4Veto"),                                       
  jetTag = cms.InputTag("leadingJet"),  
  leadingElTag = cms.InputTag("leadingEle"),                                            
  metTag = cms.InputTag("pfMet"),
  #triggerTag = cms.InputTag("TriggerResults","","REDIGI311X"),
  triggerTag = cms.InputTag("TriggerResults","","HLT"),
  drMin = cms.double(0.0),
  looseEleSelection = cms.string(ELE_BASE + " && " + ELE_ID_LOOSE_2011),
  tightEleSelection = cms.string(ELE_BASE + " && " + ELE_MERGE_ID + " && " + ELE_MERGE_ISO + " && " + ELE_MERGE_CONV + " && " + ELE_MERGE_IP),
)

process.looseElePlusJet = cms.EDProducer("CandViewCombiner",
        decay = cms.string('leadingJet looseEleBooster'),
        cut = cms.string(''),
        checkCharge = cms.bool(False),
        checkOverlap = cms.bool(False),
)

process.treeL = cms.EDFilter("ProbeTreeProducer",
 src = cms.InputTag("looseElePlusJet"),
 variables = cms.PSet(
   lepPt = cms.string("daughter(1).pt()"),
   lepEta = cms.string("daughter(1).eta()"),
   lepPhi = cms.string("daughter(1).phi()"),
   lepIsEB = cms.string("daughter(1).isEB()"),
   lepIsEE = cms.string("daughter(1).isEE()"),
   jetPt = cms.string("daughter(0).pt()"),
   dR = cms.string("deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi)"),
   mt = cms.string("daughter(1).userFloat('mt')"),
   pfmet = cms.string("daughter(1).userFloat('pfmet')"),
   mass = cms.string("daughter(1).userFloat('pairMass')"),
   lLepPt = cms.string("daughter(1).userFloat('leadingLepPt')"),
 ),
 flags = cms.PSet(
   isTight = cms.string("daughter(1).userFloat('passTight')"),
   isLoose = cms.string("daughter(1).userFloat('passLoose')"),
   passEl8 = cms.string("daughter(1).userInt('passEl8')"),
   passEl17 = cms.string("daughter(1).userInt('passEl17')"),
   ),
)





process.TFileService = cms.Service("TFileService",fileName = cms.string("tree.root"))

process.out = cms.OutputModule("PoolOutputModule", 
                               outputCommands =  cms.untracked.vstring(), 
                               fileName = cms.untracked.string('output.root') 
                               )

process.out.outputCommands =  cms.untracked.vstring(
    'drop *',
    #'keep *_slimPatJetsTriggerMatch_*_*',
    #'keep *_boostedEles_*_*',
    #'keep *_onlyIDEles_*_Demo',
    #'keep *_tightEles_*_Demo',
    #'keep *_looseEleBooster_*_*',
    #'keep *_myFirstNtuple_*_*'
    'keep *_diLeptons4Veto_*_*'
)


# ===== TRIGGER FILTER
from HLTrigger.HLTfilters.triggerResultsFilter_cfi import triggerResultsFilter
process.hltFilter = triggerResultsFilter.clone(
    l1tResults = '',
    #hltResults = cms.InputTag( "TriggerResults", "", "REDIGI311X"),
    hltResults = cms.InputTag( "TriggerResults", "", "HLT"),
    throw = True,
    triggerConditions = [ 'HLT_Ele8_CaloIdL_CaloIsoVL_v*',
                          'HLT_Ele17_CaloIdL_CaloIsoVL_v*',
                          #'HLT_Ele8_v*',
                          ]
    )


# ======== For Filtering muons in MuEnriched QCD sample
process.genMuFilter = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("prunedGen"),
        cut = cms.string("abs(pdgId) == 13 && pt >15 && abs(eta)<2.5 && status == 1"),  
        filter = cms.bool(True),
)



# ==== TRIGGER REPORT
process.MessageLogger.categories.append('HLTrigReport')
process.load( "HLTrigger.HLTanalyzers.hlTrigReport_cfi" )
process.hlTrigReport.HLTriggerResults   = cms.InputTag("TriggerResults", "", "REDIGI311X")
process.hlTrigReport.ReferencePath      = cms.untracked.string( "HLTriggerFinalPath" )
process.hlTrigReport.ReferenceRate      = cms.untracked.double( 100.0 )






# ===== PATHS
process.p = cms.Path(
    #process.hltFilter*
    #~process.genMuFilter*
    process.realJets*
    process.tightElectrons*process.looseElectrons*
    process.diLeptons4Veto*#~process.diLeptons4VetoFilter*
    process.leadingJet*process.leadingEle*
    #process.tightEles*
    process.looseEleBooster*process.looseElePlusJet*process.treeL
)




#process.report = cms.EndPath( process.hlTrigReport )
process.outpath    = cms.EndPath(process.out)


#Options
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

