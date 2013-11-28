import FWCore.ParameterSet.Config as cms

process = cms.Process("MuFakeRate")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/mangano/productionFolder/reloadLP/forStep1/CMSSW_4_1_8_patch1/myTestOnFakes/QCD.root'
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

process.leadingMu = cms.EDFilter("LargestPtCandViewSelector",
    src = cms.InputTag("onlyIDMuons"), 
    maxNumber = cms.uint32(1)
)

#### muon pairs
process.diLeptons4Veto = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string('IdAndIpMuons@+ IdAndIpMuons@-'),
        cut = cms.string('abs(daughter(0).pdgId) == abs(daughter(1).pdgId) && (mass < 12 || abs(mass-91.1876) < 15)')
)

process.diLeptons4VetoFilter  = cms.EDFilter("CandViewCountFilter", 
 src = cms.InputTag("diLeptons4Veto"), minNumber = cms.uint32(1)
)


process.load("WWAnalysis.StandAloneBkgStudies.looseLeptonDefinitions_cff")



############ loose
process.looseMuonBooster = cms.EDProducer("LooseMuonBooster",
  muonTag = cms.InputTag("onlyIDMuons"),   
  muonPairTag = cms.InputTag("diLeptons4Veto"),                                       
  jetTag = cms.InputTag("leadingJet"),  
  leadingMuTag = cms.InputTag("leadingMu"),                                            
  metTag = cms.InputTag("pfMet"),
  #triggerTag = cms.InputTag("TriggerResults","","REDIGI311X"),
  triggerTag = cms.InputTag("TriggerResults","","HLT"),
  drMin = cms.double(0.0)
)

process.looseMuonPlusJet = cms.EDProducer("CandViewCombiner",
        decay = cms.string('leadingJet looseMuonBooster'),
        cut = cms.string(''),
        checkCharge = cms.bool(False),
        checkOverlap = cms.bool(False),
)

process.treeL = cms.EDFilter("ProbeTreeProducer",
 src = cms.InputTag("looseMuonPlusJet"),
 variables = cms.PSet(
   muPt = cms.string("daughter(1).pt()"),
   muEta = cms.string("daughter(1).eta()"),
   muPhi = cms.string("daughter(1).phi()"),
   jetPt = cms.string("daughter(0).pt()"),
   dR = cms.string("deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi)"),
   mt = cms.string("daughter(1).userFloat('mt')"),
   pfmet = cms.string("daughter(1).userFloat('pfmet')"),
   mass = cms.string("daughter(1).userFloat('pairMass')"),
   lMuPt = cms.string("daughter(1).userFloat('leadingMuPt')"),
 ),
 flags = cms.PSet(
   isTight = cms.string("daughter(1).userFloat('passTight')"),
   isLoose = cms.string("daughter(1).userFloat('passLoose')"),
   passMu7 = cms.string("daughter(1).userInt('passMu7')"),
   passMu8 = cms.string("daughter(1).userInt('passMu8')"),
   passMu15 = cms.string("daughter(1).userInt('passMu15')"),
   passMu24 = cms.string("daughter(1).userInt('passMu24')"),
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
    #'keep *_boostedMuons_*_*',
    #'keep *_onlyIDMuons_*_Demo',
    #'keep *_tightMuons_*_Demo',
    #'keep *_looseMuonBooster_*_*',
    #'keep *_myFirstNtuple_*_*'
    'keep *_diLeptons4Veto_*_*'
)



# ===== TRIGGER FILTER
from HLTrigger.HLTfilters.triggerResultsFilter_cfi import triggerResultsFilter
process.hltFilter = triggerResultsFilter.clone(
    l1tResults = '',
    hltResults = cms.InputTag( "TriggerResults", "", "REDIGI311X"),
    throw = True,
    triggerConditions = [ #'HLT_Mu7',
                          'HLT_Mu15_v*',
                          ]
    )


# ======== For Filtering muons in MuEnriched QCD sample
process.genMuFilter = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("prunedGen"),
        cut = cms.string("abs(pdgId) == 13 && pt >15 && abs(eta)<2.5 && status == 1"),  
        filter = cms.bool(True),
)



# ==== TRIGGER REPORT
#process.MessageLogger.categories.append('HLTrigReport')
#process.load( "HLTrigger.HLTanalyzers.hlTrigReport_cfi" )
#process.hlTrigReport.HLTriggerResults   = cms.InputTag("TriggerResults", "", "REDIGI311X")
#process.hlTrigReport.ReferencePath      = cms.untracked.string( "HLTriggerFinalPath" )
#process.hlTrigReport.ReferenceRate      = cms.untracked.double( 100.0 )






# ===== PATHS
process.p = cms.Path(
    #process.hltFilter*
    #~process.genMuFilter*
    process.realJets*
    process.IdAndIpMuons*process.diLeptons4Veto*#~process.diLeptons4VetoFilter*
    process.onlyIDMuons*
    process.leadingJet*process.leadingMu*
    #process.tightMuons*
    process.looseMuonBooster*process.looseMuonPlusJet*process.treeL
)




#process.report = cms.EndPath( process.hlTrigReport )
#process.outpath    = cms.EndPath(process.out)


#Options
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

