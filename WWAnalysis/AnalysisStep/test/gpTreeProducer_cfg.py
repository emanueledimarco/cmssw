import FWCore.ParameterSet.Config as cms

process = cms.Process("EWKTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
#     fileNames = cms.untracked.vstring('file:RMMEIN.root'),
    fileNames = cms.untracked.vstring(
#        'file:/hadoop/cms/store/user/mwlebour/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WW_414_SKIM_V04_ID026_WJetsToLNuMad_Step2/a3ede9ee43e9e4198483e6ce9f0ce8dc/WJetsToLNuMad_16_1_UwO.root'
        '/store/user/mangano/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/R414_S1_V05_S2_V02_S3_V00_ID101160_ggToH160toWWto2L2Nu_Step2-v2/2ff183dd1f3dc60c128420a6b73e4e4b/ggToH160toWWto2L2Nu_2_1_hFM.root',
        '/store/user/mangano/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/R414_S1_V05_S2_V02_S3_V00_ID101160_ggToH160toWWto2L2Nu_Step2-v2/2ff183dd1f3dc60c128420a6b73e4e4b/ggToH160toWWto2L2Nu_1_1_Dfc.root',
    ),
)
process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


Template = cms.EDFilter("ProbeTreeProducer",
    cut = cms.string("leptEtaCut() && ptByPt(0) > 20 && ptByPt(1) > 10"),
    variables = cms.PSet(
        hypo = cms.string("hypo()"),
        mll  = cms.string("mll()"),
        ptll = cms.string("pTll()"),
        pt1  = cms.string("ptByPt(0)"),
        pt2  = cms.string("ptByPt(1)"),
        met  = cms.string("met('PFMET')"),
        mmet = cms.string("met('MINMET')"),
        pmet = cms.string("projMet('PFMET')"),
        pmmet = cms.string("projMet('MINMET')"),
        dphill = cms.string("dPhill()/3.14159"),
        drll   = cms.string("dRll()"),
        dphilmet = cms.string("dPhilMet('PFMET')"),
        mtw1 = cms.string("mTByPt(0,'PFMET')"),
        mtw2 = cms.string("mTByPt(1,'PFMET')"),
        mth  = cms.string("mTHiggs('PFMET')"),
        njet  = cms.string("nJets(30)"),
        cnjet = cms.string("nCentralJets(30,2.5)"),
        bestiso  = cms.string("allIsoByIso(0)/ptByIso(0)"),
        worstiso = cms.string("allIsoByIso(1)/ptByIso(1)"),
        worstJetLepPt = cms.string("max(matchedJetPt(0, 0.5)/pt(0), matchedJetPt(1, 0.5)/pt(1))"),
    ),
    flags = cms.PSet(
        bveto  = cms.string("bTaggedJetsUnder(30,2.1) == 0 && nSoftMu(3) == 0"),
        lveto = cms.string("nExtraLep == 0"),
        zveto = cms.string("abs(mll-91.1876)>15"),
        trigger = cms.string("triggerMatchingCut('MC')"),
    ),
)

process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter.triggerConditions = cms.vstring(  )
process.triggerResultsFilter.l1tResults = ''
process.triggerResultsFilter.throw = True
process.triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "REDIGI311X" )
process.elelTrigger = process.triggerResultsFilter.clone(triggerConditions = [ "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*"])
process.elmuTrigger = process.triggerResultsFilter.clone(triggerConditions = [ "HLT_Mu21_v* OR HLT_Mu5_Ele17_v* OR HLT_Mu11_Ele8_v*"])
process.muelTrigger = process.elmuTrigger.clone()
process.mumuTrigger = process.triggerResultsFilter.clone(triggerConditions = [ "HLT_Mu21_v* OR HLT_DoubleMu5_v*"])

for X in "elel", "mumu", "elmu", "muel":
    setattr(process, X+"Tree",     Template.clone(src = cms.InputTag("ww%sIPLHT"   % X)))
    setattr(process, X+"TreeNoPU", Template.clone(src = cms.InputTag("ww%sJetNoPU" % X)))
    setattr(process, X+"Path",     cms.Path(getattr(process,X+"Trigger") + getattr(process,X+"Tree")))# + getattr(process,X+"TreeNoPU")))

# process.TFileService = cms.Service("TFileService",fileName = cms.string("RMMEOUT"))
process.TFileService = cms.Service("TFileService",fileName = cms.string("tree.root"))
 
