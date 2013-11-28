import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
what = args[0] if len(args) > 0 else "Data2011B"

process = cms.Process("WGstarStep2")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 25

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


from WWAnalysis.AnalysisStep.wwElectrons_cfi import * 
ELE_FULL  =  ELE_BASE+" && "+ELE_ID_LOOSE_SMURFS+" && "+ELE_BDT_ID_SMURF+" && "+ELE_MERGE_ISO+" && "+ELE_MERGE_CONV+" && "+ELE_MERGE_IP
ELE_NOISO =  ELE_BASE+" && "+ELE_ID_LOOSE_SMURFS+" && "+ELE_BDT_ID_SMURF+" && "+                     ELE_MERGE_CONV+" && "+ELE_MERGE_IP

 
from WWAnalysis.AnalysisStep.wwMuons_cfi import *
MU_FULL  = PRESEL_MU+" && "+MUON_KINKID_CUT+" && "+ MUON_MERGE_ISO+"&& "+MUON_MERGE_IP
MU_NOISO = PRESEL_MU+" && "+MUON_KINKID_CUT+                       "&& "+MUON_MERGE_IP

merger  = cms.EDProducer("CandViewMerger",    src = cms.VInputTag())
counter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("FILL_ME"), minNumber = cms.uint32(1))

process.softMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string("pt > 5 && "+
                     "track.isNonnull && track.hitPattern.numberOfValidPixelHits > 0 && track.hitPattern.numberOfValidHits > 11 && "+
                     MUON_MERGE_IP+" && "+
                     "(numberOfMatches > 1 || muonID('TMOneStationTight'))"),
)

# make the actual IsoDeposits
process.pfCharged = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("reducedPFCands"),
    cut = cms.string("charge() != 0"),
)
process.pfNeutral = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("reducedPFCands"),
    cut = cms.string("charge() == 0"),
)
process.isoDepMuonsFromPFCharged = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("softMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        inputCandView = cms.InputTag("pfCharged"),
        Diff_z = cms.double(999.9),
        DR_Max = cms.double(0.5),
        Diff_r = cms.double(99999.99),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)
process.isoDepMuonsFromPFNeutral = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("softMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        inputCandView = cms.InputTag("pfNeutral"),
        Diff_z = cms.double(999.9),
        DR_Max = cms.double(0.5),
        Diff_r = cms.double(99999.99),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)

from WWAnalysis.AnalysisStep.isoAdding_cff import *
process.reBoostedMuons = isoAddedMuons.clone( muonTag = "softMuons" )
process.reBoostedMuons.deposits.append( cms.PSet(
            label  = cms.string("rePFIso"),
            src    = cms.InputTag("isoDepMuonsFromPFCharged"),
            deltaR = cms.double(0.3),
            vetos  = cms.vstring("softMuons:0.01"),
            weight = cms.string('1'),
            skipDefaultVeto = cms.bool(True),
            mode            = cms.string('sum'),
            ) )
process.reBoostedMuons.deposits.append( cms.PSet(
            label  = cms.string("rePFIsoN"),
            src    = cms.InputTag("isoDepMuonsFromPFNeutral"),
            deltaR = cms.double(0.3),
            vetos  = cms.vstring("Threshold(1.0)"),
            weight = cms.string('1'),
            skipDefaultVeto = cms.bool(True),
            mode            = cms.string('sum'),
            ) )

#process.isoMuons = process.softMuons.clone(src = "reBoostedMuons", cut = "pt > 3 && "+MU_NOISO +" && userFloat('rePFIso')/pt < 0.2")
process.isoMuons = process.softMuons.clone(src = "reBoostedMuons", cut = "pt > 5")
process.isoElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string("pt > 5 && "+ELE_FULL)
)

process.isoElectronsForW = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string("pt > 15 && "+ELE_FULL)
)
process.isoMuonsForW = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("isoMuons"),
    cut = cms.string("pt > 15 && "+MU_FULL),
)

process.isoLeptons = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("isoMuons"), cms.InputTag("isoElectrons"))
)
process.isoLeptonsForW = merger.clone(src = [cms.InputTag("isoMuonsForW"), cms.InputTag("isoElectronsForW")])

process.leadingMuonFilter =  cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("isoMuonsForW"),
    minNumber = cms.uint32(1)
)
process.leadingElectronFilter = counter.clone(src = "isoElectronsForW")
process.tripleIsoMuFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("isoMuons"),
    minNumber = cms.uint32(3)
)
process.doubleIsoMuFilter = process.tripleIsoMuFilter.clone(minNumber = 2)

process.gammaStarMM = cms.EDProducer("CandViewShallowCloneCombiner",
    #decay = cms.string("isoMuons@+ isoMuons@-"),
    decay = cms.string("isoMuons isoMuons"),
    cut   = cms.string("mass < 200 && "+
                       #"daughter(0).pt >  8 && "+
                       "daughter(0).pt >= daughter(1).pt && "+
                       "deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.01"),
    checkCharge = cms.bool(False),
)
#process.gammaStarEE = process.gammaStarMM.clone(decay = "isoElectrons@+ isoElectrons@-")
#process.gammaStarEE = process.gammaStarMM.clone(decay = "isoElectrons isoElectrons")
#process.gammaStarLL = cms.EDProducer("CandViewMerger",
#    src = cms.VInputTag(cms.InputTag("gammaStarEE"), cms.InputTag("gammaStarMM"))
#)
#process.gammaStarLLFilter = counter.clone(src = "gammaStarMM")
process.gammaStarMMFilter = counter.clone(src = "gammaStarMM")
#process.gammaStarEEFilter = process.gammaStarLLFilter.clone(src = "gammaStarEE")

MT="sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi)))"
process.wToMuNu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("isoMuonsForW@+ pfMet"),
    cut   = cms.string("30 < "+MT+" < 120 && daughter(0).pt > 15"),
)
process.wToElNu = process.wToMuNu.clone(decay = "isoElectronsForW@+ pfMet")
process.wToLNu = merger.clone(src = [cms.InputTag("wToMuNu"), cms.InputTag("wToElNu")])
process.wToLNuFilter = counter.clone(src = "wToLNu")

NO_OVERLAP = "deltaR(daughter(0).daughter(0).eta,daughter(0).daughter(0).phi,daughter(1).daughter(%d).eta,daughter(1).daughter(%d).phi) > 0.02"
process.wgStar = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("wToLNu gammaStarMM"),
    cut   = cms.string("%s && %s" % (NO_OVERLAP % (0,0) , NO_OVERLAP % (1,1))),
    checkCharge = cms.bool(False),
)
process.wgStarFilter = counter.clone(src = "wgStar")

process.LeptonReco = cms.Sequence(
    ( process.softMuons + 
      process.pfCharged +
      process.pfNeutral ) *
    ( process.isoDepMuonsFromPFCharged *
      process.isoDepMuonsFromPFNeutral *
      process.reBoostedMuons *
      process.isoMuons +
      process.isoElectrons ) *
    ( process.isoMuonsForW + process.isoElectronsForW ) *
    ( process.isoLeptonsForW +
      process.isoLeptons )
)
process.CombinatorialReco = cms.Sequence(
    #( process.gammaStarMM + process.gammaStarEE ) *
    #  process.gammaStarLL *
      process.gammaStarMM *
    ( process.wToMuNu + process.wToElNu ) *
      process.wToLNu +
      process.wgStar 
)

process.reco_WGstar_filters = cms.Sequence(
      process.LeptonReco +
      #process.tripleIsoLeptonFilter + 
      #process.leadingIsoLeptonFilter +
      process.CombinatorialReco +
      process.gammaStarMMFilter +
      process.wToLNuFilter +
      process.wgStarFilter 
)

process.jetsNotLeptons = cms.EDProducer("CandViewCleaner",
    srcObject = cms.InputTag("slimPatJetsTriggerMatch"),
    srcObjectsToRemove = cms.VInputTag("isoLeptonsForW"),
    deltaRMin = cms.double(0.3),
    srcObjectSelection = cms.string("pt > 30"),
    srcObjectsToRemoveSelection = cms.string("pt > 20"), # don't remove jets which overlap with soft leptons
)
process.njets = cms.EDProducer("CandMultiplicityCounter", 
    probes  = cms.InputTag("wgStar"),
    objects = cms.InputTag("jetsNotLeptons"),
)
process.nbjets = cms.EDProducer("CandMultiplicityCounter",
    probes  = cms.InputTag("wgStar"),
    objects = cms.InputTag("slimPatJetsTriggerMatch"),
    objectSelection = cms.string("pt > 10 && bDiscriminator('trackCountingHighEffBJetTags') > 2.1 && abs(userFloat('dz')) < 2.0 && abs(userFloat('dzPV')) < 0.1"),
)

process.nleps = cms.EDProducer("CandMultiplicityCounter",
    probes  = cms.InputTag("wgStar"),
    objects = cms.InputTag("isoLeptons"),
)
process.ncands = cms.EDProducer("CandMultiplicityCounter",
    probes  = cms.InputTag("wgStar"),
    objects = cms.InputTag("wgStar"),
)
process.ncandsOS = cms.EDProducer("CandMultiplicityCounter",
    probes  = cms.InputTag("wgStar"),
    objects = cms.InputTag("wgStar"),
    objectSelection = cms.string("daughter(1).daughter(0).charge != daughter(1).daughter(1).charge")
)
process.ncandsMOS = cms.EDProducer("CandMultiplicityCounter",
    probes  = cms.InputTag("wgStar"),
    objects = cms.InputTag("wgStar"),
    objectSelection = cms.string("daughter(1).daughter(0).charge != daughter(1).daughter(1).charge && abs(daughter(0).daughter(0).pdgId) == 13")
)


process.chmet = cms.EDProducer("ChargedMETComputerForTnP",
    objects = cms.InputTag("wgStar"),
    pfCands = cms.InputTag("reducedPFCands"),
    pfCandSelection = cms.string("charge() != 0"),
)

process.wgTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("wgStar"),
    cut = cms.string(""),
    variables = cms.PSet(
        # compound object variables
        wg_pt = cms.string("pt"),
        wg_m  = cms.string("mass"),
        wg_mt = cms.string("sqrt(max(0, pow(hypot(daughter(1).pt, daughter(1).mass) + hypot(daughter(0).pt, 80.4), 2) - pow(daughter(1).px+daughter(0).px,2) - pow(daughter(1).py+daughter(0).py,2) ))"),
        w_pt  = cms.string("daughter(0).pt"),
        w_mt  = cms.string("sqrt(2*daughter(0).daughter(0).pt*daughter(0).daughter(1).pt*(1 - cos(daughter(0).daughter(0).phi - daughter(0).daughter(1).phi)))"),
        z_pt  = cms.string("daughter(1).pt"),
        z_m   = cms.string("daughter(1).mass"),
        z_y   = cms.string("daughter(1).rapidity"),
        # lepton kinematic variables
        l1_pt  = cms.string("daughter(0).daughter(0).pt"),
        l1_eta = cms.string("daughter(0).daughter(0).eta"),
        l1_phi = cms.string("daughter(0).daughter(0).phi"),
        l1_q   = cms.string("daughter(0).daughter(0).charge"),
        l1_id  = cms.string("daughter(0).daughter(0).pdgId"),
        l2_pt  = cms.string("daughter(1).daughter(0).pt"),
        l2_eta = cms.string("daughter(1).daughter(0).eta"),
        l2_phi = cms.string("daughter(1).daughter(0).phi"),
        l2_q   = cms.string("daughter(1).daughter(0).charge"),
        l2_id  = cms.string("daughter(1).daughter(0).pdgId"),
        l3_pt  = cms.string("daughter(1).daughter(1).pt"),
        l3_eta = cms.string("daughter(1).daughter(1).eta"),
        l3_phi = cms.string("daughter(1).daughter(1).phi"),
        l3_q   = cms.string("daughter(1).daughter(1).charge"),
        l3_id  = cms.string("daughter(1).daughter(1).pdgId"),
        met    = cms.string("daughter(0).daughter(1).pt"),
        # other variables
        rel_iso_wE = cms.string("daughter(0).daughter(0).masterClone.userFloat('eleSmurfPF04')/daughter(0).daughter(0).pt"),
        rel_iso_wM = cms.string("daughter(0).daughter(0).masterClone.userFloat('muSmurfPF' )/daughter(0).daughter(0).pt"),
        worse_rel_iso_z  = cms.string("max(daughter(1).daughter(0).masterClone.userFloat('rePFIso')/daughter(1).daughter(0).pt, "+
                                          "daughter(1).daughter(1).masterClone.userFloat('rePFIso')/daughter(1).daughter(1).pt)"),
        worse_rel_isoN_z = cms.string("max(daughter(1).daughter(0).masterClone.userFloat('rePFIsoN')/daughter(1).daughter(0).pt, "+
                                          "daughter(1).daughter(1).masterClone.userFloat('rePFIsoN')/daughter(1).daughter(1).pt)"),
        worse_rel_isoA_z = cms.string("max(daughter(1).daughter(0).masterClone.userFloat('rePFIsoN')/daughter(1).daughter(0).pt + "+
                                          "daughter(1).daughter(0).masterClone.userFloat('rePFIso' )/daughter(1).daughter(0).pt,  "+
                                          "daughter(1).daughter(1).masterClone.userFloat('rePFIsoN')/daughter(1).daughter(1).pt + "+
                                          "daughter(1).daughter(1).masterClone.userFloat('rePFIso' )/daughter(1).daughter(1).pt)"),
        sum_iso_z  = cms.string("daughter(1).daughter(0).masterClone.userFloat('rePFIso') + "+
                                "daughter(1).daughter(1).masterClone.userFloat('rePFIso')"),
        sum_isoN_z = cms.string("daughter(1).daughter(0).masterClone.userFloat('rePFIsoN') + "+
                                "daughter(1).daughter(1).masterClone.userFloat('rePFIsoN')"),
        worse_rel_iso_3 = cms.string("max(max(daughter(0).daughter(0).masterClone.userFloat('muSmurfPF')/daughter(0).daughter(0).pt,  "+
                                             "daughter(0).daughter(0).masterClone.userFloat('eleSmurfPF04')/daughter(0).daughter(0).pt),"+
                                         "max(daughter(1).daughter(0).masterClone.userFloat('rePFIso')/daughter(1).daughter(0).pt, "+
                                             "daughter(1).daughter(1).masterClone.userFloat('rePFIso')/daughter(1).daughter(1).pt))"),
        # lepton id variables
        #l2_muid_soft = cms.string("daughter(1).daughter(0).masterClone.muonID('TMOneStationTight')"),
        l2_muid_tm2  = cms.string("daughter(1).daughter(0).masterClone.numberOfMatches"),
        #l2_muid_tkh  = cms.string("daughter(1).daughter(0).masterClone.track.hitPattern.numberOfValidHits"),
        #l3_muid_soft = cms.string("daughter(1).daughter(1).masterClone.muonID('TMOneStationTight')"),
        l3_muid_tm2  = cms.string("daughter(1).daughter(1).masterClone.numberOfMatches"),
        #l3_muid_tkh  = cms.string("daughter(1).daughter(1).masterClone.track.hitPattern.numberOfValidHits"),
        l2_muid_tip = cms.string("abs(daughter(1).daughter(0).masterClone.userFloat('tip'))"),
        l3_muid_tip = cms.string("abs(daughter(1).daughter(1).masterClone.userFloat('tip'))"),
        l1_muid_lip = cms.string("abs(daughter(1).daughter(0).masterClone.userFloat('dzPV'))"),
        l2_muid_lip = cms.string("abs(daughter(1).daughter(0).masterClone.userFloat('dzPV'))"),
        l3_muid_lip = cms.string("abs(daughter(1).daughter(1).masterClone.userFloat('dzPV'))"),
        # external variables
        # external variables
        njets  = cms.InputTag("njets"),
        nleps  = cms.InputTag("nleps"),
        nbjets = cms.InputTag("nbjets"),
        ncands = cms.InputTag("ncands"),
        ncandsOS  = cms.InputTag("ncandsOS"),
        ncandsMOS = cms.InputTag("ncandsMOS"),
        chmet  = cms.InputTag("chmet"),
        # scale factor (for MC)
        scalef = cms.string("1"),
    ),
    flags = cms.PSet(
        OS  = cms.string("daughter(1).daughter(0).charge() != daughter(1).daughter(1).charge()"),
        SS  = cms.string("daughter(1).daughter(0).charge() == daughter(1).daughter(1).charge()"),
        EMM = cms.string("abs(daughter(0).daughter(0).pdgId) == 11 && abs(daughter(1).daughter(0).pdgId) == 13"),
        MMM = cms.string("abs(daughter(0).daughter(0).pdgId) == 13 && abs(daughter(1).daughter(0).pdgId) == 13"),
        #EEE = cms.string("abs(daughter(0).daughter(0).pdgId) == 11 && abs(daughter(1).daughter(0).pdgId) == 11"),
        #MEE = cms.string("abs(daughter(0).daughter(0).pdgId) == 13 && abs(daughter(1).daughter(0).pdgId) == 11"),
    ),
    addRunLumiInfo = cms.bool(True)
)

process.treeVariables = cms.Sequence(
    process.jetsNotLeptons +
    process.njets +
    process.nbjets +
    process.nleps  +
    process.ncands +
    process.ncandsOS +
    process.ncandsMOS +
    process.chmet
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("prunedGen"),                                                                 
    printP4 = cms.untracked.bool(True),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex = cms.untracked.bool(True),
    status = cms.untracked.vint32( 3 )
)
process.printList = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),
  printVertex = cms.untracked.bool(False),
  src = cms.InputTag("prunedGen")
)


process.Skim_WGstar = cms.Path(
    process.reco_WGstar_filters
)
process.Analysis_WGstar = cms.Path(
    process.reco_WGstar_filters +
    process.treeVariables       +
    process.wgTree
)
#process.Debug_WGstar = cms.Path(
#    process.reco_WGstar_filters +
#    process.printTree +
#    process.printList
#)
#process.maxEvents.input = 50


process.TFileService = cms.Service("TFileService", fileName = cms.string("wgStarPlots.root"))
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("wgStar.root"),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_*_*_WGstar",
        "keep *_reBoostedElectrons_*_*",
        "keep *_reBoostedMuons_*_*",
        "keep *_isoElectrons_*_*",
        "keep *_isoMuons_*_*",
        "keep recoCompositeCandidates_*_*_WGstar",
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("Skim_WGstar")),
)
#process.end = cms.EndPath(process.out)


scalef = {
    'WGstarMM'  : 0.00914,
    'Wgamma'    : 0.05182,
    'Wjets'     : 0.38636,
    'WZ'        : 0.00149,
    'ZZ'        : 0.00449,
    'TT'        : 0.01145,
    'TW'        : 0.00969,
    'TbarW'     : 0.02439,
}
process.source.fileNames = [ 'file:lmumuSkim_%s.root' % what ]
if "2011" in what:
    from glob import glob
    process.source.fileNames = [ 'file:'+x for x in glob("lmumuSkim_"+what+"*.root") ]

process.out.fileName = "wgStar_%s.root" % what
process.TFileService.fileName = "wgStarPlots_%s.root" % what

if "2011" not in what: process.wgTree.variables.scalef = str(4.63 * scalef[what])
