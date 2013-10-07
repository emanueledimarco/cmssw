import FWCore.ParameterSet.Config as cms

process = cms.Process("WGstar")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.EventContent.EventContent_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

isMC = False
#isMC = True
#process.GlobalTag.globaltag = 'START42_V13::All' if isMC else 'GR_R_42_V19::All'

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('RMMEFN'))
MC = [
	'file:/hadoop/cms/store/user/mangano/WGstarToLNu2Mu_TuneZ2_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID085_WGstarToLNu2Mu/05d66fb6a9a6a82265cbb7f3f7cbb2e8/WGstarToLNu2Mu_10_1_Ev3.root',
	'file:/hadoop/cms/store/user/mangano/WGstarToLNu2Mu_TuneZ2_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID085_WGstarToLNu2Mu/05d66fb6a9a6a82265cbb7f3f7cbb2e8/WGstarToLNu2Mu_1_1_GXY.root',
	'file:/hadoop/cms/store/user/mangano/WGstarToLNu2Mu_TuneZ2_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID085_WGstarToLNu2Mu/05d66fb6a9a6a82265cbb7f3f7cbb2e8/WGstarToLNu2Mu_2_2_lqU.root',
	'file:/hadoop/cms/store/user/mangano/WGstarToLNu2Mu_TuneZ2_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID085_WGstarToLNu2Mu/05d66fb6a9a6a82265cbb7f3f7cbb2e8/WGstarToLNu2Mu_2_3_hD3.root',
	'file:/hadoop/cms/store/user/mangano/WGstarToLNu2Mu_TuneZ2_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID085_WGstarToLNu2Mu/05d66fb6a9a6a82265cbb7f3f7cbb2e8/WGstarToLNu2Mu_3_1_hsU.root',
]
DATA = [
	'file:/hadoop/cms/store/user/jfernan2/DoubleMu/R42X_S1_V08_S2_V05_S3_V11_ID143a_DoubleMuon2011Bv1a/47b630e2469c9e3550fd039775e94457//DoubleMuon2011Bv1a_109_1_05o.root',
	'file:/hadoop/cms/store/user/jfernan2/DoubleMu/R42X_S1_V08_S2_V05_S3_V11_ID143a_DoubleMuon2011Bv1a/47b630e2469c9e3550fd039775e94457//DoubleMuon2011Bv1a_111_1_qx7.root',
	'file:/hadoop/cms/store/user/jfernan2/DoubleMu/R42X_S1_V08_S2_V05_S3_V11_ID143a_DoubleMuon2011Bv1a/47b630e2469c9e3550fd039775e94457//DoubleMuon2011Bv1a_112_5_1mt.root',
	'file:/hadoop/cms/store/user/jfernan2/DoubleMu/R42X_S1_V08_S2_V05_S3_V11_ID143a_DoubleMuon2011Bv1a/47b630e2469c9e3550fd039775e94457//DoubleMuon2011Bv1a_114_1_BxT.root',
	'file:/hadoop/cms/store/user/jfernan2/DoubleMu/R42X_S1_V08_S2_V05_S3_V11_ID143a_DoubleMuon2011Bv1a/47b630e2469c9e3550fd039775e94457//DoubleMuon2011Bv1a_116_3_bjA.root',
]
process.source.fileNames = MC if isMC else DATA
#process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


from WWAnalysis.AnalysisStep.wwElectrons_cfi import * 
ELE_FULL  =  ELE_BASE+" && "+ELE_ID_LOOSE_SMURFS+" && "+ELE_BDT_ID_SMURF+" && "+ELE_MERGE_ISO+" && "+ELE_MERGE_CONV+" && "+ELE_MERGE_IP

 
from WWAnalysis.AnalysisStep.wwMuons_cfi import *
MU_FULL  = PRESEL_MU+" && "+MUON_KINKID_CUT+" && "+ MUON_MERGE_ISO+"&& "+MUON_MERGE_IP

process.softMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string("pt > 3 && "+
                     "track.isNonnull && "+
                     "track.hitPattern.numberOfValidPixelHits > 0 && track.hitPattern.numberOfValidHits > 10 && "+
                      MUON_MERGE_IP),
)

process.leadElectronFilter = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string("pt > 20 && "+ELE_FULL),
    filter = cms.bool(True),
)
process.leadMuonFilter = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string("pt > 20 && "+MU_FULL),
    filter = cms.bool(True),
)
process.triMuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("softMuons"),
    minNumber = cms.uint32(3)
)
process.diMuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("softMuons"),
    minNumber = cms.uint32(2)
)

process.pfMet20Filter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("pfMet"),
    cut = cms.string("pt > 20"),
    filter = cms.bool(True),
)


# make the actual IsoDeposits
process.pfCharged = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("reducedPFCands"),
    cut = cms.string("charge() != 0"),
)
process.isoDepMuonsFromPF = cms.EDProducer("CandIsoDepositProducer",
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

from WWAnalysis.AnalysisStep.isoAdding_cff import *
process.reBoostedMuons = isoAddedMuons.clone( muonTag = "softMuons" )
process.reBoostedMuons.deposits.append( cms.PSet(
            label  = cms.string("rePFIso"),
            src    = cms.InputTag("isoDepMuonsFromPF"),
            deltaR = cms.double(0.3),
            vetos  = cms.vstring("softMuons:0.01"),
            weight = cms.string('1'),
            skipDefaultVeto = cms.bool(True),
            mode            = cms.string('sum'),
            ) )

process.isoMuons = process.softMuons.clone(src = "reBoostedMuons", cut = "userFloat('rePFIso')/max(10,pt) < 1")

process.triIsoMuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("isoMuons"),
    minNumber = cms.uint32(3)
)
process.diIsoMuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("isoMuons"),
    minNumber = cms.uint32(2)
)


process.LeptonReco = cms.Sequence(
    process.softMuons 
)

process.LeptonReIso = cms.Sequence(
    process.pfCharged *
    process.isoDepMuonsFromPF *
    process.reBoostedMuons *
    process.isoMuons
)


process.Skim_MuDiMuMET_NoIso = cms.Path(
      process.pfMet20Filter +
      process.LeptonReco +
      process.leadMuonFilter + process.triMuonFilter
)
process.Skim_ElDiMuMET_NoIso = cms.Path(
      process.pfMet20Filter +
      process.LeptonReco +
      process.leadElectronFilter + process.diMuonFilter
)
process.Skim_MuDiMuMET_Iso = cms.Path(
      process.pfMet20Filter +
      process.LeptonReco +
      process.leadMuonFilter + process.triMuonFilter +
      process.LeptonReIso +
      process.triIsoMuonFilter
)
process.Skim_ElDiMuMET_Iso = cms.Path(
      process.pfMet20Filter +
      process.LeptonReco +
      process.leadElectronFilter + process.diMuonFilter +
      process.LeptonReIso +
      process.diIsoMuonFilter
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("lmumuSkim.root"),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_*_*_WGstar",
        "keep *_TriggerResults_*_WGstar",
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
)
process.end = cms.EndPath(process.out)
process.out.SelectEvents.SelectEvents += [ 'Skim_MuDiMuMET_Iso' ]
process.out.SelectEvents.SelectEvents += [ 'Skim_ElDiMuMET_Iso' ]
#process.out.SelectEvents.SelectEvents += [ 'Skim_MuDiMuMET_NoIso' ]
#process.out.SelectEvents.SelectEvents += [ 'Skim_ElDiMuMET_NoIso' ]
    
