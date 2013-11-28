import FWCore.ParameterSet.Config as cms

process = cms.Process("Skim4L")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.EventContent.EventContent_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

isMC = True
#isMC = True
#process.GlobalTag.globaltag = 'START42_V13::All' if isMC else 'GR_R_42_V19::All'

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('RMMEFN'))
MC = [
	'/store/user/mwlebour/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID037_DY50toLLMadD6T/26a07df0670d12c6eef310bf324d6d6d/DY50toLLMadD6T_62_0_TQq.root',
	'/store/user/mwlebour/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID037_DY50toLLMadD6T/26a07df0670d12c6eef310bf324d6d6d/DY50toLLMadD6T_63_0_9QG.root',
	'/store/user/mwlebour/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID037_DY50toLLMadD6T/26a07df0670d12c6eef310bf324d6d6d/DY50toLLMadD6T_64_0_rU0.root',
	'/store/user/mwlebour/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID037_DY50toLLMadD6T/26a07df0670d12c6eef310bf324d6d6d/DY50toLLMadD6T_65_0_LvQ.root',
	'/store/user/mwlebour/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/R42X_S1_V09_S2_V06_S3_V13_ID037_DY50toLLMadD6T/26a07df0670d12c6eef310bf324d6d6d/DY50toLLMadD6T_66_0_0NE.root',
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

combiner = cms.EDProducer("CandViewMerger", src = cms.VInputTag())
counter  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("FIXME"), minNumber = cms.uint32(1))
filter   = cms.EDFilter("CandViewSelector", src = cms.InputTag("FIXME"), cut = cms.string(""), filter = cms.bool(True))

process.recMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string("pt > 2 && track.isNonnull && (isGlobalMuon || numberOfMatches > 0) && abs(userFloat('dzPV')) < 0.1"),
)
process.recEl = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string("pt > 3 && abs(userFloat('dzPV')) < 0.1"),
)
process.recMM = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("recMu@+ recMu@-"),
    cut = cms.string("mass > 2")
)
process.recEE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("recEl@+ recEl@-"),
    cut = cms.string("mass > 2")
)
process.recLL = combiner.clone(src = ["recMM", "recEE"])
process.recLLLL = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("recLL recLL"),
    cut = cms.string("daughter(0).daughter(0).pt != daughter(1).daughter(0).pt && daughter(0).daughter(1).pt != daughter(1).daughter(1).pt && "+
                     "daughter(0).daughter(1).pt != daughter(1).daughter(0).pt && daughter(0).daughter(0).pt != daughter(1).daughter(1).pt")
)
process.recRecoSeq = cms.Sequence(process.recMu + process.recEl + process.recMM + process.recEE + process.recLL + process.recLLLL)
process.recAny = counter.clone(src = 'recLLLL')

process.Skim_Rec4L = cms.Path(
    process.recRecoSeq +
    process.recAny 
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("fourEllSkim.root"),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_*_*_Skim4L",
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
)
process.end = cms.EndPath(process.out)
process.out.SelectEvents.SelectEvents += [ 'Skim_Rec4L' ]
    
