import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

isMC = True
#isMC = True
#process.GlobalTag.globaltag = 'START42_V13::All' if isMC else 'GR_R_42_V19::All'

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('RMMEFN'))
MC    = [ 'file:latinosYieldSkim.root' ]
#MC    = [ 'file:/data/gpetrucc/7TeV/hzz/step1/HToZZTo4L_M-120_Fall11S6.00215E21D5C4.root' ]
#MC    = [ 'file:/afs/cern.ch/user/g/gpetrucc/scratch0/higgs/CMSSW_4_2_7_patch2/src/WWAnalysis/SkimStep/test/latinosYieldSkim.root' ]
MCAOD = [ 'file:/data/gpetrucc/7TeV/hzz/aod/HToZZTo4L_M-120_Fall11S6.00215E21D5C4.root' ]
DATA  = []
process.source.fileNames = MC if isMC else DATA
#process.source.secondaryFileNames = cms.untracked.vstring(*MCAOD)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

combiner = cms.EDProducer("CandViewMerger", src = cms.VInputTag())
counter  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("FIXME"), minNumber = cms.uint32(1))
filter   = cms.EDFilter("CandViewSelector", src = cms.InputTag("FIXME"), cut = cms.string(""), filter = cms.bool(True))

### Load MC-related sequences:
###   genTauVeto   EDFilter selecting events with one Z decaying to taus (to be applied with a ~ to veto them)
###   gen3RecoSeq  produces a gen3LLLL = 4-lepton candidates out of status == 3 genParticles (no acceptance cuts)
###   gen1RecoSeq  produces a gen1LLLL = 4-lepton candidates out of status == 1 genParticles 
###                      muons pt > 2, |eta| < 2.6; electrons pt > 3, |eta| < 3.2
###   gen3FilterAny EDFilter to select events with one ZZ4L gen candidate at status == 3
###   gen1FilterEta254PtMin53 EDFilter to select one ZZ4L gen candidate, status == 1
##                       muons pt > 3, |eta| < 2.4, electrons pt > 5, |eta| < 2.5
process.load("WWAnalysis.AnalysisStep.zz4l.mcSequences_cff")

### Load sequence to do combinatorics of reco objects
process.load("WWAnalysis.AnalysisStep.zz4l.recoSequence_cff")
process.recMu.cut = "pt > 3 && track.isNonnull && (isGlobalMuon || numberOfMatches > 0 || userInt('pfMuId') > 0)"
process.recEl.cut = "pt > 5 && test_bit(electronID('cicVeryLoose'),0) == 1 && gsfTrack.trackerExpectedHitsInner.numberOfHits<=1"

### Load matcher from gen to reco (module genRecoMatcher)
process.load("WWAnalysis.AnalysisStep.genRecoMatcher_cfi")

### make IsoDeposits (sequence isoDepositSeq)
process.load("WWAnalysis.AnalysisStep.zz4l.isoDeposits_cff")

### make 4L candidates with isolation vetoing other leptons (module isoLLLL)
process.load("WWAnalysis.AnalysisStep.zz4l.pfIsolation_cff")

process.skimEvent4L = cms.EDProducer("SkimEvent4LProducer",
    src = cms.InputTag("isoLLLL"),
    #jets = cms.InputTag("slimPatJetsTriggerMatch"),
    jets = cms.InputTag("slimPatJets"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    isMC = cms.bool(isMC),
    mcMatch = cms.InputTag("genRecoMatcher"),
)

process.load("WWAnalysis.AnalysisStep.zz4lTree_cfi")

if isMC:
    process.begin = cms.Sequence(
        ~process.genTauVeto +
        process.gen3RecoSeq +
        process.gen1RecoSeq +
        process.gen1FilterEta254PtMin53 +
        process.recoSeq +
        process.genRecoMatcher
    )
else:
    process.begin = cms.Sequence(
        process.recoSeq
    )

process.go = cms.Path(
    process.begin +
    process.recFilterAny +
    process.isoDepositSeq +
    process.isoLLLL +
    process.skimEvent4L +
    process.zz4lTree
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("zz4lTree.root"))
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.maxEvents.input = -1
#process.source.fileNames = [ 'file:/data/gpetrucc/7TeV/hzz/step1/DY50JetsToLL.Skim4L.root' ]
