import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from glob import glob
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames += ['file:%s'%x for x in glob('/home/avartak/CMS/Higgs/CMSSW_4_2_8/src/WWAnalysis/SkimStep/skims/gghzz4l130/*.root')]
#process.source.fileNames = [ 
#    'root://pcmssd12//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/DoubleElectron2012/hzz4lSkim_42_1_kGb.root',
#    'root://pcmssd12//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/DoubleMu2012/hzz4lSkim_37_1_htl.root'
#]
#process.source.fileNames = [ 
#    'root://pcmssd12//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/ggH_120_52X/hzz4lSkim_101_2_OyG.root',
#    'root://pcmssd12//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/ggH_120_52X/hzz4lSkim_11_2_c7l.root',
#]
process.source.fileNames = [
    #'file:/afs/cern.ch/work/g/gpetrucc/HZZ/CMSSW_5_2_4_patch4/src/WWAnalysis/SkimStep/test/hzz4lSkim.MC_Synch.root'
    'file:/afs/cern.ch/work/g/gpetrucc/CMSSW_4_2_8_patch7/src/WWAnalysis/SkimStep/test/hzz4lSkim.4Sync_NoSmear.root'
]

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("WWAnalysis.AnalysisStep.hzz4l_selection_cff")
process.load("WWAnalysis.AnalysisStep.zz4l.reSkim_cff")
process.load("WWAnalysis.AnalysisStep.zz4l.mcSequences_cff")
process.load("WWAnalysis.AnalysisStep.fourLeptonBlinder_cfi")
process.load("WWAnalysis.AnalysisStep.zz4lTree_cfi")

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *

### Alexey sync customizations
process.skim40.cut = 'mass > 40 && abs(daughter(0).pdgId) == abs(daughter(1).pdgId)'
SEL_ZZ4L_PRE_BEST_1 = "lByPt(0).pt > 20 && lByPt(1).pt > 10 && lByPt(0).userFloat('pfCombRelIso04EACorr') < 0.2 && lByPt(1).userFloat('pfCombRelIso04EACorr') < 0.2"
SEL_ZZ4L_PRE_BEST_2 = "40 < mz(0) < 120 && 4 < mz(1) < 120"
SEL_ZZ4L_PRE_BEST_3 = "nGoodPairs(\"mass > 4\", 1) >= 6"
SEL_ZZ4L_PRE_BEST_4 = "mass > 100"
SEL_ANY_Z = "mass > 0 && daughter(0).pdgId = - daughter(1).pdgId" 


### =========== BEGIN COMMON PART ==============

## 1) DEFINE LOOSE LEPTONS 
process.looseMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string(MUID_LOOSE),
)

process.looseElNoClean = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string(ELID_LOOSE),
)

## 2) LEPTONS CLEANING

process.looseEl = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("looseElNoClean"),
    preselection = cms.string(""),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src = cms.InputTag("boostedMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection = cms.string("userInt('pfMuId')>0"),
           deltaR  = cms.double(0.05),
           checkRecoComponents = cms.bool(False),
           pairCut  = cms.string(""),
           requireNoOverlaps  = cms.bool(True),
        )
    ),
    finalCut = cms.string(''),
)

process.looseLep = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("looseEl", "looseMu")
)

## 3) DEFINE GOOD LEPTONS 

process.goodMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("looseMu"),
    cut = cms.string(MUID_GOOD),
)

process.goodEl = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("looseEl"),
    cut = cms.string(ELID_GOOD),
)

process.goodLep = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("goodMu", "goodEl")
)
process.countGoodLep = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodLep"),
    minNumber = cms.uint32(4),
)


## 4) MAKE Z CANDIDATES

process.goodLL = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goodLep@+ goodLep@-"),
    cut = cms.string(SEL_ANY_Z)
)
process.zll = cms.EDProducer("SkimEvent2LProducer",
    src = cms.InputTag("goodLL"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    doMassRes = cms.bool(False),
)

### =========== BEGIN SIGNAL PART ==============
# Here starts the ZZ step
process.zz = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("zll zll"),
    cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.01 &&"+
                     "deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.01 &&"+
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.01 &&"+
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.01"),
)
process.oneZZ = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zz"),
    minNumber = cms.uint32(1),
)


process.skimEvent4LNoArb = cms.EDProducer("SkimEvent4LProducer",
    src = cms.InputTag("zz"),
    reducedPFCands = cms.InputTag("reducedPFCands"),
    jets = cms.InputTag("slimPatJets"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    isMC = cms.bool(False),
    mcMatch = cms.InputTag(""),
    doswap = cms.bool(True), # sort the two Z's to ensure that Z1 is closest to the nominal mass
    doMELA = cms.bool(True),
    melaQQZZHistos = cms.string("WWAnalysis/AnalysisStep/data/QQZZ8DTemplatesNotNorm.root"),
    doMassRes = cms.bool(True),
)

process.zz4lTreeNoArb = process.zz4lTree.clone(src = cms.InputTag("skimEvent4LNoArb"))

process.selectedZZs1 = cms.EDFilter("SkimEvent4LSelector",
    src = cms.InputTag("skimEvent4LNoArb"),
    cut = cms.string(SEL_ZZ4L_PRE_BEST_1),
    filter = cms.bool(True),
)
process.selectedZZs2 = process.selectedZZs1.clone(
    src = "selectedZZs1",
    cut = SEL_ZZ4L_PRE_BEST_2,
)
process.selectedZZs3 = process.selectedZZs1.clone(
    src = "selectedZZs2",
    cut = SEL_ZZ4L_PRE_BEST_3,
)
process.selectedZZs4 = process.selectedZZs1.clone(
    src = "selectedZZs3",
    cut = SEL_ZZ4L_PRE_BEST_4,
)

process.best4Lpass1 = cms.EDProducer("SkimEvent4LSorterWithTies",
    src = cms.InputTag("selectedZZs4"),
    sortBy = cms.string(SEL_ZZ4L_ARBITRATION_1),
    sortAscending = cms.bool(False),
    maxNumber = cms.uint32(1),
    keepTies = cms.bool(True),
)
process.best4L  = cms.EDProducer("SkimEvent4LSorterWithTies",
    src = cms.InputTag("best4Lpass1"),
    sortBy = cms.string(SEL_ZZ4L_ARBITRATION_2),
    sortAscending = cms.bool(False),
    maxNumber = cms.uint32(1),
    keepTies = cms.bool(False),
)

process.best4LTree = process.zz4lTree.clone(src = cms.InputTag("best4L"))


### =========== CONTROL SELECTION ==========================

process.muonTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("looseMu"),
    sortDescendingBy = cms.string("pt"),
    cut = cms.string(""),
    variables   = cms.PSet(
       pt     = cms.string("pt"),
       eta    = cms.string("eta"),
       phi    = cms.string("phi"),
       pfIsoChHad04      = cms.string("userFloat('pfChHadIso04')"),
       pfIsoNHad04_NoEA  = cms.string("userFloat('muonPFIsoNHad04')"),
       pfIsoPhot04_NoEA  = cms.string("userFloat('muonPFIsoPhoton04')"),
       pfIsoEAtot        = cms.string("userFloat('pfCombIso04EACorrEAtot')"),
       pfIsoComb04EACorr = cms.string("userFloat('pfCombIso04EACorr')"),
       bdtIso = cms.string("userFloat('bdtIso')"),
       sip  = cms.string("userFloat('sip')"),
       rho  = cms.string("userFloat('rhoMu')"),
    ),
    flags = cms.PSet(
       newID  = cms.string("userInt('newID')"),
       prlID = cms.string("userInt('prlID')"),
    ),
    addRunLumiInfo = cms.bool(True),
)

process.electronTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("looseEl"),
    sortDescendingBy = cms.string("pt"),
    cut = cms.string(""),
    variables   = cms.PSet(
       pt     = cms.string("pt"),
       eta    = cms.string("eta"),
       phi    = cms.string("phi"),
       pfIsoChHad04      = cms.string("userFloat('pfChHadIso04')"),
       pfIsoNHad04_NoEA  = cms.string("userFloat('electronPFIsoNHad04')"),
       pfIsoPhot04_NoEA  = cms.string("userFloat('electronPFIsoPhoton04')"),
       pfIsoEApho        = cms.string("userFloat('pfCombIso04EACorrEApho')"),
       pfIsoEAneuHad     = cms.string("userFloat('pfCombIso04EACorrEAneuHad')"),
       pfIsoComb04EACorr = cms.string("userFloat('pfCombIso04EACorr')"),
       bdtID  = cms.string("userFloat('bdtID')"),
       bdtIso = cms.string("userFloat('bdtIso')"),
       sip  = cms.string("userFloat('sip')"),
       rho  = cms.string("userFloat('rhoEl')"),
    ),
    flags = cms.PSet(
       newID = cms.string("userInt('newID')"),
       prlID = cms.string("userInt('prlID')"),
    ),
    addRunLumiInfo = cms.bool(True),
)

process.selectedZ = cms.EDFilter("SkimEvent2LSelector",
    src = cms.InputTag("bestZ1"),
    cut = cms.string(SEL_BEST_Z1),
    filter = cms.bool(True),
)

# Setting up paths

process.common = cms.Sequence(
    process.fourLeptonBlinder +
    process.reboosting +
    process.reskim + 
    process.looseMu +
    process.looseElNoClean + process.looseEl +
    process.looseLep +
    process.goodMu +
    process.goodEl +
    process.goodLep +
    process.goodLL +
    process.zll 
)

process.zzPath = cms.Sequence(
    process.common +
    process.muonTree + process.electronTree +
    process.countGoodLep +
    process.zz     +
    process.oneZZ  +
    process.skimEvent4LNoArb +  process.zz4lTreeNoArb +
    process.selectedZZs1 +
    process.selectedZZs2 +
    process.selectedZZs3 +
    process.selectedZZs4 +
    process.best4Lpass1  +
    process.best4L       +  process.best4LTree 
)


from WWAnalysis.AnalysisStep.zz4l.mcSequences_cff import makeMCSplitPath

### Count events with generator leptons in the pt/eta region
process.ZZ_GenPtEta_Any   = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny +                      process.gen3FilterEta254PtMin5 )
process.ZZ_GenPtEta_4Mu   = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny + process.gen3ZZ4M   + process.gen3FilterEta254PtMin5 )
process.ZZ_GenPtEta_4E    = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny + process.gen3ZZ4E   + process.gen3FilterEta254PtMin5 )
process.ZZ_GenPtEta_2E2Mu = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny + process.gen3ZZ2E2M + process.gen3FilterEta254PtMin5 )
### Normal sequence
process.ZZ_Any   = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny +                      process.zzPath )
process.ZZ_4Mu   = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny + process.gen3ZZ4M   + process.zzPath )
process.ZZ_4E    = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny + process.gen3ZZ4E   + process.zzPath )
process.ZZ_2E2Mu = cms.Path(  process.gen3RecoSeq + process.gen3FilterAny + process.gen3ZZ2E2M + process.zzPath )

process.TFileService = cms.Service("TFileService", fileName = cms.string("hzzTreeMCSync.root"))
