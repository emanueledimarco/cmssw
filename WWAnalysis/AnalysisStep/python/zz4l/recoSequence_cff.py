import FWCore.ParameterSet.Config as cms


cleanedEl = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("boostedElectrons"),
    preselection = cms.string(""),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src = cms.InputTag("boostedMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection = cms.string("pt > 2 && (isGlobalMuon || numberOfMatches > 1)"),
           deltaR  = cms.double(0.01),
           checkRecoComponents = cms.bool(False),
           pairCut  = cms.string(""),
           requireNoOverlaps  = cms.bool(True),
        )
    ),
    finalCut = cms.string(''),
)


recMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string(""),
)

recEl = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("cleanedEl"),
    cut = cms.string(""),
)

recMM = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("recMu@+ recMu@-"),
    cut = cms.string("mass > 0")
)
recEE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("recEl@+ recEl@-"),
    #decay = cms.string("cleanedEl@+ cleanedEl@-"),
    cut = cms.string("mass > 0")
)
recLL = cms.EDProducer("CandViewMerger", 
    src = cms.VInputTag("recMM", "recEE")
)

recE1Filter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("cleanedEl"),
    minNumber = cms.uint32(2)
)


recE2Filter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("recEl"),
    minNumber = cms.uint32(2)
)

recLLLLdups = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("recLL recLL"),
    cut = cms.string("daughter(0).daughter(0).pt != daughter(1).daughter(0).pt && daughter(0).daughter(1).pt != daughter(1).daughter(1).pt && "+
                     "daughter(0).daughter(1).pt != daughter(1).daughter(0).pt && daughter(0).daughter(0).pt != daughter(1).daughter(1).pt")
)

recLLLL = cms.EDProducer("CompositeCandidateDuplicateCleaner",
    src = cms.InputTag("recLLLLdups"),
    algo = cms.string("BestZ1Mass"),
)

recoSeq = cms.Sequence(
    cleanedEl + 
    recMu + recEl + 
    recMM + recEE + 
    #recE1Filter + recE2Filter + 
    recLL + 
    recLLLLdups +
    recLLLL
)

from WWAnalysis.AnalysisStep.zz4l.ptEtaFilters import *
#### Simple filter, select ZZ4L events
recFilterAny = cms.EDFilter("CandViewCountFilter", 
    src = cms.InputTag("recLLLL"), 
    minNumber = cms.uint32(1)
)

#### Acceptance filter (eta 2.5/2.4, pt 5/3)
recFilterEta254PtMin53 = cms.EDFilter("CandViewSelector", 
    src = cms.InputTag("recLLLL"), 
    cut = cms.string(adaEtaFilter(2.5,2.4)+" && "+adaPtMinFilter(5,3)),
    filter = cms.bool(True)
)

#### Flavour Filters
recZZ4E = recFilterEta254PtMin53.clone(cut = "abs(daughter(0).daughter(0).pdgId()) == 11 && abs(daughter(1).daughter(0).pdgId()) == 11")
recZZ4M = recFilterEta254PtMin53.clone(cut = "abs(daughter(0).daughter(0).pdgId()) == 13 && abs(daughter(1).daughter(0).pdgId()) == 13")
recZZ2E2M = recFilterEta254PtMin53.clone(cut = "abs(daughter(0).daughter(0).pdgId()) != abs(daughter(1).daughter(0).pdgId())")
