import FWCore.ParameterSet.Config as cms

genTauVeto = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("prunedGen"),
    cut = cms.string("status == 3 && 22 <= abs(pdgId) <= 23 && abs(daughter(0).pdgId) == 15"),
    filter = cms.bool(True),
)
gen3Mu = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("prunedGen"),
    cut = cms.string("status == 3 && abs(pdgId) == 13 && abs(motherRef.pdgId) == 23 && abs(motherRef.motherRef.pdgId) == 25"),
)
gen3El = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("prunedGen"),
    cut = cms.string("status == 3 && abs(pdgId) == 11 && abs(motherRef.pdgId) == 23 && abs(motherRef.motherRef.pdgId) == 25"),
)
gen3MM = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("gen3Mu@+ gen3Mu@-"),
    cut = cms.string("mass > 0")
)
gen3EE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("gen3El@+ gen3El@-"),
    cut = cms.string("mass > 0")
)
gen3LL = cms.EDProducer("CandViewMerger", src = cms.VInputTag("gen3MM", "gen3EE"))

gen3LLLL = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("gen3LL gen3LL"),
    cut = cms.string("daughter(0).daughter(0).pt != daughter(1).daughter(0).pt && daughter(0).daughter(1).pt != daughter(1).daughter(1).pt && "+
                     "daughter(0).daughter(1).pt != daughter(1).daughter(0).pt && daughter(0).daughter(0).pt != daughter(1).daughter(1).pt")
)

gen3RecoSeq = cms.Sequence(gen3Mu + gen3El + gen3MM + gen3EE + gen3LL + gen3LLLL)

gen1Mu = gen3Mu.clone(cut = 'status == 1 && abs(pdgId) == 13')
gen1El = gen3El.clone(cut = 'status == 1 && abs(pdgId) == 11')
gen1MM = gen3MM.clone(decay = 'gen1Mu@+ gen1Mu@-')
gen1EE = gen3MM.clone(decay = 'gen1El@+ gen1El@-')
gen1LL = gen3LL.clone(src = ["gen1MM", "gen1EE"])
gen1LLLL = gen3LLLL.clone(decay = "gen1LL gen1LL")
gen1RecoSeq = cms.Sequence(gen1Mu + gen1El + gen1MM + gen1EE + gen1LL + gen1LLLL)


from WWAnalysis.AnalysisStep.zz4l.ptEtaFilters import *
#### Simple filter, select ZZ4L events
gen3FilterAny = cms.EDFilter("CandViewCountFilter", 
    src = cms.InputTag("gen3LLLL"), 
    minNumber = cms.uint32(1)
)
gen1FilterAny = gen3FilterAny.clone(src = "gen1LLLL")

#### Acceptance filter (eta 2.5/2.4, pt 5/3)
gen1FilterEta254PtMin5 = cms.EDFilter("CandViewSelector", 
    src = cms.InputTag("gen1LLLL"), 
    cut = cms.string(adaEtaFilter(2.5,2.4)+" && "+adaPtMinFilter(5,5)),
    filter = cms.bool(True)
)
#### Acceptance filter (eta 2.5/2.4, pt 5/3)
gen3FilterEta254PtMin5 = cms.EDFilter("CandViewSelector", 
    src = cms.InputTag("gen3LLLL"), 
    cut = cms.string(adaEtaFilter(2.5,2.4)+" && "+adaPtMinFilter(5,5)),
    filter = cms.bool(True)
)
#### Geometic-only acceptance
gen1FilterEta254 = gen1FilterEta254PtMin5.clone(cut = adaEtaFilter(2.5,2.4))
gen3FilterEta254 = gen3FilterEta254PtMin5.clone(cut = adaEtaFilter(2.5,2.4))

#### Flavour Filters
gen1ZZ4E = gen1FilterEta254PtMin5.clone(cut = "abs(daughter(0).daughter(0).pdgId()) == 11 && abs(daughter(1).daughter(0).pdgId()) == 11")
gen1ZZ4M = gen1FilterEta254PtMin5.clone(cut = "abs(daughter(0).daughter(0).pdgId()) == 13 && abs(daughter(1).daughter(0).pdgId()) == 13")
gen1ZZ2E2M = gen1FilterEta254PtMin5.clone(cut = "abs(daughter(0).daughter(0).pdgId()) != abs(daughter(1).daughter(0).pdgId())")

#### Flavour Filters
gen3ZZ4E = gen3FilterEta254PtMin5.clone(cut = "abs(daughter(0).daughter(0).pdgId()) == 11 && abs(daughter(1).daughter(0).pdgId()) == 11")
gen3ZZ4M = gen3FilterEta254PtMin5.clone(cut = "abs(daughter(0).daughter(0).pdgId()) == 13 && abs(daughter(1).daughter(0).pdgId()) == 13")
gen3ZZ2E2M = gen3FilterEta254PtMin5.clone(cut = "abs(daughter(0).daughter(0).pdgId()) != abs(daughter(1).daughter(0).pdgId())")

def makeMCSplitPath(process, name, seq, status=1, tauVeto=True, genFilters="none"):  
    " genFilters can be 'none', 'eta', 'ptEta' "
    seqs = { 1:process.gen1RecoSeq, 3:process.gen3RecoSeq }
    ptFilt  = { 1:process.gen1FilterEta254PtMin5, 3:process.gen3FilterEta254PtMin5 }
    etaFilt = { 1:process.gen1FilterEta254, 3:process.gen3FilterEta254 }
    anyFilt = { 1:process.gen1FilterAny, 3:process.gen3FilterAny }
    filters = { 1:{'4E':gen1ZZ4E, '4M':gen1ZZ4M, '2E2M':gen1ZZ2E2M },
                3:{'4E':gen3ZZ4E, '4M':gen3ZZ4M, '2E2M':gen3ZZ2E2M } }
    seq0 = seqs[status]
    if tauVeto: seq0 = ~process.genTauVeto + seq0
    if genFilters == "ptEta": seq0 += ptFilt[status]
    elif genFilters == "eta": seq0 += etaFilt[status]
    elif genFilters == "any": seq0 += anyFilt[status]
    elif genFilters != "none": raise RuntimeError, "Bogus genFilter"
    ## Part below is a bit convoluted, but the cmssw sequences play tricks if they're not in the process
    setattr(process, name+"_All_GSeq",  cms.Sequence(seq0))
    setattr(process, name+"_4E_GSeq",   cms.Sequence(getattr(process,name+"_All_GSeq") + filters[status]['4E'  ]))
    setattr(process, name+"_2E2M_GSeq", cms.Sequence(getattr(process,name+"_All_GSeq") + filters[status]['2E2M']))
    setattr(process, name+"_4M_GSeq",   cms.Sequence(getattr(process,name+"_All_GSeq") + filters[status]['4M'  ]))
    for X in "_All _4E _2E2M _4M".split():
        setattr(process, name+X, cms.Path(getattr(process,name+X+"_GSeq") + seq))
