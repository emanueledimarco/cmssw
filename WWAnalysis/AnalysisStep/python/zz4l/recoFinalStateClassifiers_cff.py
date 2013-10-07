import FWCore.ParameterSet.Config as cms

z1EE = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("bestZ"),
    cut = cms.string("abs(daughter(0).pdgId) == 11"),
    filter = cms.bool(True),
)
z1MM = z1EE.clone(
    cut = "abs(daughter(0).pdgId) == 13"
)

# for 3-path logic
oneZEE = cms.EDFilter("SkimEvent2LSelector",
    src = cms.InputTag("zll"),
    cut = cms.string("abs(daughter(0).pdgId) == 11"),
    filter = cms.bool(True),
)
oneZMM = oneZEE.clone(cut = "abs(daughter(0).pdgId) == 13")
bestZMM = cms.EDFilter("SkimEvent2LSorter",
    src = cms.InputTag("oneZMM"),
    sortBy = cms.string("abs(mass - 91.188)"),
    sortAscending = cms.bool(True),
    maxNumber = cms.uint32(1),
)
bestZEE = bestZMM.clone(src = "oneZEE")
selectedZEE = oneZEE.clone(src = "bestZEE", cut = "40 < mass < 120")
selectedZMM = oneZMM.clone(src = "bestZMM", cut = "40 < mass < 120")


countGoodLep4M = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodMu"),
    minNumber = cms.uint32(4),
)
countGoodLep4E = countGoodLep4M.clone(
    src = cms.InputTag("goodEl"),
)
countGoodLep1E = countGoodLep4E.clone(minNumber = 1)
countGoodLep1M = countGoodLep4M.clone(minNumber = 1)
countGoodLep2E = countGoodLep4E.clone(minNumber = 2)
countGoodLep2M = countGoodLep4M.clone(minNumber = 2)
countGoodLep3E = countGoodLep4E.clone(minNumber = 3)
countGoodLep3M = countGoodLep4M.clone(minNumber = 3)
countGoodLep2E2M = cms.Sequence(countGoodLep2E + countGoodLep2M)
countGoodLep3EM  = cms.Sequence(countGoodLep1E + countGoodLep1M)
countLooseLep1M = countGoodLep1M.clone(src = 'looseMu')
countLooseLep1E = countGoodLep1E.clone(src = 'looseEl')
countLooseLep2M = countGoodLep2M.clone(src = 'looseMu')
countLooseLep2E = countGoodLep2E.clone(src = 'looseEl')
countLooseLep3M = countGoodLep3M.clone(src = 'looseMu')
countLooseLep3E = countGoodLep3E.clone(src = 'looseEl')
countLooseLep4M = countGoodLep4M.clone(src = 'looseMu')
countLooseLep4E = countGoodLep4E.clone(src = 'looseEl')
countLooseLep2E2M = cms.Sequence(countLooseLep2E + countLooseLep2M)
countLooseLep3EM  = cms.Sequence(countLooseLep1E + countLooseLep1M)

zz4E = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("best4L"),
    cut = cms.string("abs(daughter(0).daughter(0).pdgId) == 11 && abs(daughter(1).daughter(0).pdgId) == 11"),
    filter = cms.bool(True),
)
zz4M = zz4E.clone(
    cut = cms.string("abs(daughter(0).daughter(0).pdgId) == 13 && abs(daughter(1).daughter(0).pdgId) == 13"),
)
zz2E2M = zz4E.clone(
    cut = cms.string("abs(daughter(0).daughter(0).pdgId) != abs(daughter(1).daughter(0).pdgId)"),
)

def makeSplittedPaths4L(process, name, triggerFilterName=None, doThreePathLogic=True):
    from PhysicsTools.PatAlgos.tools.helpers import listModules, cloneProcessingSnippet
    path = getattr(process, name)
    modules = listModules(path)
    seq4E   = cms.Sequence(sum(modules[1:], modules[0]))
    seq4M   = cms.Sequence(sum(modules[1:], modules[0]))
    seq2E2M = cms.Sequence(sum(modules[1:], modules[0]))
    if process.bestZ in modules:
        seq4E.replace(process.bestZ, process.bestZ + process.z1EE)
        seq4M.replace(process.bestZ, process.bestZ + process.z1MM)
    if process.countGoodLep3 in modules:
        seq4E.replace(process.countGoodLep3, process.countGoodLep3E)
        seq4M.replace(process.countGoodLep3, process.countGoodLep3M)
        seq2E2M.replace(process.countGoodLep3, process.countGoodLep3 + countGoodLep3EM)
    if process.countGoodLep in modules:
        seq4E.replace(process.countGoodLep, process.countGoodLep4E)
        seq4M.replace(process.countGoodLep, process.countGoodLep4M)
        seq2E2M.replace(process.countGoodLep, process.countGoodLep2E2M)
    if process.countLooseLep3 in modules:
        seq4E.replace(process.countLooseLep3, process.countLooseLep3E)
        seq4M.replace(process.countLooseLep3, process.countLooseLep3M)
        seq2E2M.replace(process.countLooseLep3, process.countLooseLep3 + countLooseLep3EM)
    if process.countLooseLep in modules:
        seq4E.replace(process.countLooseLep, process.countLooseLep4E)
        seq4M.replace(process.countLooseLep, process.countLooseLep4M)
        seq2E2M.replace(process.countLooseLep, process.countLooseLep2E2M)
    if process.best4L in modules:
        seq4E.replace(process.best4L, process.best4L + process.zz4E)
        seq4M.replace(process.best4L, process.best4L + process.zz4M)
        seq2E2M.replace(process.best4L, process.best4L + process.zz2E2M)
    setattr(process, name+"_4E", cms.Path(seq4E))
    setattr(process, name+"_4M", cms.Path(seq4M))
    setattr(process, name+"_2E2M", cms.Path(seq2E2M))
    ## ---------- the three-path logic comes here ----------------------
    if doThreePathLogic:
        modules4E = listModules(getattr(process, name+"_4E"))
        modules4M = listModules(getattr(process, name+"_4M"))
        seq4E_3P = cms.Sequence(sum(modules4E[1:], modules4E[0]))
        seq4M_3P = cms.Sequence(sum(modules4M[1:], modules4M[0]))
        ## ---- Change trigger filter ----
        if triggerFilterName:
            seq4E_3P.replace(getattr(process,triggerFilterName), getattr(process,triggerFilterName+"_El"))
            seq4M_3P.replace(getattr(process,triggerFilterName), getattr(process,triggerFilterName+"_Mu"))
        ## ---- Change Z1 definition -----
        if process.oneZ in modules:
            seq4E_3P.replace(process.oneZ, process.oneZEE + process.bestZEE + process.selectedZEE)
            seq4M_3P.replace(process.oneZ, process.oneZMM + process.bestZMM + process.selectedZMM)
        for m in [process.bestZ, process.selectedZ1]: 
            if m in modules4E: seq4E_3P.remove(m)
            if m in modules4M: seq4M_3P.remove(m)
        ## ---- Change combinatoric -----
        process.zz3pathLogic4M = process.zz.clone(decay = "selectedZMM oneZMM" if process.zz.decay.value() == "selectedZ1 zll" else "oneZMM oneZMM")
        process.zz3pathLogic4E = process.zz.clone(decay = "selectedZEE oneZEE" if process.zz.decay.value() == "selectedZ1 zll" else "oneZEE oneZEE")
        process.oneZZ3pathLogic4E = process.oneZZ.clone(src = "zz3pathLogic4E")
        process.oneZZ3pathLogic4M = process.oneZZ.clone(src = "zz3pathLogic4M")
        seq4E_3P.replace(process.zz, process.zz3pathLogic4E)
        seq4M_3P.replace(process.zz, process.zz3pathLogic4M)
        seq4E_3P.replace(process.oneZZ, process.oneZZ3pathLogic4E)
        seq4M_3P.replace(process.oneZZ, process.oneZZ3pathLogic4M)
        process.skimEvent4LNoArb3pathLogic4E = process.skimEvent4LNoArb.clone(src = "zz3pathLogic4E")
        process.skimEvent4LNoArb3pathLogic4M = process.skimEvent4LNoArb.clone(src = "zz3pathLogic4M")
        process.zzCombinatoric_4E = cloneProcessingSnippet(process,process.zzCombinatoric,"3pathLogic4E")
        process.zzCombinatoric_4M = cloneProcessingSnippet(process,process.zzCombinatoric,"3pathLogic4M")
        modulesZZC = [process.skimEvent4LNoArb] + listModules(process.zzCombinatoric)
        for m in modulesZZC:
            if m in modules4E: seq4E_3P.replace(m, getattr(process,m.label_() + "3pathLogic4E"))
            if m in modules4M: seq4M_3P.replace(m, getattr(process,m.label_() + "3pathLogic4M"))
        if process.best4Lpass1.src.value() == "skimEvent4LNoArb":
            process.best4Lpass13pathLogic4M.src = "skimEvent4LNoArb3pathLogic4M"
            process.best4Lpass13pathLogic4E.src = "skimEvent4LNoArb3pathLogic4E"
        if process.selectedZZs1.src.value() == "skimEvent4LNoArb":
            process.selectedZZs13pathLogic4M.src = "skimEvent4LNoArb3pathLogic4M"
            process.selectedZZs13pathLogic4E.src = "skimEvent4LNoArb3pathLogic4E"
        for m in [process.z1EE, process.z1MM, process.zz4E, process.zz4M, process.zz4lTree, process.zz4lTreeNoArb]: 
            if m in modules4E: seq4E_3P.remove(m)
            if m in modules4M: seq4M_3P.remove(m)
        setattr(process, name+"_4E_3Path", cms.Path(seq4E_3P))
        setattr(process, name+"_4M_3Path", cms.Path(seq4M_3P))
