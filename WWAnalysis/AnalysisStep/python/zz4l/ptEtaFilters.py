def etaFilter(x):
    return " && ".join(["abs(daughter(%d).daughter(%d).eta) < %g" % (i,j,x) for (i,j) in [ (0,0), (0,1), (1,0), (1,1) ]])
def ptMinFilter(x):
    return "min(min(daughter(0).daughter(0).pt,daughter(0).daughter(1).pt), min(daughter(1).daughter(0).pt,daughter(1).daughter(1).pt)) > %g" % x
def adaPtMinFilter(ptEle,ptMu):
    (ptMin,ptMax,id) = (ptMu,ptEle,13) if ptMu <= ptEle else (ptEle,ptMu,11)
    DAU="daughter(%d).daughter(%d)"; ijseq = [ (0,0), (0,1), (1,0), (1,1) ]
    return " && ".join([("( "+DAU+".pt > %g || abs("+DAU+".pdgId) == %d && "+DAU+".pt > %g)") % (i,j,ptMax,i,j,id,i,j,ptMin) for (i,j) in ijseq ])
def adaEtaFilter(etaEle,etaMu):
    (etaMin,etaMax,id) = (etaMu,etaEle,11) if etaMu <= etaEle else (etaEle,ptMu,13)
    DAU="daughter(%d).daughter(%d)"; ijseq = [ (0,0), (0,1), (1,0), (1,1) ]
    return " && ".join([("( abs("+DAU+".eta) < %g || abs("+DAU+".pdgId) == %d && abs("+DAU+".eta) < %g)") % (i,j,etaMin,i,j,id,i,j,etaMax) for (i,j) in ijseq ])
