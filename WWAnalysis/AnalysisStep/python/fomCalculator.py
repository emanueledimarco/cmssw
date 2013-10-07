import math
import ROOT
# ROOT.gROOT.SetBatch(True)
from ROOT import RooStats, RooFit
ROOT.RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

def initializeWorkspace(name="w"):
    w = ROOT.RooWorkspace(name,name)
    w.factory("N[0,100]")
    w.factory("B[1]")
    w.factory("dB[0.35]")
    w.factory("S[1]")
    w.factory("dS[0.1]")
    w.factory("r[1,0,50]")
    w.factory("PROD::nuisancePdf(Gaussian::nuisPdfS(thetaS[-5,5],0,1), Gaussian::nuisPdfB(thetaB[-5,5],0,1))")
    w.factory("expr::expB('B*pow(1+dB,thetaB)', B, dB, thetaB)")
    w.factory("expr::expS('S*pow(1+dS,thetaS)', S, dS, thetaS)")
    w.factory("PROD::model_s(Poisson::modelStat_s(N, sum(prod(r,expS),expB), 1), nuisancePdf)")
    w.factory("PROD::model_b(Poisson::modelStat_b(N,                  expB,  1), nuisancePdf)")
    w.factory("Uniform::prior(r)")
    return w

def limitBayesian(w=initializeWorkspace(), B=1, dB=0.35, S=1, dS=0.1):
    w.var("B").setVal(B)
    w.var("dB").setVal(dB)
    w.var("S").setVal(S)
    w.var("dS").setVal(dS)
    w.var("N").setVal(B)

    obs = ROOT.RooArgSet(w.var("N"))
    data = ROOT.RooDataSet("data","data", obs)
    data.add(obs) # insert entry

    rbackup = w.var("r").getMax()

    nuisances = ROOT.RooArgSet()
    if dB != 0:
        nuisances.add(w.var("thetaB"))
        w.var("thetaB").setConstant(False)
    else:
        w.var("thetaB").setConstant(True)

    if dS != 0:
        nuisances.add(w.var("thetaS"))
        w.var("thetaS").setConstant(False)
    else: 
        w.var("thetaS").setConstant(True)
    
    poi = ROOT.RooArgSet(w.var("r"))
    model = w.pdf({True: "model_s", False: "modelStat_s"}[dS>0 or dB>0])
    prior = w.pdf("prior")
    nuisSet = {True:nuisances, False:None}[dS>0 or dB>0]
    bcalc = RooStats.BayesianCalculator(data, model, poi, prior, nuisSet)
    bcalc.SetLeftSideTailFraction(0)
    bcalc.SetConfidenceLevel(0.95)
    interval = bcalc.GetInterval()
    ret = interval.UpperLimit()
    while ret > 0.5*w.var("r").getMax():
        #print "Got ",ret,"too close to the upper bound",w.var("r").getMax()
        if (w.var("r").getMax() > 200): break
        w.var("r").setMax(w.var("r").getMax()*2)
        interval = bcalc.GetInterval()
        ret = interval.UpperLimit()
    w.var("r").setMax(rbackup)
    return ret


