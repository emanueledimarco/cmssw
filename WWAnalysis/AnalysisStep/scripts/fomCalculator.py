#!/usr/bin/env python
from ROOT import gROOT, TCanvas, TH1F, TFile, TDirectoryFile, TTree, AddressOf
from ROOT import kCyan, kGreen, kMagenta, kBlue
from ROOT import RooWorkspace, RooDataSet, RooArgSet, RooMsgService, RooStats, RooFit
import math
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

from WWAnalysis.AnalysisStep.fomCalculator import *

import sys
import math 

def main():
    w = initializeWorkspace()
#     wp95 = [0.16,0.85+1.07+0.02]
#     wp90 = [0.14,0.77]
#     wp85 = [0.12,0.70]
#     lh95 = [0.16,0.86+0.53]
#     lh90 = [0.15,0.82]
#     print limitBayesian(w,wp95[1],0.35,wp95[0],0.1) 
#     print limitBayesian(w,wp90[1],0.35,wp90[0],0.1) 
#     print limitBayesian(w,wp85[1],0.35,wp85[0],0.1) 
#     print limitBayesian(w,lh95[1],0.35,lh95[0],0.1) 
#     print limitBayesian(w,lh90[1],0.35,lh90[0],0.1) 
    print limitBayesian(w,  87.71 ,0.35,  27.62, 0.1 )
    print limitBayesian(w,  108.45,0.35,  30.08, 0.1 )


def runTest():
    w = initializeWorkspace()
    # example optimization for a cut on 'x'
    #   eff(S) = 1-x
    #   eff(B) = exp(-5*x)
    # S(X=0) = 4    
    # B(x=0) = 30  
    # dS/S = 15%, dB/B = 35%
    # let's scan X, print values of limit and save them to a tree
    S0 = 4
    B0 = 20
    dS = 0.1
    dB = 0.35
    fOut = TFile.Open("figureOfMerit.root", "RECREATE")
    tOut = TTree("fom","fom")
#     x = 0 
#     S = 0
#     B = 0
#     limit = 0
#     naive = 0
#     tOut.Branch("x", AddressOf(x), "X/D")  
#     tOut.Branch("S", AddressOf(S), "S/D") 
#     tOut.Branch("B", AddressOf(B), "B/D") 
#     tOut.Branch("limit", AddressOf(limit), "limit/D") 
#     tOut.Branch("naive", AddressOf(naive), "naive/D")
    gROOT.ProcessLine("struct vars { float x; float S; float B; float limit; float naive; };" );
    from ROOT import vars
    var = vars();
    tOut.Branch("var",AddressOf(var),"x/F:S:B:limit:naive")

    print " %6s  %6s  %6s  %6s  %6s" % (" X ", " S ", " B ", "limit", "naive")
    # Tree's don't freakin work !
    for x in range(0,4):
        var.x = float(x)/20.
        var.S = S0*(1-var.x)
        var.B = B0*math.exp(-5*var.x)
        var.limit = limitBayesian(w,var.B,dB,var.S,dS)
        var.naive = 2*math.sqrt(var.B+(var.B*dB*var.B*dB))/var.S
        print " %6.2f  %6.2f  %6.2f  %6.2f  %6.2f" % (var.x, var.S, var.B, var.limit, var.naive)
        tOut.Fill()
    tOut.Write() 
    fOut.Close()




if __name__ == "__main__":
    w=initializeWorkspace()
    if len(sys.argv) != 5: raise RuntimeError, "Usage: fomCalculator.py S dS B dB "
    (S,dS,B,dB,) = [float(x) for x in sys.argv[1:]]
    print limitBayesian(w,B,dB,S,dS)
