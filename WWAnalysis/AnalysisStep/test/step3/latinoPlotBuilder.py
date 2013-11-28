#!/usr/bin/env python
import os

from WWAnalysis.AnalysisStep.mcAnalysis import *


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] samples.txt cuts.txt plots.txt")
    addMCAnalysisOptions(parser)
    parser.add_option("-O", "--outdir",   dest="outdir",  type="string",       default="plots/",  help="path to directory where plots will be placed (plots/)") 
    (options, args) = parser.parse_args()

    from ROOT import gROOT, gInterpreter, TCanvas
    gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoPlot.C+" )
    gInterpreter.ExecuteMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoStyle.C" )
    from ROOT import LatinoPlot
    c1 = TCanvas()

    re1='(.*?)' ; re2='(\\[)' ; re3='(.*?)' ; re4='(\\])' ;
    rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
    rebin = 1

    mcAnal = MCAnalysis(args[0],options)
    cf     = CutsFile  (args[1],options)
    pf     = PlotsFile (args[2],options)
    for (name,expr,bins,xlabel,ylabel) in pf.plots():
        for i in range(0,len(cf.cuts())):
            reports = mcAnal.getPlotsForCut(name,expr,bins,cf.allCuts(),False)
            for j,whatev in enumerate(reports[reports.keys()[0]]):
                hs = LatinoPlot()
                hs.setLumi(options.lumi)
                hs.addLabel("After "+cf.cuts()[i][0])
                hs.setMass(options.mass)
                m = rg.search(xlabel)
                if m:
                    hs.setLabel(m.group(1))
                    hs.setUnits(m.group(3))
                elif xlabel != "": hs.setLabel(xlabel)
                else:              hs.setLabel(name)
                for x in reports:
                    if x == 'signal':
                        hs.setHWWHist(reports[x][j][1])
                    elif getattr(hs,'set{0}Hist'.format(x),None) != None:
                        getattr(hs,'set{0}Hist'.format(x))(reports[x][j][1])
                on = "{0}/{1}/{2}/{3}".format(options.outdir,options.mass,reports[x][j][0],name)
                if not os.access(on,os.F_OK): os.makedirs(on)
                on = "{0}/{1}/{2}/{3}/{4:0>2.0f}".format(options.outdir,options.mass,reports[x][j][0],name,i)
                print on
                c1.SetLogy(False); hs.Draw(rebin); c1.Print(on+".lin.png"); c1.Print(on+".lin.eps");
                c1.SetLogy(True);  hs.Draw(1);     c1.Print(on+".log.png"); c1.Print(on+".log.eps");


