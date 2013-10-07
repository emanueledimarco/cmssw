#!/usr/bin/env python
import os

from WWAnalysis.AnalysisStep.mcAnalysis import *
from ROOT import gROOT, gInterpreter, TCanvas, TFile, TObject

def main():
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] samples.txt cuts.txt plots.txt")
    addMCAnalysisOptions(parser)
    parser.add_option("-O", "--outdir",   dest="outdir",  type="string",       default="plots/",  help="path to directory where plots will be placed (plots/)") 
    parser.add_option("-L", "--nolabel",  dest="nolabel",  action="store_true", default=False,  help="don't add the extra label")
    parser.add_option("-D", "--div",      dest="div",  action="store_true", default=False,  help="also plot the ratio plot")
    parser.add_option("-2", "--two",      dest="two",  action="store_true", default=False,  help="comparison mode for div (instead of Data/MC, need just two samples)")
    parser.add_option(      "--max",      dest="max",  type="float", default=-1.,  help="set the max for the histogram")
    (options, args) = parser.parse_args()

    gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/BasePlot.C+" )
    if options.div: 
        gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/DivPlot.C+" )
        from ROOT import MWLPlot as Plot
    else:
        gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoPlot.C+" )
        from ROOT import LatinoPlot as Plot

#     if options.div: c1 = TCanvas("c1","c1",550,794)
#     else          : c1 = TCanvas("c1","c1",550,600)

    mcAnal = MCAnalysis(args[0],options)
    cf     = CutsFile  (args[1],options)
    pf     = PlotsFile (args[2],options)

    # check if -2 and -D 
    if options.div and options.two:
        raise "Sorry, LatinoPlot has yet to be adapted to handle the division plot"
    # check if -2 and MCAnalysis only has two samples
    if options.two and mcAnal.size() != 2:
        raise "Sorry, in comparison mode, you can only have two samples in your mcAnalysis file"

    for (name,expr,bins,xlabel,ylabel) in pf.plots():
        for i in [len(cf.cuts())-1] if options.final else range(0,len(cf.cuts())):
            reports = mcAnal.getPlotsForCut(name,expr,bins,cf.allCuts(i),False)
            if options.two:
                drawDivisionComparison(reports)
            else:
                drawLatinoPlot(Plot,name,xlabel,cf.cuts()[i][0],reports,options,i)


def drawLatinoPlot(Plot,name,xlabel,autreLabel,reports,options,i):
    rebin = 1
    re1='(.*?)' ; re2='(\\[)' ; re3='(.*?)' ; re4='(\\])' ;
    rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
    for j,whatev in enumerate(reports[reports.keys()[0]]):

        on = "{0}/{1}".format(options.outdir,options.mass)
        if not os.access(on,os.F_OK): os.makedirs(on)

        hs = Plot()
        hs.notPrelim();
        c1 = TCanvas()
        hs.setLumi(options.lumi)
        if options.nolabel: hs.addLabel("")
        else              : hs.addLabel("After "+autreLabel)
        hs.setMass(options.mass)
        hs.setMax(options.max)
        m = rg.search(xlabel)
        if m:
            hs.setLabel(m.group(1))
            hs.setUnits(m.group(3))
        elif xlabel != "": hs.setLabel(xlabel)
        else:              hs.setLabel(name)
        f = TFile(on+"/allHists.root","UPDATE")
        for x in reports:
            if x == 'signal':
                hs.setHWWHist(reports[x][j][1])
                f.cd()
                reports[x][j][1].Write("HWW"+reports[x][j][0],TObject.kOverwrite)
            elif x == 'Data' or x == 'data':
                hs.setDataHist(reports[x][j][1])
                f.cd()
                reports[x][j][1].Write("Data"+reports[x][j][0],TObject.kOverwrite)
            else:
                hs.setMCHist(x,reports[x][j][1])
                f.cd()
                reports[x][j][1].Write(x+reports[x][j][0],TObject.kOverwrite)
        f.Close();

        if options.final: on = "{0}/{1}/{3}_{2}".format(           options.outdir,options.mass,reports[x][j][0],name)
        else:             on = "{0}/{1}/{3}_{2}_{4:0>2.0f}".format(options.outdir,options.mass,reports[x][j][0],name,i)
        c1.SetLogy(False); hs.Draw(c1,rebin,options.div); c1.Print(on+"_lin.png"); c1.Print(on+"_lin.eps");
        c1.SetLogy(True);  hs.Draw(c1,1    ,options.div); c1.Print(on+"_log.png"); c1.Print(on+"_log.eps");

if __name__ == "__main__":
    main()

