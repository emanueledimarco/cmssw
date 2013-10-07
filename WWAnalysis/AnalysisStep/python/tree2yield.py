#!/usr/bin/env python
from math import *
import re
import os.path
import ROOT
from copy import *
ROOT.gROOT.SetBatch(True)

MM=0;EE=1;EM=2;ME=3;SF=4;OF=5;ALL=6;
class CutsFile:
    def __init__(self,txtfileOrCuts,options=None):
        if type(txtfileOrCuts) == list:
            self._cuts = deepcopy(txtfileOrCuts[:])
        elif isinstance(txtfileOrCuts,CutsFile):
            self._cuts = deepcopy(txtfileOrCuts.cuts())
        else:
            self._cuts = []
            file = open(txtfileOrCuts, "r")
            if not file: raise RuntimeError, "Cannot open "+txtfileOrCuts+"\n"
            for cr,cn,cv in options.cutsToAdd:
                if re.match(cr,"entry point"): self._cuts.append((cn,cv))
            for line in file:
                if len(line.strip()) == 0 or line.strip()[0] == '#': continue
                (name,cut) = [x.strip() for x in line.split(":")]
                if name == "entry point" and cut == "1": continue
                if options.startCut and not re.search(options.startCut,name): continue
                if options.startCut and re.search(options.startCut,name): options.startCut = None
                self._cuts.append((name,cut))
                for cr,cn,cv in options.cutsToAdd:
                    if re.match(cr,name): self._cuts.append((cn,cv))
                if options.upToCut and re.search(options.upToCut,name):
                    break
            for ci in options.cutsToInvert:  self.invert(ci)
            for ci in options.cutsToExclude: self.remove(ci)
            for cr,cn,cv in options.cutsToReplace: self.replace(cr,cn,cv)
    def __str__(self):
        newstring = ""
        for cut in self._cuts:
            newstring += "{0} : {1}\n".format(cut[0],cut[1])
        return newstring[:-1]
    def remove(self,cut):
        self._cuts = [(cn,cv) for (cn,cv) in self._cuts if not re.search(cut,cn)]
        return self
    def invert(self,cut):
        for i,(cn,cv) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                if cn.startswith("not ") and re.match(r"!\(.*\)", cv):
                    self._cuts[i] = (cn[4:], cv[2:-1])
                else:
                    self._cuts[i] = ("not "+cn, "!("+cv+")")
        return self
    def replace(self,cut,newname,newcut):       
        for i,(cn,cv) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                self._cuts[i] = (newname, newcut)
        return self
    def cuts(self):
        return self._cuts[:]
    def nMinusOne(self):
        return CutsFile(self.nMinusOneCuts())
    def nMinusOneCuts(self):
        ret = []
        for cn,cv in self._cuts:
            nm1 = " && ".join("(%s)" % cv1 for cn1,cv1 in self._cuts if cn1 != cn)
            ret.append(("all but "+cn, nm1))
        return ret
    def allCuts(self,n=-1):
        return " && ".join("(%s)" % x[1] for x in (self._cuts[0:n+1] if n != -1 and n+1 < len(self._cuts) else self._cuts))
    def addAfter(self,cut,newname,newcut):
        for i,(cn,cv) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                self._cuts.insert(i+1,(newname, newcut))
                break
        return self
    def insert(self,index,newname,newcut):
        self._cuts.insert(index,(newname, newcut))
        return self
    def add(self,newname,newcut):
        self._cuts.append((newname,newcut))
        return self

class PlotsFile:
    def __init__(self,txtfileOrPlots,options=None):
        if type(txtfileOrPlots) == list:
            self._plots = txtfileOrPlots[:]
        else:
            self._plots = []
            file = open(txtfileOrPlots, "r")
            if not file: raise RuntimeError, "Cannot open "+txtfileOrPlots+"\n"
            for line in file:
                if re.match("\\s*#.*", line): continue
                if ("%d" in line): line = line % options.mass
                words = [x.strip() for x in line.split(":")]
                (name,expr,bins) = [ x for x in words[0:3] ]
                xlabel = words[3] if len(words) > 3 else ""
                ylabel = words[4] if len(words) > 4 else "Entries"
                self._plots.append((name,expr,bins,xlabel,ylabel))
    def plots(self):
        return self._plots[:]
    def __str__(self):
        newstring = ""
        for plot in self._plots:
            newstring += "{0} : {1} : {2}\n".format(plot[0],plot[1],plot[2])
        return newstring[:-1]

class TreeToYield:
    def __init__(self,root,options,scaleFactor=1.0):
        self._fname = root
        self._name  = os.path.basename(root).replace("tree_","").replace(".root","").replace("latino_","")
        self._tfile = ROOT.TFile.Open(root)
        if not self._tfile: raise RuntimeError, "Cannot open %s\n" % root
        self._options = options
        t = self._tfile.Get(options.tree)
        if not t: raise RuntimeError, "Cannot find tree %s in file %s\n" % (options.tree, root)
        self._tree  = t
        self._chancuts = []
        for i,h in enumerate(("mumu","elel","elmu","muel")):
            self._chancuts.append((h,"channel==%i"%i))
        self._weight  = (options.weight and '2011A' not in self._name and '2011B' not in self._name )
        self._weightString  = options.weightString
        self._scaleFactor = scaleFactor
        self._treeMVA = None
        self._treeEff = None
        if options.mva: self.attachMVA(options.mva)
        if options.eff: self.attachEff(options.eff)
        if options.keysHist:  ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    def setScaleFactor(self,scaleFactor):
        self._scaleFactor = scaleFactor
    def getScaleFactor(self):
        return self._scaleFactor
    def attachMVA(self,name):
        if self._treeMVA: # must first detach
            self._tree.RemoveFriend(self._treeMVA)
            self._treeMVA = None
        self._fnameMVA = self._fname.replace(".root","_"+name+".root")
        self._tfileMVA = ROOT.TFile.Open(self._fnameMVA)
        if not self._tfileMVA: raise RuntimeError, "Cannot open %s\n" % self._fnameMVA
        t = self._tfileMVA.Get('latinobdt')
        if not t: raise RuntimeError, "Cannot find tree %s in file %s\n" % (name, self._fnameMVA)
        self._treeMVA = t
        self._tree.AddFriend(self._treeMVA)
    def attachEff(self,name):
        if self._treeEff: # must first detach
            self._tree.RemoveFriend(self._treeEff)
            self._treeEff = None
        self._fnameEff = self._fname.replace(".root","."+name+".root")
        self._tfileEff = ROOT.TFile.Open(self._fnameEff)
        if not self._tfileEff: raise RuntimeError, "Cannot open %s\n" % self._fnameEff
        t = self._tfileEff.Get(name)
        if not t: raise RuntimeError, "Cannot find tree %s in file %s\n" % (name, self._fnameEff)
        self._treeEff = t
        self._tree.AddFriend(self._treeEff)
    def name(self):
        return self._name
    def getYields(self,cuts):
        report = []; cut = ""
        cutseq = [ ['entry point','1'] ]
        sequential = False
        if self._options.nMinusOne: 
            cutseq = cuts.nMinusOneCuts()
            cutseq += [ ['all',cuts.allCuts()] ]
            sequential = False
        elif self._options.final:
            cutseq += [ ['all', cuts.allCuts()] ]
        else:
            cutseq += cuts.cuts();
            sequential = True
        for cn,cv in cutseq:
            if sequential:
                if cut: cut += " && "
                cut += "(%s)" % cv
            else:
                cut = cv
            report.append((cn,self._getYields(cut)))
        return report
    def prettyPrint(self,report):

        hyposToExclude = [] if self._options.flavors else ['sf','of']
        hyposToPrint = [ hypo for (hypo,a,b) in report[0][1] if hypo not in hyposToExclude ] 
        nch  = len(hyposToPrint)

        # maximum length of the cut descriptions
        clen = max([len(cut) for cut,yields in report]) + 3
        cfmt = "%%-%ds" % clen;

        fmtlen = 12
        nfmtL = "    %8d"
        nfmtS = "    %8.3f" if self._weight else nfmtL

        if self._options.errors:
            nfmtS+=u"%8.3f"
            nfmtL+=u"%8.3f"
            fmtlen+=8
        if self._options.fractions:
            nfmtS+="%7.1f%%"
            nfmtL+="%7.1f%%"
            fmtlen+=8

        print "cut".center(clen),
        for hypo in hyposToPrint:
            print ("    "+hypo).center(fmtlen),
        print ""
        print "-"*((fmtlen+1)*nch+clen)
        for i,(cut,yields) in enumerate(report):
            print cfmt % cut,
            for j,(hypo,nev,err) in enumerate(yields):
                if hypo not in hyposToPrint: continue
                den = report[i-1][1][j][1] if i>0 else 0
                fraction = nev/float(den) if den > 0 else 1
                if self._options.nMinusOne: 
                    fraction = report[-1][1][j][1]/nev if nev > 0 else 1
                toPrint = (nev,)
                if self._options.errors:    toPrint+=(err,)
                if self._options.fractions: toPrint+=(fraction*100,)
                if self._weight and nev < 1000: print nfmtS % toPrint,
                else                          : print nfmtL % toPrint,
            print ""
    def getPlots(self,plots,cut):
        ret = [ [name, self.getPlots(name,expr,bins,cut)] for (name,expr,bins) in plots.plots()]
        return ret
    def getPlots(self,name,expr,bins,cut):
        plots = [ [k,self._getPlot(self._tree,name+"_"+k,expr,bins,("("+cut+")&&(%s)")%c)] for (k,c) in self._chancuts ]
        hall  = plots[0][1].Clone(name+"_all"); hall.Reset()
        for k,h in plots: hall.Add(h)
        all   = [ ['all', hall] ]
        if self._options.inclusive:
            plots = all
        else:
            hsf  = plots[0][1].Clone(name+"_sf"); hsf.Reset()
            hof  = plots[0][1].Clone(name+"_of"); hof.Reset()
            for k,h in [ (x,y) for x,y in plots if x == "elel" or x == "mumu"]: hsf.Add(h)
            for k,h in [ (x,y) for x,y in plots if x == "elmu" or x == "muel"]: hof.Add(h)
            plots += [ ['sf', hsf] ]
            plots += [ ['of', hof] ]
            plots += all
        return plots
    def dumpEvents(self,cut,vars=['run','lumi','event','dataset','channel']):
        for (k,c) in self._chancuts:
            print "Dump for channel ",k
            self._tree.SetScanField(100000)
            self._tree.Scan(":".join(vars), ("("+cut+")&&(%s)")%c,"colsize=20")
            print
    def getAverageWeight(self,cut):
        if not self._weight: return 1.0
        nev = 0; sumw = 0;
        for (k,c) in self._chancuts: 
            (n,sw) = self._getNumAndWeight(self._tree,("("+cut+")&&(%s)")%c)
            nev += n; sumw += sw
        if nev == 0: return 0
        return sumw/nev;
    def _getYields(self,cut):
        yields = [ [k] + self._getYield(self._tree,("("+cut+")&&(%s)")%c) for (k,c) in self._chancuts ]
        all    = [ ['all', sum(x for h,x,e in yields), sqrt(sum(e*e for h,x,e in yields))] ]
        if self._options.inclusive:
            yields = all
        else:
            sf     = [ ['sf', sum(x for h,x,e in yields if h == 'elel' or h == 'mumu'), sqrt(sum(e*e for h,x,e in yields if h == 'elel' or h == 'mumu'))] ]
            of     = [ ['of', sum(x for h,x,e in yields if h == 'muel' or h == 'elmu'), sqrt(sum(e*e for h,x,e in yields if h == 'muel' or h == 'elmu'))] ]
            yields += sf
            yields += of
            yields += all
        return yields
    def _getYield(self,tree,cut):
        if self._weight:
            histo = self._getPlot(tree,"dummy","0.5","1,0.,1.",cut)
            return [ histo.GetBinContent(1), histo.GetBinError(1) ]
        else: 
            npass = tree.Draw("1",cut,"goff");
            return [ npass, sqrt(npass) ]
    def _getNumAndWeight(self,tree,cut):
            histo = ROOT.TH1F("dummy","dummy",1,0.,1.)
            histo.Sumw2()
            nev = tree.Draw("0.5>>dummy", "(%s)*(%s)*(%s)" % (self._weightString,self._scaleFactor,cut) ,"goff")
            if nev == 0: return (0,0)
            sumw = histo.GetBinContent(1)*self._options.lumi
            histo.Delete()
            return (nev,sumw)
    def _getPlot(self,tree,name,expr,bins,cut):
            if self._weight: cut = "(%s)*(%s)*(%s)*(%s)" % (self._weightString,self._options.lumi, self._scaleFactor, cut)
            (nb,xmin,xmax) = bins.split(",")
            if ROOT.gROOT.FindObject("dummy") != None: ROOT.gROOT.FindObject("dummy").Delete()
#             histo = ROOT.TH1F("dummy","dummy",int(nb),float(xmin),float(xmax))
#             histo.Sumw2()
#             nev = tree.Draw("%s>>%s" % (expr,"dummy"), cut ,"goff")
#Implementation to use when TH1Keys is not returning the right normalization
#           if self._options.keysHist and self._scaleFactor!=1:
#               newScale = histo.Integral()
#               if ROOT.gROOT.FindObject("dummy2") != None: ROOT.gROOT.FindObject("dummy2").Delete()
#               histo2 = ROOT.TH1Keys("dummy2","dummy2",int(nb),float(xmin),float(xmax), ROOT.RooKeysPdf.MirrorBoth)
#               nev = tree.Draw("%s>>%s" % (expr,"dummy2"), cut ,"goff")
#               histo = histo2.GetHisto().Clone(name)
#               if histo.Integral() != 0: histo.Scale( newScale / histo.Integral() )
#           return histo.Clone(name)
            if ROOT.gROOT.FindObject("dummy") != None: ROOT.gROOT.FindObject("dummy").Delete()
            if self._options.keysHist and self._weight and self._options.keysMaxEntries == 0:
                histo = ROOT.TH1Keys("dummy","dummy",int(nb),float(xmin),float(xmax), ROOT.RooKeysPdf.MirrorBoth)
                nev = tree.Draw("%s>>%s" % (expr,"dummy"), cut ,"goff")
                return histo.GetHisto().Clone(name)
            else:
                histo = ROOT.TH1F("dummy","dummy",int(nb),float(xmin),float(xmax))
                histo.Sumw2()
                nev = tree.Draw("%s>>%s" % (expr,"dummy"), cut ,"goff")
                if self._weight and self._options.keysHist and nev > 0 and nev < self._options.keysMaxEntries:
                    histo.Delete()
                    histo = ROOT.TH1Keys("dummy","dummy",int(nb),float(xmin),float(xmax), ROOT.RooKeysPdf.MirrorBoth)
                    nev = tree.Draw("%s>>%s" % (expr,"dummy"), cut ,"goff")
                    return histo.GetHisto().Clone(name)
                else:
                    return histo.Clone(name)
    def __str__(self):
        mystr = ""
        mystr += str(self._fname) + '\n'
        mystr += str(self._tfile) + '\n'
        mystr += str(self._chancuts) + '\n'
        mystr += str(self._weight) + '\n'
        mystr += str(self._scaleFactor)
        return mystr

def addTreeToYieldOptions(parser):
    parser.add_option("-l", "--lumi",           dest="lumi",   type="float", default="1.0", help="Luminosity (in 1/fb)");
    parser.add_option("-w", "--weight",         dest="weight",       action="store_true", help="Use weight (in MC events)");
    parser.add_option("-W", "--weightString",   dest="weightString", type="string", default="puW*kfW*baseW*effW*triggW", help="Use weight (in MC events)");
    parser.add_option(      "--mva",            dest="mva",    help="Attach this MVA (e.g. BDT_5ch_160)");
    parser.add_option(      "--eff",            dest="eff",    help="Attach this efficiency tree");
    parser.add_option(      "--keysHist",       dest="keysHist",  action="store_true", default=False, help="Use TH1Keys to make smooth histograms");
    parser.add_option(      "--keysMaxEntries",       dest="keysMaxEntries", type="int", default="0", help="If > 0, will use TH1Keys only if the histogram has less than this number of entries");
    parser.add_option("-i", "--inclusive",  dest="inclusive", action="store_true", help="Only show totals, not each final state separately");
    parser.add_option("--comboppo",  dest="comboppo", action="store_true", help="Include a (em+me) column");
    parser.add_option("-f", "--final",  dest="final", action="store_true", help="Just compute final yield after all cuts");
    parser.add_option("-e", "--errors",  dest="errors", action="store_true", help="Include uncertainties in the reports");
    parser.add_option("-S", "--start-at-cut",   dest="startCut",   type="string", help="Run selection starting at the cut matched by this regexp, included.") 
    parser.add_option("-U", "--up-to-cut",      dest="upToCut",   type="string", help="Run selection only up to the cut matched by this regexp, included.") 
    parser.add_option("-X", "--exclude-cut", dest="cutsToExclude", action="append", default=[], help="Cuts to exclude (regexp matching cut name), can specify multiple times.") 
    parser.add_option("-I", "--invert-cut",  dest="cutsToInvert",  action="append", default=[], help="Cuts to invert (regexp matching cut name), can specify multiple times.") 
    parser.add_option("-R", "--replace-cut", dest="cutsToReplace", action="append", default=[], nargs=3, help="Cuts to invert (regexp of old cut name, new name, new cut); can specify multiple times.") 
    parser.add_option("-A", "--add-cut",     dest="cutsToAdd",     action="append", default=[], nargs=3, help="Cuts to insert (regexp of cut name after which this cut should go, new name, new cut); can specify multiple times.") 
    parser.add_option("-N", "--n-minus-one", dest="nMinusOne", action="store_true", help="Compute n-minus-one yields and plots")
    parser.add_option("-t", "--tree",          dest="tree", default='latino', help="Pattern for tree name");
    parser.add_option("-F", "--flavors",       dest="flavors",  action="store_true", default=False, help="Print the SF and OF");
    parser.add_option("-G", "--no-fractions",  dest="fractions",action="store_false", default=True, help="Don't print the fractions");

def mergeReports(reports):
    import copy
    one = copy.deepcopy(reports[0])
    for i,(c,x) in enumerate(one):
        for j,xj in enumerate(x):
            one[i][1][j][2] = pow(one[i][1][j][2], 2)
    for two in reports[1:]:
        for i,(c,x) in enumerate(two):
            for j,xj in enumerate(x):
                one[i][1][j][1] += xj[1]
                one[i][1][j][2] += pow(xj[2],2)
    for i,(c,x) in enumerate(one):
        for j,xj in enumerate(x):
            one[i][1][j][2] = sqrt(one[i][1][j][2])
    return one

def mergePlots(plots, clone=False):
    one = plots[0]
    if clone: one = [(k,h.Clone()) for k,h in plots[0]]
    for two in plots[1:]:
        for i,(k,h) in enumerate(two): 
            one[i][1].Add(h)
    return one

