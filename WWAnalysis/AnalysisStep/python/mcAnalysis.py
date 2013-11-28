#!/usr/bin/env python
from WWAnalysis.AnalysisStep.tree2yield import *

class MCAnalysis:
    def __init__(self,samples,options):
        self._foutName = options.out if options.out else samples.replace(".txt","")+".root"
        self._fout     = None
        self._options = options
        self._allData     = {}
        self._data        = []
        self._signals     = []
        self._backgrounds = [] 
        self._isSignal    = {}
        for line in open(samples,'r'):
            if re.match("\s*#.*", line): continue
            field = [f.strip() for f in line.split(':')]
            if field[0][0] == '#': continue
            rootfile = options.path+"/latino_%s.root" % field[1].strip()
            signal = ("%d" in rootfile)
            if field[0][-1] == "+": 
                signal = True
                field[0] = field[0][:-1]
            if field[0][-1] == "-": 
                signal = False
                field[0] = field[0][:-1]
            if ("%d" in rootfile): rootfile = rootfile.replace('%d',str(options.mass))
            if ("%d" in field[0]): field[0] = field[0].replace('%d',str(options.mass))
            tty = TreeToYield(rootfile, options)
            if len(field) == 3: tty.setScaleFactor(field[2])
            if signal: 
                self._signals.append(tty)
                self._isSignal[field[0]] = True
            elif field[0] == "data":
                self._data.append(tty)
            else:
                self._isSignal[field[0]] = False
                self._backgrounds.append(tty)
            if field[0] in self._allData: self._allData[field[0]].append(tty)
            else                        : self._allData[field[0]] =     [tty]
        #if len(self._signals) == 0: raise RuntimeError, "No signals!"
        #if len(self._backgrounds) == 0: raise RuntimeError, "No backgrounds!"
    def listProcesses(self):
        return self._allData.keys()[:]
    def listSignals(self):
        return [ p for p in self._allData.keys() if p != 'data' and self._isSignal[p] ]
    def listBackgrounds(self):
        return [ p for p in self._allData.keys() if p != 'data' and not self._isSignal[p] ]
    def hasProcess(self,process):
        return process in self._allData
    def scaleProcess(self,process,scaleFactor):
        for tty in self._allData[process]: tty.setScaleFactor(scaleFactor)
    def scaleUpProcess(self,process,scaleFactor):
        for tty in self._allData[process]: 
            tty.setScaleFactor( "((%s) * (%s))" % (tty.getScaleFactor(),scaleFactor) )
    # assumes that all the processes have the same SF, oops
    def getScale(self,process):
        for tty in self._allData[process]: return tty.getScaleFactor()
    def getScales(self,process):
        return [ tty.getScaleFactor() for tty in self._allData[process] ] 
    def attachMVA(self,name):
        for p, ttys in self._allData.iteritems(): 
            for tty in ttys: tty.attachMVA(name)
    def attachEff(self,name):
        for p, ttys in self._allData.iteritems(): 
            for tty in ttys: tty.attachEff(name)
    def getYields(self,cuts,process=None,subprocess=None,nodata=False,makeSummary=False,subprocesses=None):
        ret = { }
        allSig = []; allBg = []
        for key in self._allData:
            if key == 'data' and nodata: continue
            if process != None and key != process: continue
            if subprocess != None and sum([tty.name() == subprocess for tty in self._allData[key]]) == 0: continue
            if subprocesses:
                reports = [ (tty.name(), self._getYields([tty],cuts,subprocess)) for tty in self._allData[key] ]
                ret[key] = mergeReports([y for x,y in reports])
                for name, rep in reports: ret['.'+name] = rep
            else:
                ret[key] = self._getYields(self._allData[key],cuts,subprocess)
            if key != 'data':
                if self._isSignal[key]: allSig.append(ret[key])
                else: allBg.append(ret[key])
        if makeSummary:
            if self._signals and not ret.has_key('signal') and len(allSig) > 0:
                ret['signal'] = mergeReports(allSig)
            if self._backgrounds and not ret.has_key('background') and len(allBg) > 0:
                ret['background'] = mergeReports(allBg)
        return ret
    def prettyPrint(self,reports,makeSummary=False):
        allSig = []; allBg = []
        for key in reports:
            if key == 'data': continue
            print "\n ==== {0} ====".format(key)
            self._allData[key][0].prettyPrint(reports[key])
            if key != 'data':
                if self._isSignal[key]: allSig.append(reports[key])
                else: allBg.append(reports[key])
        if makeSummary:
            if len(allSig)>0:
                print "\n ==== ALL SIGNALS ==== "
                self._signals[0].prettyPrint(mergeReports(allSig))
            if len(allBg)>0:
                print "\n ==== ALL BACKGROUNDS ==== "
                self._backgrounds[0].prettyPrint(mergeReports(allBg))
        if self._allData.has_key('data') and reports.has_key('data'):
            print "\n ==== DATA ==== "
            self._allData['data'][0].prettyPrint(reports['data'])
    def getPlots(self,plots,cut,makeSummary=False):
        for (name,expr,bins,xlabel,ylabel) in plots.plots():
            self.getPlotsForCut(name,expr,bins,cut,makeSummary=makeSummary)
    def getPlotsForCut(self,name,expr,bins,cut,write=True,process=None,nodata=False,makeSummary=False):
        report = {}
        allSig = []; allBg = []
        for key in self._allData:
            if key == 'data' and nodata: continue
            if process != None and key not in process: continue
            report[key] = self._getPlots(name,expr,bins,cut,self._allData[key])
            for (k,h) in report[key]:
                if write: self._fOut(key).WriteTObject(h)
            if key != 'data':
                if self._isSignal[key]: allSig.append(report[key])
                else: allBg.append(report[key])
        if makeSummary:
            if self._signals and not report.has_key('signal') and len(allSig)>0:
                report['signal'] = mergePlots(allSig)
            if self._backgrounds and not report.has_key('background') and len(allBg)>0:
                report['background'] = mergePlots(allBg)
        return report;
    def dumpEvents(self,cut,vars=['run','lumi','event','dataset','channel']):
        for tty in self._data: tty.dumpEvents(cut,vars)
    def _getYields(self,ttylist,cuts,subprocess=None):
        return mergeReports([tty.getYields(cuts) for tty in ttylist if (subprocess == None or tty.name() == subprocess)])
    def _getPlots(self,name,expr,bins,cut,ttylist):
        return mergePlots([tty.getPlots(name,expr,bins,cut) for tty in ttylist])
    def _fOutName(self):
            return self._foutName;
    def _fOut(self,dir=None):
        if dir == None:
            if not self._fout: self._fout = ROOT.TFile.Open(self._foutName, "RECREATE")
            return self._fout
        else:
            self._fOut()
            tdir = self._fout.GetDirectory(dir)
            if tdir: return tdir
            else:    return self._fout.mkdir(dir)
    def _fOut(self,dir=None):
        if dir == None:
            if not self._fout: self._fout = ROOT.TFile.Open(self._foutName, "RECREATE")
            return self._fout
        elif "/" in dir:
            pieces = dir.split("/")
            base, last = ("/".join(dir[:-1]), dir[-1])
            basedir = self._fOut(base)
            tdir = basedir.GetDirectory(dir)
            if tdir: return tdir
            else:    return basedir.mkdir(dir)
        else:
            self._fOut()
            tdir = self._fout.GetDirectory(dir)
            if tdir: return tdir
            else:    return self._fout.mkdir(dir)
    def __str__(self):
        mystr = ""
        for a in self._allData:
            mystr += str(a) + '\n' 
        for a in self._data:
            mystr += str(a) + '\n' 
        for a in self._signals:
            mystr += str(a) + '\n' 
        for a in self._backgrounds:
            mystr += str(a) + '\n'
        return mystr[:-1]

def addMCAnalysisOptions(parser,addTreeToYieldOnesToo=True):
    if addTreeToYieldOnesToo: addTreeToYieldOptions(parser)
    parser.add_option("-o", "--out",    dest="out",     help="Output file name. by default equal to input -'.txt' +'.root'");
    parser.add_option("-m", "--mass",   dest="mass",    type="int",    default="160", help="Higgs boson mass");
    parser.add_option("-P", "--path",   dest="path",    type="string", default="./",      help="path to directory with trees (./)") 

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree.root cuts.txt")
    addMCAnalysisOptions(parser)
    parser.add_option("-D", "--dump",   dest="dump",    action="store_true", help="Dump events passing selection");
    parser.add_option("-p", "--plots",  dest="plots",   type="string", metavar="FILE", help="Make the plots defined in plot file");
    parser.add_option("--process",      dest="process", type="string", default=None, help="Process to print out (default = all)");
    (options, args) = parser.parse_args()
    tty = TreeToYield(args[0],options) if ".root" in args[0] else MCAnalysis(args[0],options)
    cf  = CutsFile(args[1],options)
    for cutFile in args[2:]:
        temp = CutsFile(cutFile,options)
        for cut in temp.cuts():
            cf.add(cut[0],cut[1])
    if options.plots:
        pf = PlotsFile(options.plots, options)
        tty.getPlots(pf, cf.allCuts()) #, process=options.process)
    else:
        report = tty.getYields(cf)#, process=options.process)
        tty.prettyPrint(report)
    #tty.getPlots("gammaMRStar","gammaMRStar","200,0.,200.",cf.allCuts())
    if options.dump: tty.dumpEvents(cf.allCuts())
