# python WWAnalysis/AnalysisStep/test/zz4l/syncDump.py hzzTree.root --ebe-corr=mc53x >! dump.ext.california.53x.id1125
# python WWAnalysis/AnalysisStep/test/zz4l/syncDump.py hzzTree.root --ebe-corr=mc53x --short >! dump.short.california.53x.id1125
# python WWAnalysis/AnalysisStep/test/zz4l/syncDump.py hzzTree.root --ebe-corr=mc42x >! dump.ext.california.42x.id1125
# python WWAnalysis/AnalysisStep/test/zz4l/syncDump.py hzzTree.root --ebe-corr=mc42x --short >! dump.short.california.42x.id1125


#!/usr/bin/env python
import re, sys, json, string
from sys import stdout, stderr, exit, modules
from optparse import OptionParser
from math import *

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argvbackup = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = argvbackup[:]

parser = OptionParser(usage="usage: %prog [options] rootfile [what] \nrun with --help to get list of options")
parser.add_option("-b", "--blind",   dest="blind",  default=False, action="store_true",  help="blind OS")
parser.add_option("-M", "--max-events", dest="maxevents", default=999999, type="int", nargs=1, help="Maximum number of events to process")
parser.add_option("-m", "--mass-range", dest="massrange", default=(0,9999), type="float", nargs=2, help="Mass range")
parser.add_option("-r", "--run-range",  dest="runrange", default=(0,99999999), type="float", nargs=2, help="Run range")
parser.add_option("--dr", "--min-dr-cut",  dest="minDR", default=0, type="float", help="Run range")
parser.add_option("-t", "--tree",  dest="tree", default=None, type="string", help="Tree to use")
parser.add_option("-p", "--presel-cut",  dest="precut", default=None, type="string", help="Cut to apply before removing duplicates")
parser.add_option("-c", "--cut",  dest="cut", default=None, type="string", help="Cut to apply")
parser.add_option("-T", "--type",  dest="type", default=None, type="string", help="Type of events to select")
parser.add_option("-N", "--nicola",   dest="nicola",  default=False, action="store_true",  help="print Nicola's header row")
parser.add_option("-F", "--fudge",   dest="fudge",  default=False, action="store_true",  help="print -999 for missing variables")
parser.add_option("--id", "--event-id", dest="eventid", default=False, action="store_true", help="print run:lumi:event for skimming")
parser.add_option("-j", "--json",   dest="json",  default=None, type="string", help="JSON file to apply")
parser.add_option("-f", "--format",   dest="fmt",  default=None, type="string",  help="Print this format string")
parser.add_option("-s", "--short",   dest="short",  default=False, action="store_true",  help="print short format")
parser.add_option("--ebe-corr", dest="ebeCorr", default="reco53x", type="string",  help="Add also corrected ev-by-ev mass errors")

(options, args) = parser.parse_args(); sys.argv = []
what = args[1] if len(args) > 1 else "signal"
if what not in [ "signal", "CRss", "CRos" ]: raise RuntimeError, "Unknown what"

jsonmap = {}
if options.json:
    J = json.load(open(options.json, 'r'))
    for r,l in J.iteritems():
        jsonmap[long(r)] = l
    stderr.write("Loaded JSON %s with %d runs\n" % (options.json, len(jsonmap)))

def testJson(ev):
    try:
        lumilist = jsonmap[ev.run]
        for (start,end) in lumilist:
            if start <= ev.lumi and ev.lumi <= end:
                return True
        return False
    except KeyError:
        return False

massErrFile = None; massErrMap = {}
if options.ebeCorr:
    massErrFile = ROOT.TFile("/afs/cern.ch/user/g/gpetrucc/public/ebeOverallCorrections.HCP2012.v1.root")
    massErrMap = dict([(K,massErrFile.Get(K+"_"+options.ebeCorr)) for K in ('mu','el')])

def massErrCorr(ev):
    #period = 'reco53x' if ev.run >= 190000L else 'reco42x' 
    sumErr = ev.pho1massErr**2 + ev.pho2massErr**2
    for i in range(1,5):
        pt  = getattr(ev, 'l%dpt' % i)
        eta = abs(getattr(ev, 'l%deta' % i))
        kind = 'mu' if abs(getattr(ev, 'l%dpdgId'%i)) == 13 else 'el'
        errB = getattr(ev, 'l%dmassErr' % i)
        mapC = massErrMap[kind]
        ptBin  = min(max(1,mapC.GetXaxis().FindBin(pt )), mapC.GetNbinsX());
        etaBin = min(max(1,mapC.GetYaxis().FindBin(eta)), mapC.GetNbinsY());
        sumErr += (errB * mapC.GetBinContent(ptBin,etaBin))**2
    return sqrt(sumErr)

def deltaR(eta1,phi1,eta2,phi2):
    dphi=phi1-phi2;
    while dphi >  pi: dphi -= pi
    while dphi < -pi: dphi += pi
    return hypot(dphi, eta1-eta2)
def mllCut4(ev):
    if (ev.z1mll <= 4 or ev.z2mll <= 4): return False
    if ev.qll13 == 0 and ev.mll13 <= 4: return False
    if ev.qll14 == 0 and ev.mll14 <= 4: return False
    if ev.qll23 == 0 and ev.mll23 <= 4: return False
    if ev.qll24 == 0 and ev.mll24 <= 4: return False
    return True
def minDeltaRZZ(ev):
    return min([deltaR(ev.l1eta,ev.l1phi,ev.l3eta,ev.l3phi),
                deltaR(ev.l2eta,ev.l2phi,ev.l3eta,ev.l3phi),
                deltaR(ev.l2eta,ev.l2phi,ev.l4eta,ev.l4phi),
                deltaR(ev.l1eta,ev.l1phi,ev.l4eta,ev.l4phi)])
class Event:
    def __init__(self,tree,index):
        self.tree = tree
        self.index = index
        self.tree.GetEntry(index); self.tree.index = index
        self.id = (tree.run, tree.lumi, tree.event)
        self.type = "4l"
        if abs(self.l1pdgId) == 13 and abs(self.l3pdgId) == 13: self.type = "4mu"
        if abs(self.l1pdgId) == 13 and abs(self.l3pdgId) == 11: self.type = "2mu2e"
        if abs(self.l1pdgId) == 11 and abs(self.l3pdgId) == 13: self.type = "2e2mu"
        if abs(self.l1pdgId) == 11 and abs(self.l3pdgId) == 11: self.type = "4e"
        self.leptype = self.type
        if self.pho1pt > 0 and self.pho2pt > 0: self.type += "+2g"
        elif self.pho1pt + self.pho2pt > 0:     self.type += "+g"
        self.massErrCorr = massErrCorr(self) if options.ebeCorr else 0
        self.etajj = abs(self.jet1eta-self.jet2eta)
        self.fishjj = 1.92e-4 * self.mjj + 0.18 * self.etajj
        ## zero-suppressed variables
        self.Jet1pt = self.jet1pt if self.njets30 >= 1 else -1.
        self.Jet2pt = self.jet2pt if self.njets30 >= 2 else -1.
        self.mJJ    = self.mjj    if self.njets30 >= 2 else -1.
        self.etaJJ  = self.etajj  if self.njets30 >= 2 else -1.
        self.fishJJ = self.fishjj if self.njets30 >= 2 else -1.
    def __getattr__(self,attr):
        if attr in self.__dict__:
            return self.__dict__[attr]
        global options
        if self.tree.index != self.index: 
            self.tree.GetEntry(self.index)
            self.tree.index = self.index
        if options.fudge:
            try:
                return getattr(self.tree, attr)
            except AttributeError:
                return -9.99
        else:
            return getattr(self.tree, attr)
    def __getitem__(self,attr):
        return self.__getattr__(attr)
    def __lt__(self,other): return self.id <  other.id
    def __le__(self,other): return self.id <= other.id
    def __gt__(self,other): return self.id >  other.id
    def __ge__(self,other): return self.id >= other.id

class BaseDumper:
    def __init__(self,options=None):
        self.options = options
        #self.format = string.Template(options.fmt) if options.fmt else None
        #self.format = string.Formatter(options.fmt) if options.fmt else None
    def preselect(self,ev):
        return True
    def accept(self,ev):
        return True
    def __call__(self,ev):
        if not self.accept(ev): return
        if self.options.fmt: 
            print string.Formatter().vformat(options.fmt,[],ev)
            #print self.format.substitute(ev)
            return

        #compute lepton isolation
        l1iso = ev.l1pfIsoComb04EACorr/ev.l1pt
        if (abs(ev.l1pdgId) == 13) : l1iso = ev.l1pfCombRelIso04dBCorr
        l2iso = ev.l2pfIsoComb04EACorr/ev.l2pt
        if (abs(ev.l2pdgId) == 13) : l2iso = ev.l2pfCombRelIso04dBCorr
        l3iso = ev.l3pfIsoComb04EACorr/ev.l3pt
        if (abs(ev.l3pdgId) == 13) : l3iso = ev.l3pfCombRelIso04dBCorr
        l4iso = ev.l4pfIsoComb04EACorr/ev.l4pt
        if (abs(ev.l4pdgId) == 13) : l4iso = ev.l4pfCombRelIso04dBCorr

        if options.short : 
            print "%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.2f:%.0f:%.2f:%.2f:%.2f:%.3f:%.3f" % (ev.run, ev.lumi, ev.event, ev.mass, ev.z1mass, ev.z2mass, ev.massErr, ev.massErrCorr, ev.ME_SMH_ZZ, ev.pt, ev.njets30, ev.Jet1pt, ev.Jet2pt, ev.mJJ, ev.etaJJ, ev.fishJJ)
        else :
            print "%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.2f:%.1f:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f" % (ev.run, ev.lumi, ev.event, ev.mass, ev.z1mass, ev.z2mass, ev.massErr, ev.massErrCorr, ev.ME_SMH_ZZ, ev.pt, ev.njets30,  ev.Jet1pt, ev.Jet2pt, ev.mJJ, ev.etaJJ, ev.fishJJ, ev.ME_SMH_0minus, ev.ME_SMH_0hplus, ev.ME_SMH_1plus, ev.ME_SMH_1minus, ev.ME_SMH_2mplus_gg, ev.ME_SMH_2mplus_qq)

class SignalDumper(BaseDumper):
    def __init__(self,options=None):
        BaseDumper.__init__(self,options)
    def preselect(self,ev):
        if ev.run < self.options.runrange[0] or ev.run > self.options.runrange[1]: return False
        if options.precut and not eval(options.precut, globals(), {'ev':ev}): return False
        return True
    def accept(self,ev):
        if ev.mass < self.options.massrange[0] or ev.mass > self.options.massrange[1]: return False
        if options.blind: 
            if (110 <= ev.mass and ev.mass <= 140) or (ev.mass > 300): return False
        if options.type and (options.type not in ev.type): return False
        if options.cut and not eval(options.cut, globals(), {'ev':ev}): return False
        return True

class ControlDumper(BaseDumper):
    def __init__(self,what,options=None):
        BaseDumper.__init__(self,options)
        self.what = what
    def preselect(self,ev):
        if ev.run < self.options.runrange[0] or ev.run > self.options.runrange[1]: return False
        if options.precut and not eval(options.precut, globals(), {'ev':ev}): return False
        if "CRss" in self.what and ev.l3pdgId != ev.l4pdgId: return False
        if "CRos" in self.what and ev.l3pdgId == ev.l4pdgId: return False
        ev.minDRZZ = minDeltaRZZ(ev)
        if options.minDR > 0 and ev.minDRZZ < options.minDR: return False
        return True
    def accept(self,ev):
        if ev.mass < self.options.massrange[0] or ev.mass > self.options.massrange[1]: return False
        if ev.z2mass <= 12 or ev.z2mass >= 120: return False
        if not mllCut4(ev): return False
        pts = [ ev.l1pt, ev.l2pt, ev.l3pt, ev.l4pt ]; pts.sort()
        if pts[0] <= 20 or pts[1] <= 10: return False
        if options.blind and self.what == "CRos": 
            if (110 <= ev.mass and ev.mass <= 140) or (ev.mass > 300): return False
        if options.type and (options.type not in ev.type): return False
        if options.cut and not eval(options.cut, globals(), {'ev':ev}): return False
        return True

file = ROOT.TFile.Open(args[0])
treename = "zz4lTree/probe_tree" if what == "signal" else "anyZxTree/probe_tree"
if options.tree: treename = options.tree
tree = file.Get(treename)
dump = SignalDumper(options) if what == "signal" else ControlDumper(what,options)
events = {}
for i in xrange(tree.GetEntries()):
    if i > options.maxevents: break
    ev = Event(tree,i)
    if options.json and not testJson(ev): continue
    if dump.preselect(ev):
        if ev.id in events: 
            ov = events[ev.id]
            if ev.l3pt + ev.l4pt > ov.l3pt + ov.l4pt: 
                events[ev.id] = ev
        else: 
            events[ev.id] = ev
evsorted = events.values(); evsorted.sort()
for ev in evsorted: dump(ev)
