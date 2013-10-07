
import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
TYPE = args[0]
FUNC = args[1]
MIN_RUN = int(args[2]) if len(args) > 2 else 1
MAX_RUN = int(args[3]) if len(args) > 3 else 999999

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )


if TYPE   == 'MC'         : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_Zmumu_42X_MC_preTotal_HALF.root"
if TYPE   == 'MC41X'      : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_Zmumu_MC_41X.root"
elif TYPE == 'DATA_2011A' : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_SingleMu_Run2011A_Total.root"
elif TYPE == 'DATA_2011B' : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_SingleMu_Run2011B_Total.root"
elif TYPE == 'May10' : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_SingleMu_Run2011A-May10_Total.root"
elif TYPE == 'PRv4'  : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_SingleMu_Run2011A-v4_Total.root"
elif TYPE == '05Aug' : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_SingleMu_Run2011A-Aug05_Total.root"
elif TYPE == 'PRv6'  : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_SingleMu_Run2011A-v6_Total.root"
elif TYPE == 'PRv6v2': INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_SingleMu_Run2011A-v6v2_Total.root"

if TYPE == "DATA_2011B" and MIN_RUN>10:   MIN_RUN = max(MIN_RUN,175832)
elif TYPE == "DATA_2011A" and MIN_RUN>10: MAX_RUN = min(MAX_RUN, 175770)

## From Thea
singleMuDataPaths = cms.vstring(
    "1-163261:HLT_Mu15_v*",
    "163262-164237:HLT_Mu24_v*",
    "165085-173235:HLT_Mu30_v*",
    "173236-175972:HLT_Mu40_v*",
    "175973-999999:HLT_Mu40_eta2p1_v*",
    "163262-167043:HLT_IsoMu17_v*",
    "167044-167913:HLT_IsoMu17_eta2p1_v*",
    "167914-175910:HLT_IsoMu20_v*",
    "175911-175921:HLT_IsoMu24_v*",
    "175922-176544:HLT_IsoMu24_eta2p1_v*",
    "176545-999999:HLT_IsoMu30_eta2p1_v*",
)
## Json to ignore dummy runs
import json
bigjson = json.load(open("Cert_160404-179431_7TeV_PromptReco_Collisions11_JSON.txt","r"))
## Create set of run ranges and triggers
triggersWithRuns = []; allendpoints = []
for x in singleMuDataPaths:
    (start,end,trigger) = x.replace('-',':').replace("_v*","").split(":")
    triggersWithRuns.append([int(start), int(end), trigger, "Iso" in trigger, "eta2p1" in trigger])
    allendpoints.append([int(start), len(triggersWithRuns)-1, True])
    allendpoints.append([int(end),   len(triggersWithRuns)-1, False])
allendpoints.sort()
runRanges = []; lastrun = allendpoints[0][0]; wasstart = True;
for endpoint, index, starting in allendpoints[1:]:
    if endpoint <= lastrun+1: continue
    found = False
    for run in range(lastrun,endpoint+1):
        if str(run) in bigjson:
            found = True;
            break
    if not found:
        print "Run range [%d, %d] does not contain collision runs." % (lastrun,endpoint)
        continue
    center = 0.5*(endpoint + lastrun)
    triggers = []
    for (start,end,trigger,iso,barrel) in triggersWithRuns:
        if center >= start and center <= end: triggers.append(trigger)
    if len(triggers):
        runRanges.append((lastrun +0.5-wasstart , endpoint  +0.5-starting, triggers))
    lastrun = endpoint; wasstart = starting
from pprint import pprint
pprint(runRanges)
#raise RuntimeError, "STC"

process.TnP_MuonID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    InputFileNames = cms.vstring(INPUTFILE),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeMuMu"),
    OutputFileName = cms.string("REPLACE_ME"),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "76", "106", "GeV/c^{2}"),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV"),        
        run = cms.vstring("Run range","0","999999",""),
    ),

    Categories = cms.PSet(
        passID       = cms.vstring("passID","dummy[pass=1,fail=0]"),
        passID_kink  = cms.vstring("passID_kink","dummy[pass=1,fail=0]"),
        passIP  = cms.vstring("passIP","dummy[pass=1,fail=0]"),
        passTM  = cms.vstring("passTM","dummy[pass=1,fail=0]"),
        passIso = cms.vstring("passIso","dummy[pass=1,fail=0]"),
        tag_HLT_SingleMu_LEG   = cms.vstring("tag_HLT_SingleMu_LEG","dummy[pass=1,fail=0]"),
    ),
    PDFs = cms.PSet(
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(20),
    Efficiencies = cms.PSet(), # will be filled later
)

ONE_BIN = cms.PSet(
    pt     = cms.vdouble( 10, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 40, 110),
    abseta = cms.vdouble( 0.0, 1.2, 2.1, 2.4),
    tag_HLT_SingleMu_LEG = cms.vstring("pass"),
)

if MIN_RUN < 10: runRanges = [ runRanges[MIN_RUN] ]
for (start,end,triggers) in runRanges:
    if ((MIN_RUN > 10 and end < MIN_RUN) or start > MAX_RUN): continue
    OUTPUTFILE = "MuonsTrigger_%s_TnP_Z_%s_SingleMu_PtEtaBins_MuFromTk_GiovanniBins_RunRange%d-%d.root" % (FUNC,TYPE,start+0.5,end-0.5)
    TnP_MuonID = process.TnP_MuonID.clone(  OutputFileName = cms.string(OUTPUTFILE) )
    EXPR = "+".join(triggers)
    for t in triggers: setattr(TnP_MuonID.Categories, t, cms.vstring(t,"dummy[pass=1,fail=0]"))
    TnP_MuonID.Efficiencies.PASSING_all = cms.PSet(
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = ONE_BIN.clone(
            passIso     = cms.vstring("pass"),
            passID_kink = cms.vstring("pass"),
            passIP      = cms.vstring("pass"),
            passTM      = cms.vstring("pass"),
            run         = cms.vdouble(start, end)

        ),
        BinToPDFmap = cms.vstring(FUNC),
    )
    if triggers == [ "HLT_Mu15" ]: # fix tag trigger
        del TnP_MuonID.Efficiencies.PASSING_all.BinnedVariables.tag_HLT_SingleMu_LEG
        TnP_MuonID.Categories.HLT_Mu15_bis1 = cms.vstring("Mu15","dummy[pass=1,fail=0]")
        TnP_MuonID.Categories.tag_HLT_Mu15_bis1 = cms.vstring("Mu15","dummy[pass=1,fail=0]")
        TnP_MuonID.Efficiencies.PASSING_all.BinnedVariables.tag_HLT_Mu15_bis1 = cms.vstring("pass")
        triggers = [ "HLT_Mu15_bis1" ]
    elif triggers == [ "HLT_Mu24", "HLT_IsoMu17" ]:
        del TnP_MuonID.Efficiencies.PASSING_all.BinnedVariables.tag_HLT_SingleMu_LEG
        TnP_MuonID.Categories.tag_HLT_Mu24 = cms.vstring("Mu15","dummy[pass=1,fail=0]")
        TnP_MuonID.Efficiencies.PASSING_all.BinnedVariables.tag_HLT_Mu24 = cms.vstring("pass")
    if len(triggers) > 1:
        TnP_MuonID.Expressions = cms.PSet(theOR = cms.vstring("theOR", EXPR, *triggers))
        TnP_MuonID.Cuts        = cms.PSet(theORcut = cms.vstring("theORcut", "theOR", "0.5"))
        TnP_MuonID.Efficiencies.PASSING_all.EfficiencyCategoryAndState = cms.vstring("theORcut", "above")
    else:
        TnP_MuonID.Efficiencies.PASSING_all.EfficiencyCategoryAndState = cms.vstring(triggers[0], "pass")
    print "Will process range (%d, %d) for %s" % (start, end, triggers)
    setattr(process,"TnP_MuonID_%d_%d" %(start+0.5,end-0.5), TnP_MuonID)
    setattr(process,"p_%d_%d" %(start+0.5,end-0.5), cms.Path(TnP_MuonID))
