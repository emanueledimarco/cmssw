
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
        ## Cuts Sequence
        passID       = cms.vstring("passID","dummy[pass=1,fail=0]"),
        passID_kink  = cms.vstring("passID_kink","dummy[pass=1,fail=0]"),
        passIP  = cms.vstring("passIP","dummy[pass=1,fail=0]"),
        passIso = cms.vstring("passIso","dummy[pass=1,fail=0]"),

        passTM  = cms.vstring("passTM","dummy[pass=1,fail=0]"),
        passGlb = cms.vstring("passGlb","dummy[pass=1,fail=0]"),

        tag_HLT_SingleMu_LEG   = cms.vstring("tag_HLT_SingleMu_LEG","dummy[pass=1,fail=0]"),
        tag_HLT_Mu17leg_Mu8_bis1 = cms.vstring("tag_HLT_Mu17leg_Mu8_bis1","dummy[pass=1,fail=0]"),
        tag_HLT_Mu13leg_Mu8_bis1 = cms.vstring("tag_HLT_Mu13leg_Mu8_bis1","dummy[pass=1,fail=0]"),

        #HLT_DoubleMu7_bis1   = cms.vstring("HLT_DoubleMu7","dummy[pass=1,fail=0]"),
        HLT_Mu13leg_Mu8_bis1 = cms.vstring("HLT_Mu13leg_Mu8","dummy[pass=1,fail=0]"),
        HLT_Mu13_Mu8leg_bis1 = cms.vstring("HLT_Mu13_Mu8leg","dummy[pass=1,fail=0]"),
        HLT_Mu17leg_Mu8_bis1 = cms.vstring("HLT_Mu17leg_Mu8","dummy[pass=1,fail=0]"),
        HLT_Mu17_Mu8leg_bis1 = cms.vstring("HLT_Mu17_Mu8leg","dummy[pass=1,fail=0]"),
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
    pt     = cms.vdouble( 10, 11.5, 13, 14.5, 16, 18, 20, 30, 40, 110),
    abseta = cms.vdouble( 0.0, 1.2, 2.1, 2.4),
    tag_HLT_SingleMu_LEG = cms.vstring("pass"),
)


runRanges = [
    #   start    end       pass probe1   
    #(     1.0, 164237.5, "HLT_DoubleMu7"),
    (165084.5, 178419.5, "HLT_Mu13leg_Mu8" ),
    (165084.5, 178419.5, "HLT_Mu13_Mu8leg"),
    (178419.5, 999999.9, "HLT_Mu17leg_Mu8"),
    (178419.5, 999999.9, "HLT_Mu17_Mu8leg"),
]
if MIN_RUN < 10: runRanges = [ runRanges[MIN_RUN] ]
for (start, end, trigger) in runRanges:
    if TYPE == "DATA_2011A":
        if start > 175770: continue
        if end   > 175770: end = 175770
    if TYPE == "DATA_2011B":
        if end   < 175832: continue
        if start < 175832: start = 175832
    if ((MIN_RUN > 10 and end < MIN_RUN) or start > MAX_RUN): continue
    OUTPUTFILE = "MuonsTrigger_%s_TnP_Z_%s_PtEtaBins_MuFromTk_GiovanniBins_RunRange%d-%d.root" % (FUNC,trigger,start+0.5,end-0.5)
    TnP_MuonID = process.TnP_MuonID.clone(  OutputFileName = cms.string(OUTPUTFILE) )
    TnP_MuonID.Efficiencies.PASSING_all = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring(trigger+"_bis1", "pass"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = ONE_BIN.clone(
            passIso     = cms.vstring("pass"),
            passID_kink = cms.vstring("pass"),
            passIP      = cms.vstring("pass"),
            run         = cms.vdouble(start, end)
        ),
        BinToPDFmap = cms.vstring(FUNC),
    )
    setattr(process,"TnP_MuonID_%d_%d_%s" %(start+0.5,end-0.5,trigger), TnP_MuonID)
    setattr(process,"p_%d_%d_%s" %(start+0.5,end-0.5,trigger), cms.Path(TnP_MuonID))
