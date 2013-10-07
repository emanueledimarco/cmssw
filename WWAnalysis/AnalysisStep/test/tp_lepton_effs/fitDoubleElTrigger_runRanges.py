
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


if TYPE == 'DATA_2011A' : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_DoubleElectron_Run2011A_Total.root"
elif TYPE == 'DATA_2011B'  : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_DoubleElectron_Run2011B_Total.root"
elif TYPE == 'DATA_2011BC' : INPUTFILE = "/data/gpetrucc/7TeV/tnp/tnpZ_PD_DoubleElectron_Run2011B_Total.crop.root"

if "DATA_2011B" in TYPE and MIN_RUN>10:   MIN_RUN = max(MIN_RUN, 175832)
elif "DATA_2011A" in TYPE and MAX_RUN>10: MAX_RUN = min(MAX_RUN, 175770)

process.TnP_MuonID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    InputFileNames = cms.vstring(INPUTFILE),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeElElBDT"),
    OutputFileName = cms.string("REPLACE_ME"),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "76", "106", "GeV/c^{2}"),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV"),        
        run = cms.vstring("Run range","0","999999",""),
    ),

    Categories = cms.PSet(
        ## Cuts Sequence
        passAllBDT = cms.vstring("passID","dummy[pass=1,fail=0]"),
        passIso    = cms.vstring("passIso","dummy[pass=1,fail=0]"),
        passConvR  = cms.vstring("passIso","dummy[pass=1,fail=0]"),
        passIP     = cms.vstring("passIso","dummy[pass=1,fail=0]"),
        # This is also needed
        passLooseFake = cms.vstring("passIso","dummy[pass=1,fail=0]"),
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
    passAllBDT  = cms.vstring("pass"),
    passIso     = cms.vstring("pass"),
    passIP      = cms.vstring("pass"),
    passConvR   = cms.vstring("pass"),
    passLooseFake = cms.vstring("pass"),
    #pt     = cms.vdouble( 10, 13, 16, 18, 20, 24, 30, 40, 110),
    pt     = cms.vdouble( 10, 11.5, 13, 14, 15, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 23, 24, 25, 26, 28, 30, 40, 110),
    abseta = cms.vdouble( 0.0, 1.5, 2.5),
)

runRanges = [
    #   start    end       pass probe1          den name         den flags   
    #(     1.0, 164237.5, "HLT_DoubleMu7"),
    (100000.5, 170052.5, "HLT_Ele17_Ele8_Ele17Leg",        ),
    (100000.5, 170052.5, "HLT_Ele17_Ele8_Ele8Leg",         ),
    (170052.5, 999999.9, "HLT_Ele17_Ele8_Tighter_Ele17Leg",),
    (170052.5, 999999.9, "HLT_Ele17_Ele8_Tighter_Ele8Leg", ),
]
denominators = [
    #("SC8",  "HLT_Ele17_SC8"),
    ("Ele8", "HLT_Ele17_Ele8"),
    #("SC17", "HLT_Ele32_SC17"),
]
if MIN_RUN < 10: runRanges = [ runRanges[MIN_RUN] ]
if MAX_RUN < 10: denominators = [ denominators[MAX_RUN] ]
for (start, end, trigger) in runRanges:
    if TYPE == "DATA_2011A":
        if start > 175770: continue
        if end   > 175770: end = 175770
    if TYPE == "DATA_2011B":
        if end   < 175832: continue
        if start < 175832: start = 175832
    if ((MIN_RUN > 10 and end < MIN_RUN) or (MAX_RUN > 10 and start > MAX_RUN)): continue
    for (dname, dtrig) in denominators:
        print "Fitting trigger %s from %s in run range [%d, %d]" % (trigger, dname, start+0.5,end-0.5)
        OUTPUTFILE = "ElectronTriggerV2_%s_TnP_Z_%s_from%s_PtEtaBins_GiovanniBins_RunRange%d-%d.root" % (FUNC,trigger,dname,start+0.5,end-0.5)
        TnP_MuonID = process.TnP_MuonID.clone(  OutputFileName = cms.string(OUTPUTFILE) )
        # create flags
        for X in (trigger, "tag_"+dtrig+"_TagLeg", dtrig+"_ProbeLeg"):
            setattr(TnP_MuonID.Categories, X, cms.vstring(X,"dummy[pass=1,fail=0]"))
        # set numerator
        TnP_MuonID.Efficiencies.PASSING_all = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(trigger, "pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = ONE_BIN.clone(
                run         = cms.vdouble(start, end)
            ),
        )
        # set denominator 
        setattr(TnP_MuonID.Efficiencies.PASSING_all.BinnedVariables, "tag_"+dtrig+"_TagLeg",   cms.vstring("pass"))
        setattr(TnP_MuonID.Efficiencies.PASSING_all.BinnedVariables,        dtrig+"_ProbeLeg", cms.vstring("pass"))
        # set fit function
        if FUNC != 'None': 
            TnP_MuonID.Efficiencies.PASSING_all.BinToPDFmap = cms.vstring(FUNC)
        # make analyzer and path
        setattr(process,"TnP_MuonID_%d_%d_%s_%s" %(start+0.5,end-0.5,trigger,dname), TnP_MuonID)
        setattr(process,"p_%d_%d_%s_%s" %(start+0.5,end-0.5,trigger,dname), cms.Path(TnP_MuonID))

