import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

choiceID=None; postfix=""
if len(args) > 1:
    choiceID = args[1]
    postfix  = "."+choiceID


process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TnP_EleID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    InputFileNames = cms.vstring('tnpZ_Data.root'),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeElEl"),
    OutputFileName = cms.string("TnP_Z_EleID_%s%s.root" % (scenario,postfix)),

    Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "70", "120", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe #eta", "-3", "3", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "3", ""),
        nVtx = cms.vstring("N(prim. vtx)", "0", "6", ""),
        pair_nJet = cms.vstring("N(jets)",        "0", "10", ""),
        run    = cms.vstring("Run number", "0", "9999999", ""),
    ),

    Categories = cms.PSet(
        passID    = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        passIso   = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        passIP    = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        passConvR = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        HLT_Ele8_CaloIdL_CaloIsoVL = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[2,1,3])",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        emanuelesMagicFunction = cms.vstring(
            ("EXPR::signal('(@1<@0)*exp(-(@0-@1)*(@0-@1)/(@2*@2 + @3*(@0-@1)*(@0-@1))) + (@1>=@0)*exp(-(@0-@1)*(@0-@1)/(@4*@4 + @5*(@0-@1)*(@0-@1)))',"+
                           "mean[91.2, 80.0, 100.0],mass, sigmaL[2.3, 0.5, 10.0],alphaL[0.23],sigmaR[2.3, 0.5, 10.0],alphaR[0.2,0,3])"),
            "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
            "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
            "efficiency[0.5,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),

    Efficiencies = cms.PSet(), # will be filled later
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(80),
)

if scenario == "data_all":
    process.TnP_EleID.InputFileNames = [ "tnpEle2011Av1.root", "tnpEle2011Av2.root" ]
elif scenario == "mc":
    process.TnP_EleID.InputFileNames = [ "tnpEleMC.root" ]

CONSTRAINTS = cms.PSet(
    pt        = cms.vdouble(20, 100),
    pair_nJet = cms.vdouble(-0.5,0.5)
)
AVERAGE = cms.PSet(
    CONSTRAINTS,
    abseta = cms.vdouble( 0,1.479, 2.5)
)
PT_ETA_BINS = AVERAGE.clone( 
    pt = cms.vdouble(15, 20, 30, 40, 50, 100)
)
VTX_BINS = AVERAGE.clone(
    nVtx = cms.vdouble(*[2*x-0.5 for x in range(8)])
)
JET_BINS = AVERAGE.clone(
    pair_nJet = cms.vdouble(*[x-0.5 for x in range(8)])
)

ETA_FINE_BINS   = cms.PSet(CONSTRAINTS, eta = cms.vdouble(-2.5,-2.0,-1.479,-0.7,0,0.7,1.479,2.0,2.5))

EFFS = {  #Name:NUM:DEN
   'allid':(['passID','passIso','passConvR','passIP'],[]), 
   'idIpConv':(['passID','passIP','passConvR'],[]),  
   'iso':(['passIso'],['passID','passIP','passConvR']),  
   'hlt':(['HLT_Ele8_CaloIdL_CaloIsoVL'],['passID','passIso','passConvR','passIP']), # HLT
}
ALLBINS = []
ALLBINS += [("pt_eta",PT_ETA_BINS)]
ALLBINS += [("average",AVERAGE)]
ALLBINS += [("vtx",VTX_BINS)]
ALLBINS += [("njet", JET_BINS)]
ALLBINS += [("eta", ETA_FINE_BINS)]
#if "data" in scenario: ALLBINS += [ ("run", RUN_BINS) ]

for EN,E in EFFS.items(): 
    if choiceID != None and EN != choiceID: continue
    for  BN,BV in ALLBINS:
        # Num
        NUM = cms.vstring()
        for X in E[0]: NUM += [ X, "pass" ]
        # Den
        DEN = BV.clone()
        for X in E[1]: setattr(DEN, X, cms.vstring("pass"))
        # Shape
        #S = "voigtPlusExpo" if "pt_eta" in BN else "vpvPlusExpo"
        S = "emanuelesMagicFunction"
        # Eff
        setattr(process.TnP_EleID.Efficiencies, EN+"_"+BN, cms.PSet(
            EfficiencyCategoryAndState = NUM,
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = DEN,
            BinToPDFmap = cms.vstring(S)
        ))
        if "fast" in scenario: 
            del getattr(process.TnP_EleID.Efficiencies, EN+"_"+BN).BinToPDFmap

process.p = cms.Path(process.TnP_EleID)

