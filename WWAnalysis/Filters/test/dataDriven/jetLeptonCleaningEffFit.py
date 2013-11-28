import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitEleID_ICHEP.py <scenario>
### scenarios:
###   - data_all:    will fit tnpJPsi_Data.root with bins suitable for the current data
###   - datalike_mc: will fit tnpJPsi_{JPsiMuMu,ppMuX}_Spring10_0.117pb.root MC but
###                  with same config as data

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TnP_JetCleaning_Ele = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    InputFileNames = cms.vstring('tnpJetCleaning.root'),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeElEl"),
    OutputFileName = cms.string("TnP_Z_JetCleaning_Ele_%s.root" % scenario),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60", "120", "GeV/c^{2}"),
        pt   = cms.vstring("Probe p_{T}",  "0",  "1000", "GeV/c"),
        eta  = cms.vstring("Probe #eta",   "-3", "3", ""),
        nVtx = cms.vstring("N(prim. vtx)", "0", "6", ""),
        nJet     = cms.vstring("N(jets)",        "0", "10", ""),
        nJetNoPU = cms.vstring("N(jets w/o PU)", "0", "10", ""),
    ),

    Cuts = cms.PSet(
        jetVeto     = cms.vstring("Jet Veto", "nJet",     "0.5"),
        jetVetoNoPU = cms.vstring("Jet Veto", "nJetNoPU", "0.5"),
    ),

    Categories = cms.PSet(
        passSafe  = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        passAggro = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        passSafeNoPU  = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        passAggroNoPU = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
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

    ),

    Efficiencies = cms.PSet(), # will be filled later
)

if scenario == "mc_powheg":
    process.TnP_JetCleaning_Ele.InputFileNames = [ "tnpJetCleaning_MC_powheg.root" ]

process.TnP_JetVeto_Ele = process.TnP_JetCleaning_Ele.clone(
    InputDirectoryName = "tpTreeElElOnePair",
    OutputFileName     = "TnP_Z_JetVeto_Ele_%s.root" % scenario,
)


AVERAGE = cms.PSet(
    pt  = cms.vdouble(  20, 100),
    eta = cms.vdouble(-2.5, 2.5),
)
AVERAGE_JV = AVERAGE.clone(nJetNoPU = cms.vdouble(-0.5,0.5))
ETA_BINS = cms.PSet(
    pt     = cms.vdouble(20, 100),
    eta = cms.vdouble(-2.5,-1.5,1.5,2.5)
)
PT_BINS = cms.PSet(
    pt     = cms.vdouble(10, 20, 30, 45, 70, 100),
)
VTX_BINS = cms.PSet(
    nVtx = cms.vdouble(*[x+0.5 for x in range(8)]),
)
VTX_JV = VTX_BINS.clone(nJetNoPU = cms.vdouble(-0.5,0.5))
JET_BINS = cms.PSet(
    nJet = cms.vdouble(*[x-0.5 for x in range(8)]),
)
JET_NOPU_BINS = cms.PSet(
    nJetNoPU = cms.vdouble(*[x-0.5 for x in range(8)]),
)

IDS=["passSafe", "passAggro", "passSafeNoPU", "passAggroNoPU"]
ALLBINS=[("eta",ETA_BINS),("pt",PT_BINS),("vtx",VTX_BINS),("vtxJetVeto",VTX_JV),("jet",JET_BINS),("jetNoPU",JET_NOPU_BINS),("averageJetVeto",AVERAGE_JV),("average",AVERAGE)]
for T in IDS: 
    for  X,B in ALLBINS:
        if T.find("jet") != -1 and X.find("jet") != -1: continue
        if T.find("jet") != -1 and X.find("Jet") != -1: continue
        for S in ["voigtPlusExpo", "vpvPlusExpo"]:
        #for S in [""]:
            M = cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(T,"pass" if T.find("jetVeto") == -1 else "below"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = B,    
            )
            if S != "": setattr(M,"BinToPDFmap", cms.vstring(S))
            setattr(process.TnP_JetCleaning_Ele.Efficiencies, T+"_"+X+"_"+S, M)

process.TnP_JetCleaning_Mu = process.TnP_JetCleaning_Ele.clone(
    InputDirectoryName = "tpTreeMuMu",
    OutputFileName     = "TnP_Z_JetCleaning_Mu_%s.root" % scenario,
)

IDS=["jetVeto","jetVetoNoPU"]
ALLBINS=[("vtx",VTX_BINS),("average",AVERAGE)]
for T in IDS: 
    for  X,B in ALLBINS:
        if T.find("jet") != -1 and X.find("jet") != -1: continue
        if T.find("jet") != -1 and X.find("Jet") != -1: continue
        #for S in ["vpvPlusExpo"]:
        for S in [""]:
            M = cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(T,"pass" if T.find("jetVeto") == -1 else "below"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = B,    
            )
            if S != "": setattr(M,"BinToPDFmap", cms.vstring(S))
            setattr(process.TnP_JetVeto_Ele.Efficiencies, T+"_"+X+"_"+S, M)


process.TnP_JetVeto_Mu = process.TnP_JetVeto_Ele.clone(
    InputDirectoryName = "tpTreeMuMuOnePair",
    OutputFileName     = "TnP_Z_JetVeto_Mu_%s.root" % scenario,
)



#process.pEleC = cms.Path(process.TnP_JetCleaning_Ele)
#process.pMuC  = cms.Path(process.TnP_JetCleaning_Mu)
process.pEleV = cms.Path(process.TnP_JetVeto_Ele)
process.pMuV  = cms.Path(process.TnP_JetVeto_Mu)


