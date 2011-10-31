import FWCore.ParameterSet.Config as cms
import re
import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]

process = cms.Process("STEP3")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'file:/nfs/bluearc/group/hww/S2/R42X_S1_V06_S2_V02_S3_V05/DYMuMu.root'
    ),
)
process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("WWAnalysis.AnalysisStep.step3_cff")
from WWAnalysis.AnalysisStep.step3_cff import * # get also functions

if len(args) == 0: args = [ 'vbfToH160toWWto2L2Nu', 101160, 0.003621062529384, 'false']
if len(args) != 4: raise RuntimeError, "step3.py dataset id json (for data) or step3.py dataset id scalefactor (for MC)"
## step3.py dataset id json   for data
## step3.py dataset id scalef for MC
dataset = ['MC','ggToH160toWWto2L2Nu']; id = 101160; 
scalef  = 0.003621062529384;
json    = None
mhiggs  = 0
from WWAnalysis.AnalysisStep.fourthScaleFactors_cff import *
fourthGenSF = 1
fermiSF = 1
puStudy = False ## set to true to add 16, yes 16 different PU possibilities
IsoStudy = False ## Set to True to get isolation variables (and a tree build only after ID+CONV+IP, without isolation)
                 ## Note: works only if running also the step2
# from WWAnalysis.AnalysisStep.scaleFactors_cff import *
# if args[1] in dataSamples or args[1] in data42xSamples:
if args[0].find('2011') != -1: args[0] = args[0][ : args[0].find('2011') ]
if args[0] in [ 'SingleElectron', 'DoubleElectron', 'SingleMuon', 'DoubleMuon', 'MuEG']:
    dataset = [args[0]]; id = args[1]
    json    = args[2]
    scalef = 1
else:
    dataset = ['MC', args[0]]; id = args[1];
    scalef  = float(args[2])
    m = re.match("ggToH(\\d+)to.*", args[0])
    n = re.match("vbfToH(\\d+)to.*", args[0])
    if m: 
        mhiggs = int(m.group(1))
        fourthGenSF = fourthGenScales[int(m.group(1))]
        fermiSF = 0
    elif n: 
        mhiggs = -1*int(n.group(1))
        fermiSF = fermiPhobicScales[int(n.group(1))]
process.step3Tree.cut = process.step3Tree.cut.value().replace("DATASET", dataset[0])
process.step3Tree.variables.trigger  = process.step3Tree.variables.trigger.value().replace("DATASET",dataset[0])
process.step3Tree.variables.dataset = str(id)

if dataset[0] == "MC":
#     process.step3Tree.eventWeight = cms.InputTag("mcWeight");
#     process.mcWeight.baseW= scalef
    process.step3Tree.variables.baseW = "%.12f" % scalef
    if mhiggs != 0:
        process.step3Tree.variables.fourW = "%.12f" % fourthGenSF
        process.step3Tree.variables.fermiW = "%.12f" % fermiSF
    else:
        process.step3Tree.variables.fourW = "1"
        process.step3Tree.variables.fermiW = "1"
    if mhiggs <=0:
        process.step3Tree.variables.kfW = cms.string("1")
    else:
        process.higgsPt.inputFilename = "HiggsAnalysis/HiggsToWW2Leptons/data/kfactors_Std/kfactors_mh%(mass)d_ren%(mass)d_fac%(mass)d.dat" % {"mass":abs(mhiggs)}
else:
    from FWCore.PythonUtilities.LumiList import LumiList
    import os    
    lumis = LumiList(filename = os.getenv('CMSSW_BASE')+'/src/WWAnalysis/Misc/Jsons/%s.json'%json)
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess = lumis.getCMSSWString().split(',')
    process.step3Tree.variables.baseW = "1"
    process.step3Tree.variables.fourW = "1"
    process.step3Tree.variables.fermiW = "1"
    process.step3Tree.variables.kfW = cms.string("1")
    process.step3Tree.variables.puW = cms.string("1")


# process.schedule = cms.Schedule()

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
# Scenario 1:
label = "IPMerge"; muon = "wwMuonsMergeIP"; ele = "wwEleIPMerge"; softmu = "wwMuons4Veto"; preSeq = cms.Sequence();
# label = "Scenario1"; muon = "wwMuScenario1"; ele = "wwEleScenario1"; softmu = "wwMu4VetoScenario1"; preSeq = cms.Sequence();
# Scenario 2-5: 
# label = "Scenario2"; muon = "wwMuScenario2"; ele = "wwEleScenario2"; softmu = "wwMu4VetoScenario2"; preSeq = cms.Sequence();
# label = "Scenario3"; muon = "wwMuScenario3"; ele = "wwEleScenario3"; softmu = "wwMu4VetoScenario3"; preSeq = cms.Sequence();
# label = "Scenario4"; muon = "wwMuScenario4"; ele = "wwEleScenario4"; softmu = "wwMu4VetoScenario4"; preSeq = cms.Sequence(process.eleBDTSelection);
# label = "Scenario5"; muon = "wwMuScenario5"; ele = "wwEleScenario5"; softmu = "wwMu4VetoScenario5"; preSeq = cms.Sequence();
if args[3] == 'True' or args[3] == 'true': 
    from WWAnalysis.AnalysisStep.skimEventProducer_cfi import addEventHypothesis
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","HLT")
    addEventHypothesis(process,label,muon,ele,softmu,preSeq)


for X in "elel", "mumu", "elmu", "muel":
    tree = process.step3Tree.clone(src = cms.InputTag("ww%s%s"% (X,label) ));
    seq = cms.Sequence()
    setattr(process, X+"Nvtx", process.nverticesModule.clone(probes = cms.InputTag("ww%s%s"% (X,label))))
    seq += getattr(process, X+"Nvtx")
    tree.variables.nvtx = cms.InputTag(X+"Nvtx")
    if IsoStudy: addIsoStudyVariables(process,tree)
    if dataset[0] == 'MC':
        setattr(process, X+"PuWeight",  process.puWeight.clone(src = cms.InputTag("ww%s%s"% (X,label))))
        setattr(process, X+"PuWeightA", process.puWeightA.clone(src = cms.InputTag("ww%s%s"% (X,label))))
        setattr(process, X+"PuWeightB", process.puWeightB.clone(src = cms.InputTag("ww%s%s"% (X,label))))
        tree.variables.puW  = cms.InputTag(X+"PuWeight")
        tree.variables.puAW = cms.InputTag(X+"PuWeightA")
        tree.variables.puBW = cms.InputTag(X+"PuWeightB")
        seq += getattr(process, X+"PuWeight")
        seq += getattr(process, X+"PuWeightA")
        seq += getattr(process, X+"PuWeightB")
        if puStudy: addExtraPUWeights(process,tree,X+label,seq)
        if mhiggs > 0:
            setattr(process, X+"PtWeight", process.ptWeight.clone(src = cms.InputTag("ww%s%s"% (X,label))))
            tree.variables.kfW = cms.InputTag(X+"PtWeight")
            seq += process.higgsPt
            seq += getattr(process, X+"PtWeight")
    setattr(process,X+"Tree", tree)
    seq += tree
    if args[3] == 'True' or args[3] == 'true': # path already set up
        p = getattr(process,'sel'+X+label)
        p += seq
        setattr(process,'sel'+X+label,p)
    else: # path not already set up
        setattr(process,'sel'+X+label, cms.Path(seq))

process.TFileService = cms.Service("TFileService",fileName = cms.string("tree.root"))


if IsoStudy:
  for X in "elel", "mumu", "elmu", "muel":
    getattr(process,"ww%s%s"% (X,label)).elTag = "wwEleIDMerge"
    getattr(process,"ww%s%s"% (X,label)).muTag = "wwMuonsMergeID"
    getattr(process,"%sTree"% X).cut = cms.string("!isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMax > 20 && ptMin > 10 && passesIP && nExtraLep(10) == 0")
    prepend = process.isoStudySequence + process.wwEleIDMerge + process.wwMuonsMergeID
    getattr(process,"sel%s%s"% (X,label))._seq = prepend + getattr(process,"sel%s%s"% (X,label))._seq



