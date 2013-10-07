import FWCore.ParameterSet.Config as cms
import WWAnalysis.Misc.VarParsing as opts
import re
import sys

options = opts.VarParsing('analysis')

#-------------------------------------------------------------------------------
# some basic cleanup
del options._register['filePrepend']
del options._register['totalSections']
del options._register['section']
del options._register['secondaryOutputFile']
del options._singletons['filePrepend']
del options._singletons['totalSections']
del options._singletons['section']
del options._singletons['secondaryOutputFile']
del options._register['secondaryInputFiles']
del options._lists['secondaryInputFiles']
#-------------------------------------------------------------------------------
options.register ( 'summary',
                  True,
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.bool,
                  'Print run summary')

options.register ('eventsToProcess',
				  '',
				  opts.VarParsing.multiplicity.list,
				  opts.VarParsing.varType.string,
				  'Events to process')

options.register ('skipEvents',
                  0,                                        # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.int,              # string, int, or float
                  'Number of events to skip')

options.register ('label',
				  'XXX',
				  opts.VarParsing.multiplicity.singleton,
				  opts.VarParsing.varType.string,
				  'Label')

options.register ('json',
				  'YYY',
				  opts.VarParsing.multiplicity.singleton,
				  opts.VarParsing.varType.string,
				  'Json file for data')

options.register ('id',
                  0,                                        # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.string,           # string, int, or float
                  'Dataset id')

options.register ('scale',
                  0,                                        # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.float,            # string, int, or float
                  'Scale factor')

options.register ('two',
                  True,                                     # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.bool,             # string, int, or float
                  'Make step2?')

options.register ('doTauEmbed',
                  False,                                    # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.bool,
                  'Turn on DY embedding mode (can be \'True\' or \'False\'')

options.register ('selection',
                  'TightTight',
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.string,           # string, int, or float
                  'Selection level [TightTight,LooseLoose]')

options.register ('doSameSign',
                  False,                                    # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.bool,
                  'Turn on Same Sign mode (can be \'True\' or \'False\'')

options.register ('doType01met',
                  False,                                    # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.bool,
                  'Turn on Type01 met correction Sign mode (can be \'True\' or \'False\'')

options.register ('doSusy',
                  False,                                    # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.bool,
                  'Turn on Susy MC dumper (can be \'True\' or \'False\'')

options.register ('doHiggs',
                       False,                                    # default value
                       opts.VarParsing.multiplicity.singleton,   # singleton or list
                       opts.VarParsing.varType.bool,
                       'Turn on Higgs MC mass dumper (can be \'True\' or \'False\'')

options.register ('doLHE',
                       False,                                    # default value
                       opts.VarParsing.multiplicity.singleton,   # singleton or list
                       opts.VarParsing.varType.bool,
                       'Turn on LHE dumper (can be \'True\' or \'False\'')

options.register ('doGen',
                       False,                                    # default value
                       opts.VarParsing.multiplicity.singleton,   # singleton or list
                       opts.VarParsing.varType.bool,
                       'Turn on gen Variables dumper (can be \'True\' or \'False\'')


#-------------------------------------------------------------------------------
# defaults
options.outputFile = 'step3.root'
options.maxEvents  = -1 #all events
#-------------------------------------------------------------------------------

options.parseArguments()

process = cms.Process("STEP3")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source('PoolSource',fileNames=cms.untracked.vstring( options.inputFiles ), skipEvents=cms.untracked.uint32( options.skipEvents ) )

process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(options.summary))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load("WWAnalysis.AnalysisStep.step3_cff")
from WWAnalysis.AnalysisStep.step3_cff import * # get also functions

def addVarFlags(s3, vars = {}, flags = {} ): 
    for (key,val) in vars.iteritems():
        setattr(s3.variables, key, val)
    for (key,val) in flags.iteritems():
        setattr(s3.flags, key, val)

def addMuVars( s3 ):
    vars = {}
    flags = {}

    muVars = dict(
        normalizedChi2               = "? {0}.isGlobalMuon ? {0}.globalTrack.normalizedChi2 : -9999",
        numberOfValidMuonHits        = "? {0}.isGlobalMuon ? {0}.globalTrack.hitPattern.numberOfValidMuonHits : -9999",
        numberOfMatches              = "? {0}.isGlobalMuon ? {0}.numberOfMatches : -9999",
        TMLastStationTight           = "? {0}.isTrackerMuon ? {0}.muonID('TMLastStationTight') : -9999",
        trkKink                      = "{0}.combinedQuality.trkKink",
        trackerLayersWithMeasurement = "{0}.innerTrack.hitPattern.trackerLayersWithMeasurement",
        numberOfValidPixelHits       = "{0}.innerTrack.hitPattern.numberOfValidPixelHits",
        trackRelErr                  = "abs({0}.track.ptError / {0}.pt)",
#         trackRelErr2                 = "abs({0}.track.ptError / {0}.track.pt)",
    )
    muFlags = dict(
        isGlobalMuon                 = "{0}.isGlobalMuon()",
        isTrakerMuon                 = "{0}.isTrackerMuon()",
        isPFMuon                     = "{0}.isPFMuon",
    )

    for i in [0,1]:
        for (name,raw) in muVars.iteritems():
            formula = ('? abs({0}.pdgId) == 13 ? '+raw+' : -9999').format('candByPt('+str(i)+')')
            vars[name+str(i+1)] = cms.string(formula)
        for (name,raw) in muFlags.iteritems():
            formula = ('? abs({0}.pdgId) == 13 ? '+raw+' : -9999').format('candByPt('+str(i)+')')
            flags[name+str(i+1)] = cms.string(formula)

    addVarFlags(s3, vars = vars, flags = flags)


doLHE            = options.doLHE
doGen            = options.doGen
doHiggs          = options.doHiggs
doSusy           = options.doSusy
doTauEmbed       = options.doTauEmbed
SameSign         = options.doSameSign  

id = 0
json    = None
mhiggs  = 0
wztth = False
dy = False
from WWAnalysis.AnalysisStep.fourthScaleFactors_cff import *
fourthGenSF = 1
fermiSF = 1
puStudy = False  ## set to true to add 16, yes 16 different PU possibilities
IsoStudy = False ## Set to True to get isolation variables (and a tree build only after ID+CONV+IP, without isolation)
                 ## Note: works only if running also the step2
Summer11 = False # set to true if you need to run the Summer11 (changes the PU distro)
Fall11   = False # set to true if you need to run the Fall11   (changes the PU distro)
                 # if both false, it means it is a sample Summer12 !

label = options.label

if '2011' in label: label = label[:label.find('2011')]
if '2012' in label: label = label[:label.find('2012')]
if label in  [ 'SingleElectron', 'DoubleElectron', 'SingleMuon', 'DoubleMuon', 'MuEG']:
    dataset = [label]
    id      = options.id
    json    = options.json
    scalef  = 1
    doPDFvar = False
    doGen = false
    doLHE = false


elif doTauEmbed == True:
    dataset = ["AllEmbed"]
    id      = options.id
    json    = options.json
    scalef  = 1
    doPDFvar = False
    doGen = false
    doLHE = false


# if args[0].find('2011') != -1: args[0] = args[0][ : args[0].find('2011') ]
# if args[0].find('2012') != -1: args[0] = args[0][ : args[0].find('2012') ]
# if args[0] in [ 'SingleElectron', 'DoubleElectron', 'SingleMuon', 'DoubleMuon', 'MuEG']:
#     dataset = [args[0]]; id = args[1]
#     json    = args[2]
#     scalef = 1
else:
    dataset = ['MC', label];
    id = options.id;
    scalef  = options.scale
    dowztth = re.match("wzttH*", label)
    m = re.match("ggToH(\\d+)to.*", label)
    n = re.match("vbfToH(\\d+)to.*", label)
    r = re.match("Graviton2PM*", label)
    s = re.match("Higgs0M*", label)
    t = re.match("SMH125*", label)
    doPDFvar = True
    doGen = True
    if m:
        mhiggs = int(m.group(1))
        fourthGenSF = fourthGenScales[int(m.group(1))]
        fermiSF = 0
    elif n:
        mhiggs = -1*int(n.group(1))
        fermiSF = fermiPhobicScales[int(n.group(1))]
    elif 'DY' in label and ('ElEl' in label or 'MuMu' in label):
        dy = True
    elif dowztth:
        wztth = True
    if m or n or dowztth or r or s or t:
        doHiggs = True

process.step3Tree.cut = process.step3Tree.cut.value().replace("DATASET", dataset[0])
process.step3Tree.variables.trigger  = process.step3Tree.variables.trigger.value().replace("DATASET",dataset[0])
idn = re.sub('[^0-9]','',id)
process.step3Tree.variables.dataset = str(idn)
# process.step3Tree.variables.dataset = str(id)

# addMuVars(process.step3Tree)

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
    process.step3Tree.variables.trpu = cms.string("1")
    process.step3Tree.variables.itpu = cms.string("1")
    process.step3Tree.variables.ootpup1 = cms.string("1")
    process.step3Tree.variables.ootpum1 = cms.string("1")
    process.step3Tree.variables.puW = cms.string("1")
    process.step3Tree.variables.puAW = cms.string("1")
    process.step3Tree.variables.puBW = cms.string("1")


# process.schedule = cms.Schedule()
process.load("WWAnalysis.AnalysisStep.hww_reboosting_cff")
if doPDFvar: process.slimPatJetsTriggerMatch.isData=  cms.untracked.bool(False)

process.preSkim = cms.Path(process.reboosting)

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")

if options.selection == 'TightTight':
    label = "Scenario6"; muon = "wwMuScenario6"; ele = "wwEleScenario6"; softmu = "wwMu4VetoScenario6"; preSeq = cms.Sequence();
elif options.selection == 'LooseLoose':
    label = "Scenario7"; muon = "wwMuScenario7"; ele = "wwEleScenario5"; softmu = "wwMu4VetoScenario6"; preSeq = cms.Sequence();
else:
    raise ValueError('selection must be either TightTight or LooseLoose') 

# step 2 (begin)
if options.two: # path already set up
    from WWAnalysis.AnalysisStep.skimEventProducer_cfi import addEventHypothesis
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","HLT")
    if doTauEmbed == True:
        process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","EmbeddedRECO")
        process.skimEventProducer.mcGenWeightTag = cms.InputTag("generator:minVisPtFilter")
    addEventHypothesis(process,label,muon,ele,softmu,preSeq)

for X in "elel", "mumu", "elmu", "muel":
    if (wztth == True) or (doPDFvar == True):
        getattr(process,"ww%s%s"% (X,label)).mcGenEventInfoTag = "generator"
        getattr(process,"ww%s%s"% (X,label)).genParticlesTag = "prunedGen"

    if doSusy == True :
        getattr(process,"ww%s%s"% (X,label)).genParticlesTag = "prunedGen"

    if doHiggs == True :
        getattr(process,"ww%s%s"% (X,label)).genParticlesTag = "prunedGen"

    if doLHE == True :
        getattr(process,"ww%s%s"% (X,label)).mcLHEEventInfoTag = "source"

    if doGen == True :
        getattr(process,"ww%s%s"% (X,label)).genParticlesTag = "prunedGen"
        getattr(process,"ww%s%s"% (X,label)).genMetTag = "genMetTrue"
	getattr(process,"ww%s%s"% (X,label)).genJetTag = cms.InputTag("selectedPatJets","genJets","Yield")

    if id in ["036", "037", "037c0", "037c1", "037c2", "037c3", "037c4", "037c5", "037c6", "037c7", "037c8", "037c9", "042", "043", "045", "046" ]: # DY-Madgraph sample
        getattr(process,"ww%s%s"% (X,label)).genParticlesTag = "prunedGen"

# step 2 (end)

for X in "elel", "mumu", "elmu", "muel", "ellell":
    tree = process.step3Tree.clone(src = cms.InputTag("ww%s%s"% (X,label) ));
    seq = cms.Sequence()
    setattr(process, X+'TreeSequence', seq)
    setattr(process, X+"Nvtx", process.nverticesModule.clone(probes = cms.InputTag("ww%s%s"% (X,label))))
    seq += getattr(process, X+"Nvtx")
    tree.variables.nvtx = cms.InputTag(X+"Nvtx")
    if IsoStudy: addIsoStudyVariables(process,tree)
    if doLHE: addLHEVariables(process,tree)
    if doGen: addGenVariables(process,tree)

    if dataset[0] == 'MC':
        setattr(process, X+"NPU",  process.nPU.clone(src = cms.InputTag("ww%s%s"% (X,label))))
        if Summer11:
            setattr(process, X+"PuWeight",  process.puWeightS4AB.clone(src = cms.InputTag("ww%s%s"% (X,label))))
            setattr(process, X+"PuWeightA", process.puWeightS4A.clone (src = cms.InputTag("ww%s%s"% (X,label))))
            setattr(process, X+"PuWeightB", process.puWeightS4B.clone (src = cms.InputTag("ww%s%s"% (X,label))))
        elif Fall11:
            setattr(process, X+"PuWeight",  process.puWeightS6AB.clone(src = cms.InputTag("ww%s%s"% (X,label))))
            setattr(process, X+"PuWeightA", process.puWeightS6A.clone (src = cms.InputTag("ww%s%s"% (X,label))))
            setattr(process, X+"PuWeightB", process.puWeightS6B.clone (src = cms.InputTag("ww%s%s"% (X,label))))
        else :
            setattr(process, X+"PuWeight",  process.puWeightS7AB.clone(src = cms.InputTag("ww%s%s"% (X,label)), nTrueInt = cms.bool(True)))
            setattr(process, X+"PuWeightA", process.puWeightS7A.clone (src = cms.InputTag("ww%s%s"% (X,label)), nTrueInt = cms.bool(True)))
            setattr(process, X+"PuWeightB", process.puWeightS7B.clone (src = cms.InputTag("ww%s%s"% (X,label)), nTrueInt = cms.bool(True)))
        tree.variables.trpu     = cms.InputTag(X+"NPU:tr")
        tree.variables.itpu     = cms.InputTag(X+"NPU:it")
        tree.variables.ootpum1  = cms.InputTag(X+"NPU:m1")
        tree.variables.ootpup1  = cms.InputTag(X+"NPU:p1")
        tree.variables.puW  = cms.InputTag(X+"PuWeight")
        tree.variables.puAW = cms.InputTag(X+"PuWeightA")
        tree.variables.puBW = cms.InputTag(X+"PuWeightB")
        seq += getattr(process, X+"NPU")
        seq += getattr(process, X+"PuWeight")
        seq += getattr(process, X+"PuWeightA")
        seq += getattr(process, X+"PuWeightB")
        if puStudy: addExtraPUWeights(process,tree,X+label,seq)
        if dy:
            setattr(process, X+"DYWeight", process.dyWeight.clone(src = cms.InputTag("ww%s%s"% (X,label))))
            tree.variables.kfW = cms.InputTag(X+"DYWeight")
            seq += getattr(process, X+"DYWeight")
        elif mhiggs > 0:
            setattr(process, X+"PtWeight", process.ptWeight.clone(src = cms.InputTag("ww%s%s"% (X,label))))
            tree.variables.kfW = cms.InputTag(X+"PtWeight")
            seq += process.higgsPt
            seq += getattr(process, X+"PtWeight")

        if id in ["036", "037", "037c0", "037c1", "037c2", "037c3", "037c4", "037c5", "037c6", "037c7", "037c8", "037c9", "042", "043", "045", "046" ]: # DY-Madgraph sample
            tree.variables.mctruth = cms.string("getFinalStateMC()")

    if doTauEmbed == True:
        tree.variables.mctruth = cms.string("mcGenWeight()")

    if wztth == True:
        tree.variables.mctruth    = cms.string("mcHiggsProd()")
        tree.variables.mcHWWdecay = cms.string("getWWdecayMC()")

    if doSusy == True :
        tree.variables.susyMstop = cms.string("getSusyStopMass()")
        tree.variables.susyMLSP  = cms.string("getSusyLSPMass()")

    if doHiggs == True :
        tree.variables.MHiggs      = cms.string("getHiggsMass()")
        tree.variables.PtHiggs     = cms.string("getHiggsPt()")
        tree.variables.HEPMCweight = cms.string("HEPMCweight()")

    if doPDFvar == True :
        tree.variables.pdfscalePDF  = cms.string("getPDFscalePDF()")
        tree.variables.pdfx1  = cms.string("getPDFx1()")
        tree.variables.pdfx2  = cms.string("getPDFx2()")
        tree.variables.pdfid1  = cms.string("getPDFid1()")
        tree.variables.pdfid2  = cms.string("getPDFid2()")
        tree.variables.pdfx1PDF  = cms.string("getPDFx1PDF()")
        tree.variables.pdfx2PDF  = cms.string("getPDFx2PDF()")

    setattr(process,X+"Tree", tree)
    seq += tree
    if options.two: # path already set up
        p = getattr(process,'sel'+X+label)
        p += seq
        setattr(process,'sel'+X+label,p)
    else: # path not already set up
        setattr(process,'sel'+X+label, cms.Path(seq))

process.TFileService = cms.Service("TFileService",fileName = cms.string(options.outputFile))


if IsoStudy:
  for X in "elel", "mumu", "elmu", "muel", "ellell":
    getattr(process,"ww%s%s"% (X,label)).elTag = "wwEleIDMerge"
    getattr(process,"ww%s%s"% (X,label)).muTag = "wwMuonsMergeID"
    getattr(process,"%sTree"% X).cut = cms.string("!isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMax > 20 && ptMin > 10 && passesIP && nExtraLep(10) == 0")
    prepend = process.isoStudySequence + process.wwEleIDMerge + process.wwMuonsMergeID
    getattr(process,"sel%s%s"% (X,label))._seq = prepend + getattr(process,"sel%s%s"% (X,label))._seq


if SameSign:
  for X in "elel", "mumu", "elmu", "muel", "ellell":
    getattr(process,"%sTree"% X).cut = cms.string("q(0)*q(1) > 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMax > 20 && ptMin > 10")




