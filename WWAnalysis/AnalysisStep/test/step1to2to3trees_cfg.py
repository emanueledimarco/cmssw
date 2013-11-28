import FWCore.ParameterSet.Config as cms

process = cms.Process("Step2Lite")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

process.GlobalTag.globaltag = 'RMMEGlobalTag'
#process.GlobalTag.globaltag = 'START311_V2::All'
#process.GlobalTag.globaltag = 'GR_R_42_V14::All'
#process.GlobalTag.globaltag = 'GR_R_311_V4::All'

is41XRelease = RMME41X
#is41XRelease = True

process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:RMMEFN'),
#   fileNames = cms.untracked.vstring('file:/data/gpetrucc/7TeV/hwww/SingleMuon2011Av2_R414_S1_V06_S2_V02_S3_V00_35_1_3c4.root'),
#   fileNames = cms.untracked.vstring('file:/data/gpetrucc/7TeV/hwww/WJetsToLNu_TuneZ2_7TeV_R414_S1_V06.root'),
    inputCommands = cms.untracked.vstring( "keep *" )
)
##############################################

#from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/nov12Skim/WWTo2L2NuZ26noPU/*.root') ]
#process.source.fileNames += [ 'file:%s'%x for x in glob('/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.root') ]

# process.source.fileNames = ['file:hwwSkim.root']

# import commands
# myDir = "/castor/cern.ch/user/m/mwlebour/hww/trees/WW_414_SKIM_V00/mwlebour/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WW_414_SKIM_V00_ID026_WJetsToLNuMad_Xfer/2bd6c1b7661d2d45dc2e9f49eb14f680/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")


##############################################
process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")

process.onlyHiggsGen = cms.EDProducer( "GenParticlePruner",
    src = cms.InputTag("prunedGen"),
    select = cms.vstring(
        "drop  *  ",
        "keep pdgId =   {h0}",
    )
)

isMC = RMMEMC
#isMC = True
#isMC = False

if isMC:
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","REDIGI311X")
else:
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","HLT")

# 0
process.wwelmu0 = process.skimEventProducer.clone(hypoType = "WWELMU", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmuel0 = process.skimEventProducer.clone(hypoType = "WWMUEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwelel0 = process.skimEventProducer.clone(hypoType = "WWELEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmumu0 = process.skimEventProducer.clone(hypoType = "WWMUMU", muTag="wwMuMatch", elTag = "wwEleMatch")

if not is41XRelease:
    process.wwelmu0.vtxTag        = cms.InputTag("goodPrimaryVerticesWithBS")
    process.wwmuel0.vtxTag        = cms.InputTag("goodPrimaryVerticesWithBS")
    process.wwelel0.vtxTag        = cms.InputTag("goodPrimaryVerticesWithBS")
    process.wwmumu0.vtxTag        = cms.InputTag("goodPrimaryVerticesWithBS")


# LHT
process.wwelmuIDLHT   = process.wwelmu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwmuelIDLHT   = process.wwmuel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwelelIDLHT   = process.wwelel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwmumuIDLHT   = process.wwmumu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )

process.wwelmuISOLHT  = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwmuelISOLHT  = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwelelISOLHT  = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwmumuISOLHT  = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )

process.wwelmuCONVLHT = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwmuelCONVLHT = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwelelCONVLHT = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwmumuCONVLHT = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )


process.skimMuMuTillConv = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumuCONVLHT"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2"),                                   
)
process.skimMuElTillConv = process.skimMuMuTillConv.clone( src = "wwmuelCONVLHT" )
process.skimElMuTillConv = process.skimMuMuTillConv.clone( src = "wwelmuCONVLHT" )
process.skimElElTillConv = process.skimMuMuTillConv.clone( src = "wwelelCONVLHT" )


#noIso = RMMEISO
noIso = True
if noIso: 
    process.wwEleCONVLHT.src = "wwEleIDLHT"
    for X in "elmu", "muel", "mumu", "elel":
        getattr(process, "ww%sCONVLHT"%X).muTag = "wwMuonsID"

process.start = cms.Sequence(
    ( process.wwElectronSequence + 
      process.wwMuonSequence ) *
    ( process.wwmumuCONVLHT +
      process.wwmuelCONVLHT +
      process.wwelmuCONVLHT +
      process.wwelelCONVLHT ) 
)


# In order to use the offline vertices with BS constratint everywhere 
# if not is41XRelease:
#     import PhysicsTools.PatAlgos.tools.helpers as configtools
#     configtools.massSearchReplaceAnyInputTag(process.wwElectronSequence,
#                                              cms.InputTag("goodPrimaryVertices"), 
#                                              cms.InputTag("goodPrimaryVerticesWithBS"),False)
#     configtools.massSearchReplaceAnyInputTag(process.wwMuonSequence,
#                                              cms.InputTag("goodPrimaryVertices"), 
#                                              cms.InputTag("goodPrimaryVerticesWithBS"),False)
pds = "RMMEPDS"; id = int('RMMEID'.replace("A","").replace("B","")); 
#pds = 'MC'; id = 101160; 
#pds = 'MC'; id = 26; 
#pds = 'SingleMuon'; id = 81; 
scalef  = RMMESF
#scalef  = 0.003621062529384;
#mhiggs  = RMMEMHIGGS 
mhiggs  = 0 

pds = pds.replace("2011Av2","").replace("2011A","")

process.load("WWAnalysis.AnalysisStep.step3_cff")

process.step3Tree.cut = "nLep >=2 && " + process.step3Tree.cut.value().replace("DATASET", pds)
process.step3Tree.variables.dataset = str(id)

if pds == "MC":
    from WWAnalysis.AnalysisStep.pileupReweighting_cfi import reWeightVector
    process.step3Tree.eventWeight = cms.InputTag("mcWeight");
    process.mcWeight.baseWeight = scalef
    if mhiggs > 0:
        process.higgsPt.inputFilename = "HiggsAnalysis/HiggsToWW2Leptons/data/kfactors_Std/kfactors_mh%{mass}d_ren%{mass}d_fac%{mass}d.dat" %mass=mhiggs
        process.mcWeight.ptWeight = cms.InputTag("higgsPt")
        process.weightPath = cms.Path(process.onlyHiggsGen + rocess.higgsPt + process.mcWeight)
    else:
        process.weightPath = cms.Path(process.onlyHiggsGen + process.mcWeight)
else:
    process.step3Tree.addRunLumiInfo = cms.addRunLumiInfo = cms.bool(True)

for X in "elel", "mumu", "elmu", "muel":
    setattr(process, X+"Nvtx", process.nverticesModule.clone(probes = cms.InputTag("ww%sCONVLHT"% X)))
    tree = process.step3Tree.clone(src = cms.InputTag("ww%sCONVLHT"% X));
    tree.variables.nvtx = cms.InputTag(X+"Nvtx")
    setattr(process, X+"Tree", tree)
    setattr(process, X+"Path", cms.Path(process.start + getattr(process,X+"Nvtx") + tree))# + getattr(process,X+"TreeNoPU")))

process.TFileService = cms.Service("TFileService",fileName = cms.string("tree.root"))
