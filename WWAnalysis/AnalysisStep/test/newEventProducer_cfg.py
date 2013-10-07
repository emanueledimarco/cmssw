import FWCore.ParameterSet.Config as cms

process = cms.Process("WWEvents")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')

# isMC = RMMEMC
isMC = True
# isMC = False

# process.GlobalTag.globaltag = 'RMMEGlobalTag'
# process.GlobalTag.globaltag = 'GR_R_42_V14::All'
process.GlobalTag.globaltag = 'START311_V2::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:1160.root'),
#     fileNames = cms.untracked.vstring('file:RMMEFN'),
#     fileNames = cms.untracked.vstring('file:WZ.1.root'),
#     fileNames = cms.untracked.vstring('file:/nfs/bluearc/group/skims/hww/mergedSelV1/ggToH160toWWto2L2Nu.root'),
#     fileNames = cms.untracked.vstring('file:/nfs/bluearc/group/skims/hww/R42X_S1_V04_S2_V00_S3_V00/SingleMuon2011A.step1.root'),
#     fileNames = cms.untracked.vstring('file:/nfs/bluearc/group/skims/hww/R42X_S1_V04_S2_V00_S3_V00/DYtoMuMu.step1.root'),
#     fileNames = cms.untracked.vstring('file:temp.root'),
#     fileNames = cms.untracked.vstring('file:/nfs/bluearc/group/skims/hww/R42X_S1_V04_S2_V00_S3_V00/GluGlu.FirstFile.step1.root'),
    inputCommands = cms.untracked.vstring( "keep *" )
)

# from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/nov12Skim/WWTo2L2NuZ26noPU/*.root') ]

# import commands
# myDir = "/castor/cern.ch/user/m/mwlebour/hww/trees/WW_414_SKIM_V00/mwlebour/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WW_414_SKIM_V00_ID026_WJetsToLNuMad_Xfer/2bd6c1b7661d2d45dc2e9f49eb14f680/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RMMEFN'),
#     fileName = cms.untracked.string('hypoEvents.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep GenEventInfoProduct_generator__*',
        'keep patJets_slimPatJetsTriggerMatch__Yield',
        'keep edmTriggerResults_*_*_*',
        'keep *_goodPrimaryVertices_*_Yield',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_prunedGen_*_*',
        'keep *_goodPrimaryVertices_*_*',
        #PatLetpons
        'keep patMuons_boostedMuons__*',
        'keep patElectrons_boostedElectrons__*',
        # Misc
        'keep *_addPileupInfo_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( )),
    dropMetaData = cms.untracked.string('ALL'),
)
process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule()


##############################################

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
from WWAnalysis.AnalysisStep.skimEventProducer_cfi import addEventHypothesis

process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","HLT")

# process.out and process.schedule need to be defined already
# addEventHypothesis(process,"0","wwMuMatch","wwEleMatch")
# for x in ['ID','ISO','CONV','IP']:
#     addEventHypothesis(process,"{0}Merge" .format(x),"wwMuonsMerge{0}".format(x),"wwEle{0}Merge" .format(x))
#     addEventHypothesis(process,"{0}Merge2".format(x),"wwMuonsMerge{0}".format(x),"wwEle{0}Merge2".format(x))
addEventHypothesis(process,"IPMerge","wwMuonsMergeIP","wwEleIPMerge")

process.schedule.append(process.e)

