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

is41XRelease = RMME41X
#is41XRelease = False

process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:RMMEFN'),
#    fileNames = cms.untracked.vstring('file:/home/mangano/productionFolder/May10ReRecoStep1/CMSSW_4_2_3_patch2/myTest/latinosYieldSkim.root '),
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

if isMC:
    process.genPath = cms.Path(process.onlyHiggsGen)
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","REDIGI311X")
else:
    process.genPath = cms.Path()
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","HLT")


# 0
process.wwelmu0 = process.skimEventProducer.clone(hypoType = "WWELMU", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmuel0 = process.skimEventProducer.clone(hypoType = "WWMUEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwelel0 = process.skimEventProducer.clone(hypoType = "WWELEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmumu0 = process.skimEventProducer.clone(hypoType = "WWMUMU", muTag="wwMuMatch", elTag = "wwEleMatch")

# if not is41XRelease:
#     process.wwelmu0.vtxTag        = cms.InputTag("offlinePrimaryVerticesWithBS")
#     process.wwmuel0.vtxTag        = cms.InputTag("offlinePrimaryVerticesWithBS")
#     process.wwelel0.vtxTag        = cms.InputTag("offlinePrimaryVerticesWithBS")
#     process.wwmumu0.vtxTag        = cms.InputTag("offlinePrimaryVerticesWithBS")


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


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RMMEFN'),
#     fileName = cms.untracked.string('step2Lite.root'),
     dropMetaData = cms.untracked.string('ALL'),
     outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_onlyHiggsGen_*_*',
        'keep GenEventInfoProduct_generator__*',
        'keep patJets_slimPatJetsTriggerMatch__Yield',
        'keep PileupSummaryInfos_addPileupInfo__*',
        'keep recoPFMETs_pfMet__*',
        'keep recoMETs_tcMet__*',
        'keep edmTriggerResults_*_*_*',
        'keep recoSkimEvents_ww*CONVLHT_*_*',
        'keep *_goodPrimaryVertices_*_Yield',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMuTillConv','selMuElTillConv','selElElTillConv','selMuMuTillConv' ))
)

# if is41XRelease:
#     process.out.outputCommands.append('keep *_offlinePrimaryVertices_*_Yield')
# else:
#     process.out.outputCommands.append('keep *_offlinePrimaryVerticesWithBS_*_*')


process.selMuMuTillConv = cms.Path(process.wwElectronSequence + process.wwMuonSequence 
                                  * process.wwmumu0 * process.wwmumuIDLHT * process.wwmumuISOLHT * process.wwmumuCONVLHT *
                                  process.skimMuMuTillConv)

process.selMuElTillConv = cms.Path(process.wwElectronSequence + process.wwMuonSequence 
                                  * process.wwmuel0 * process.wwmuelIDLHT * process.wwmuelISOLHT * process.wwmuelCONVLHT *
                                  process.skimMuElTillConv)

process.selElMuTillConv = cms.Path(process.wwElectronSequence + process.wwMuonSequence 
                                  * process.wwelmu0 * process.wwelmuIDLHT * process.wwelmuISOLHT * process.wwelmuCONVLHT *
                                  process.skimElMuTillConv)

process.selElElTillConv = cms.Path(process.wwElectronSequence + process.wwMuonSequence 
                                  * process.wwelel0 * process.wwelelIDLHT * process.wwelelISOLHT * process.wwelelCONVLHT *
                                  process.skimElElTillConv)



# In order to use the offline vertices with BS constratint everywhere 
# if not is41XRelease:
#     import PhysicsTools.PatAlgos.tools.helpers as configtools
#     configtools.massSearchReplaceAnyInputTag(process.wwElectronSequence,
#                                              cms.InputTag("offlinePrimaryVertices"), 
#                                              cms.InputTag("offlinePrimaryVerticesWithBS"),False)
#     configtools.massSearchReplaceAnyInputTag(process.wwMuonSequence,
#                                              cms.InputTag("offlinePrimaryVertices"), 
#                                              cms.InputTag("offlinePrimaryVerticesWithBS"),False)


process.e = cms.EndPath(process.out)


process.schedule = cms.Schedule(
    process.genPath,
    process.selMuMuTillConv,process.selMuElTillConv,process.selElMuTillConv,process.selElElTillConv,
    # end
    process.e
)


