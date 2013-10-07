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
# process.GlobalTag.globaltag = 'START311_V2::All'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:RMMEFN'),
#     fileNames = cms.untracked.vstring('file:../SkimStep/latinosYieldSkim.root'),
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


# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")


##############################################

process.skimMuMuTillMET = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumuCONVLHT"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 && pfMet() > 20"),                                   
)
process.skimMuElTillMET = process.skimMuMuTillMET.clone( src = "wwmuelCONVLHT" )
process.skimElMuTillMET = process.skimMuMuTillMET.clone( src = "wwelmuCONVLHT" )
process.skimElElTillMET = process.skimMuMuTillMET.clone( src = "wwelelCONVLHT" )


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RMMEFN'),
#     fileName = cms.untracked.string('step2Lite.root'),
     dropMetaData = cms.untracked.string('ALL'),
     outputCommands = cms.untracked.vstring(
        'keep *',
        'drop recoSkimEvents_*_*_*',
        'keep recoSkimEvents_ww*CONVLHT_*_*',
        'drop patJets_slimPatJetsTriggerMatchNoPU_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMuTillMET','selMuElTillMET','selElElTillMET','selMuMuTillMET' ))
)


process.selMuMuTillMET = cms.Path(process.skimMuMuTillMET)
process.selMuElTillMET = cms.Path(process.skimMuElTillMET)
process.selElMuTillMET = cms.Path(process.skimElMuTillMET)
process.selElElTillMET = cms.Path(process.skimElElTillMET)



process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    process.selMuMuTillMET,process.selMuElTillMET,process.selElMuTillMET,process.selElElTillMET,
    # end
    process.e
)




