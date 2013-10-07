import FWCore.ParameterSet.Config as cms

process = cms.Process("HwwAnal")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

process.GlobalTag.globaltag = 'START36_V9::All'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
#process.load("WWAnalysis.AnalysisStep.data.pass1.H160_2W_2lnu_gluonfusion_cff")
process.load("WWAnalysis.AnalysisStep.data.passREPLACE_PASS.REPLACEME_cff")
##############################################


process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(400) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwPatElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")


### HERE I SET THE SKIMEVENT PRODUCER CONFIGURATION ###
process.wwElectrons.filter = cms.bool(REPLACE_FILTER_EL)
process.wwMuons.filter = cms.bool(REPLACE_FILTER_MU)
process.hww = process.skimEventProducer.clone()
process.hww.branchAlias = cms.string("synchro")
process.hww.hypoType = cms.string("REPLACE_HYPOTYPE")
process.hww.muTag       = cms.InputTag("wwMuons")
process.hww.elTag       = cms.InputTag("wwElectrons")
process.hww.extraMuTag  = cms.InputTag("wwMuons4Veto")
process.hww.extraElTag  = cms.InputTag("wwElectrons")
process.hww.jetTag      = cms.InputTag("patJets")
process.hww.pfMetTag    = cms.InputTag("pfMet")
process.hww.tcMetTag    = cms.InputTag("tcMet")
process.hww.vtxTag      = cms.InputTag("goodPrimaryVertices")
##############################################


process.skimEventFilter = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("hww"),
   filter = cms.bool(True),
   cut = cms.string("nLep ==2"),                                   
)




process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent_REPLACEME.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_wwMuons_*_*',
        'keep *_wwMuons4Veto_*_*',  
        'keep *_wwElectrons_*_*', 
        'keep *_patJets_*_*',
        'keep *_hww_*_*',
        ),
    verbose = cms.untracked.bool(False),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'p' ))
)

  
process.p = cms.Path( 
    process.wwElectronSequence *
    process.wwMuonSequence *
    process.hww*
    process.skimEventFilter
)
process.e = cms.EndPath(process.out)
