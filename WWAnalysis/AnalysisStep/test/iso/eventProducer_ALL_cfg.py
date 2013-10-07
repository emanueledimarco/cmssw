import FWCore.ParameterSet.Config as cms

process = cms.Process("WWEvents")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

process.GlobalTag.globaltag = 'START38_V12::All'
# process.GlobalTag.globaltag = 'SET_GLOBALTAG'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
##############################################

from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/nov12Skim/WWTo2L2NuZ26noPU/*.root') ]
process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/hww/WW_39X_ISO_V01/id101160.Flat/*.root') ]

# process.source.fileNames = ['file:hwwSkim.root']


process.source.inputCommands = cms.untracked.vstring( 
    "keep *", 
    "drop *_conditionsInEdm_*_*",  
    "drop *_MEtoEDMConverter_*_*"
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")


### HERE I SET THE SKIMEVENT PRODUCER CONFIGURATION ###
process.wwElectrons.filter = cms.bool(False)
process.wwMuons.filter = cms.bool(False)
process.wwMuons4Veto.filter = cms.bool(False)

process.wwelmu = process.skimEventProducer.clone()
process.wwelmu.hypoType = cms.string("WWELMU")
process.wwmuel = process.skimEventProducer.clone()
process.wwmuel.hypoType = cms.string("WWMUEL")
process.wwelel = process.skimEventProducer.clone()
process.wwelel.hypoType = cms.string("WWELEL")
process.wwmumu = process.skimEventProducer.clone()
process.wwmumu.hypoType = cms.string("WWMUMU")


##############################################

V01 = False
V02 = False
if V01 or V02:
    process.wwelmu.jetTag = "cleanPatJets"
    process.wwmuel.jetTag = "cleanPatJets"
    process.wwelel.jetTag = "cleanPatJets"
    process.wwmumu.jetTag = "cleanPatJets"
    process.wwelmu.__delattr__("sptTag")
    process.wwmuel.__delattr__("sptTag")
    process.wwelel.__delattr__("sptTag")
    process.wwmumu.__delattr__("sptTag")
    process.wwelmu.__delattr__("spt2Tag")
    process.wwmuel.__delattr__("spt2Tag")
    process.wwelel.__delattr__("spt2Tag")
    process.wwmumu.__delattr__("spt2Tag")
    process.wwelmu.__delattr__("tagJetTag")
    process.wwmuel.__delattr__("tagJetTag")
    process.wwelel.__delattr__("tagJetTag")
    process.wwmumu.__delattr__("tagJetTag")

process.skimElMu = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmu"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 && tcMet > 20"),                                   
)
process.skimMuEl = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmuel"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 && tcMet > 20"),                                   
)
process.skimElEl = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelel"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 && tcMet > 20"),                                   
)
process.skimMuMu = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumu"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 && tcMet > 20"),                                   
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_cleanPatJets_*_*',
        'keep *_cleanPatJetsPF_*_*',
        'keep *_goodPrimaryVertices_*_*',
        'keep *_wwelmu_*_*',
        'keep *_wwelel_*_*',
        'keep *_wwmumu_*_*',
        ),
#     verbose = cms.untracked.bool(False),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMu','selMuEl','selElEl','selMuMu' ))
)


process.boostedElectrons.electronTag = "goodElectrons"
process.boostedMuons.muonTag = "goodMuons"
process.wwElectrons.cut = "pt > 10"
process.wwMuons.cut = "pt > 10"
  
process.p = cms.Path( 
    process.wwElectronSequence *
    process.wwMuonSequence *
    process.wwelmu*process.wwmuel*process.wwelel*process.wwmumu
)

process.selElMu = cms.Path(process.skimElMu)
process.selMuEl = cms.Path(process.skimMuEl)
process.selElEl = cms.Path(process.skimElEl)
process.selMuMu = cms.Path(process.skimMuMu)


process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.p,
                             process.selElMu,process.selMuEl,process.selElEl,process.selMuMu,
                             process.e
)



