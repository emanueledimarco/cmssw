import FWCore.ParameterSet.Config as cms
from glob import glob

process = cms.Process("PassesWW")

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
process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")
print glob('/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/id101160/*')
process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/id101160/*') ]

commonCuts = ('(q(0)*q(1) < 0 && ' +
              '!isSTA(0) && ' +
              '!isSTA(1) && ' +
              'ptMin() > 20.0 && ' +
              'ptMax() > 20.0 && ' +
              'leptEtaCut(2.4, 2.5) && ' +
              'hasGoodVertex() && ' +
              'abs(d0Reco(0)) < 0.02 && ' +
              'abs(d0Reco(1)) < 0.02 && ' +
              'abs(dZReco(0)) < 1.0 && ' +
              'abs(dZReco(1)) < 1.0 && ' +
              '!((abs(pdgId(0)) == 11 && allIso(0)/pt(0) >= 0.10) || ' +
              '(abs(pdgId(1)) == 11 && allIso(1)/pt(1) >= 0.10) || ' +
              '(abs(pdgId(0)) == 13 && allIso(0)/pt(0) >= 0.15) || ' +
              '(abs(pdgId(1)) == 13 && allIso(1)/pt(1) >= 0.15)) && ' +
              'passesIDV1(0) && passesIDV1(1) && ' +
              'passesConversion(0) && passesConversion(1) && ' +
              'tcMet() > 20. && ' +
              'mll() > 12. && ' +
              'nCentralJets( 25., 5., 1) == 0 && ' +
              'nSoftMu(3.) == 0 && ' +
              'nExtraLep(20.) == 0 && ' +
              'bTaggedJetsUnder( 25., 2.1) == 0 )' )

sameFlavor = ('( projTcMet() > 35. && abs(mll()-91.1876)>15 )' )
oppoFlavor = ('( projTcMet() > 20. )')


process.skimElMu = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmu"),
   filter = cms.bool(True),
   cut = cms.string('nLep >=2 && ' + commonCuts + ' && ' + oppoFlavor),                                   
)
process.skimElEl = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelel"),
   filter = cms.bool(True),
   cut = cms.string('nLep >=2 && ' + commonCuts + ' && ' + sameFlavor),                                   
)
process.skimMuMu = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwmumu"),
   filter = cms.bool(True),
   cut = cms.string('nLep >=2 && ' + commonCuts + ' && ' + sameFlavor),                                   
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimmedToWW.root'),
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
    verbose = cms.untracked.bool(False),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMu','selElEl','selMuMu' ))
)

  
process.setupJEC = cms.EDAnalyzer("SkimEventJECInitializer",
    l2File      = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2Relative.txt"),
    l3File      = cms.string("WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L3Absolute.txt"),
    resFile     = cms.string(""),  #set to nothing if running on MC
)

process.setup   = cms.Path(process.setupJEC)
process.selElMu = cms.Path(process.skimElMu)
process.selElEl = cms.Path(process.skimElEl)
process.selMuMu = cms.Path(process.skimMuMu)

process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule( process.setup,process.selElMu,process.selElEl,process.selMuMu, process.e)


