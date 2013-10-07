import FWCore.ParameterSet.Config as cms

process = cms.Process("EWKTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
#     fileNames = cms.untracked.vstring('file:RMMEIN.root'),
    fileNames = cms.untracked.vstring('file:skimmedToWW.root'),
)
process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.HwwTree = cms.EDAnalyzer("SkimTupleProducer", 
#     src = cms.untracked.InputTag("RMMEHYPO"),
    src = cms.untracked.InputTag("wwelel"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("dim01"), quantity = cms.untracked.string("mll()")),
        cms.PSet( tag = cms.untracked.string("pt1"),   quantity = cms.untracked.string("ptByIso(0)")),
        cms.PSet( tag = cms.untracked.string("pt2"),   quantity = cms.untracked.string("ptByIso(1)")),
        cms.PSet( tag = cms.untracked.string("met"), quantity = cms.untracked.string("tcMet()")),
        cms.PSet( tag = cms.untracked.string("delphil"), quantity = cms.untracked.string("dPhill()/3.13159")),
        cms.PSet( tag = cms.untracked.string("deltarll"), quantity = cms.untracked.string("dRll()")),
        cms.PSet( tag = cms.untracked.string("mtw1"), quantity = cms.untracked.string("mTByPt(0)")),
        cms.PSet( tag = cms.untracked.string("mtw2"), quantity = cms.untracked.string("mTByPt(1)")),
        cms.PSet( tag = cms.untracked.string("deltaphilmet"), quantity = cms.untracked.string("dPhilMet()")),
        cms.PSet( tag = cms.untracked.string("metdeltaphilmet"), quantity = cms.untracked.string("projTcMet()")),
        cms.PSet( tag = cms.untracked.string("ltype"), quantity = cms.untracked.string("mitType()")),
        cms.PSet( tag = cms.untracked.string("dilpt"), quantity = cms.untracked.string("pTll()")),
    ),
#     weight = cms.untracked.double(RMMEWEIGHT)
    weight = cms.untracked.double(1)
)

# process.TFileService = cms.Service("TFileService",fileName = cms.string("RMMEOUT"))
process.TFileService = cms.Service("TFileService",fileName = cms.string("mit160.1.root"))

process.p = cms.Path(process.HwwTree)

