import FWCore.ParameterSet.Config as cms

process = cms.Process("PatrickMET")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.GlobalTag.globaltag = 'START311_V2::All'

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/data/mangano/MC/Spring11/GluGluToHToWWTo2L2Nu_M-160_7TeV_Spring11_AOD.root',
        'file:/data/mangano/MC/Spring11/MinBias.root',
    )
)

process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVertices = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVertices.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)


process.badAssPFCandidates0 = cms.EDProducer('PatrickMETProducer',   ptJetCut = cms.double(0) )
process.badAssPFCandidates5 = cms.EDProducer('PatrickMETProducer',   ptJetCut = cms.double(5) )
process.badAssPFCandidates10 = cms.EDProducer('PatrickMETProducer',   ptJetCut = cms.double(10) )
process.badAssPFCandidates15 = cms.EDProducer('PatrickMETProducer',   ptJetCut = cms.double(15) )
process.badAssPFCandidates20 = cms.EDProducer('PatrickMETProducer',   ptJetCut = cms.double(20) )


process.patrickMET0 = process.pfMet.clone(src = cms.InputTag("badAssPFCandidates0"))
process.patrickMET5 = process.pfMet.clone(src = cms.InputTag("badAssPFCandidates5"))
process.patrickMET10 = process.pfMet.clone(src = cms.InputTag("badAssPFCandidates10"))
process.patrickMET15 = process.pfMet.clone(src = cms.InputTag("badAssPFCandidates15"))
process.patrickMET20 = process.pfMet.clone(src = cms.InputTag("badAssPFCandidates20"))


process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('h160.root'),
    fileName = cms.untracked.string('newMet.root'),
    outputCommands =  cms.untracked.vstring(
        'drop *',
        #'keep *_*_*_PatrickMET',
        'keep recoPFMETs_*__*',
        )
)

  
process.p = cms.Path(#process.offlinePrimaryVertices*
                     process.badAssPFCandidates0*
                     process.badAssPFCandidates5*
                     process.badAssPFCandidates10*
                     process.badAssPFCandidates15*
                     process.badAssPFCandidates20*
                     process.patrickMET0*
                     process.patrickMET5*
                     process.patrickMET10*
                     process.patrickMET15*
                     process.patrickMET20
)

process.e = cms.EndPath(process.out)
