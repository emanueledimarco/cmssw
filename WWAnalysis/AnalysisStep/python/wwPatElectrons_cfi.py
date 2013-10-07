import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *
from PhysicsTools.PatAlgos.tools.coreTools import *

from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff import *
patTrigger.processName  = 'HLT8E29'
patTriggerEvent.processName = 'HLT8E29'

electronTriggerMatchHLTEle15 = cms.EDProducer("PATTriggerMatcherDRLessByR",
    src     = cms.InputTag( "patElectrons" ),
    matched = cms.InputTag( "patTrigger" ),
    pathNames      = cms.vstring( '*' ),
    maxDeltaR = cms.double( 0.1 ),
    andOr          = cms.bool( False ),
    filterIds      = cms.vint32( 0 ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterLabels   = cms.vstring( 'hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter' ),
    collectionTags = cms.vstring( '*' ),
    maxDPtRel = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False ),
)

cleanPatElectronsTriggerMatch.src = 'patElectrons'
cleanPatElectronsTriggerMatch.matches = [
    cms.InputTag('electronTriggerMatchHLTEle15')
]

wwPatElectronSequence = cms.Sequence(
    patTrigger *
    electronTriggerMatchHLTEle15 *
    cleanPatElectronsTriggerMatch
)




def setPatElectronsData( process ) :
    process.patTrigger.processName  = 'HLT'
    process.patTriggerEvent.processName = 'HLT'
    #removeMCMatching(process, ['Electrons'])


