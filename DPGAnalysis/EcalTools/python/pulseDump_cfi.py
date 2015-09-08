import FWCore.ParameterSet.Config as cms

pulseDump = cms.EDAnalyzer("PulseDump",
    ChannelStatusToBeExcluded = cms.vstring(   'kNoisy',
                                               'kNNoisy',
                                               'kFixedG6',
                                               'kFixedG1',
                                               'kFixedG0',
                                               'kNonRespondingIsolated',
                                               'kDeadVFE',
                                               'kDeadFE',
                                               'kNoDataNoTP',),
                            
    # reco flags association to DB flag
    flagsMapDBReco = cms.PSet(
        kGood  = cms.vstring('kOk','kDAC','kNoLaser','kNoisy'),
        kNoisy = cms.vstring('kNNoisy','kFixedG6','kFixedG1'),
        kNeighboursRecovered = cms.vstring('kFixedG0',
                                                                                   'kNonRespondingIsolated',
                                                                                   'kDeadVFE'),
        kTowerRecovered = cms.vstring('kDeadFE'),
        kDead           = cms.vstring('kNoDataNoTP')
        ),
    minAmplitudeBarrel = cms.double(42),
    minAmplitudeEndcap = cms.double(78),
    pedestalAnalysis = cms.bool(False)
)
