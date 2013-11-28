import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *
from WWAnalysis.AnalysisStep.zz4l.fsr_cff import *
from WWAnalysis.Tools.fsrPhotons_cff import fsrPhotonPFIsoChHad03pt02, fsrPhotonPFIsoNHad03, fsrPhotonPFIsoPhoton03

fsrPhotonPFIsoChHad03pt02.leptonLabel = 'boostedFsrPhotons'
fsrPhotonPFIsoChHad03pt02.pfLabel = 'reducedPFNoPUCands'
fsrPhotonPFIsoNHad03.leptonLabel = 'boostedFsrPhotons'
fsrPhotonPFIsoNHad03.pfLabel = 'reducedPFNoPUCands'
fsrPhotonPFIsoPhoton03.leptonLabel = 'boostedFsrPhotons'
fsrPhotonPFIsoPhoton03.pfLabel = 'reducedPFNoPUCands'


boostedFsrPhotons = cms.EDProducer("PatPFParticleUserFloatAdder",
    src = cms.InputTag("boostedFsrPhotons","","SKIM"),
    valueMaps = cms.PSet(
        fsrPhotonPFIsoChHad03pt02 = cms.InputTag("fsrPhotonPFIsoChHad03pt02"),
        fsrPhotonPFIsoNHad03 = cms.InputTag("fsrPhotonPFIsoNHad03"),
        fsrPhotonPFIsoPhoton03 = cms.InputTag("fsrPhotonPFIsoPhoton03"),
    ),
    variables = cms.PSet(
        fsrPhotonPFIsoChHadPU03pt02 = cms.string("9.7/16*userFloat('fsrPhotonPFIsoChHadPU04')"), 
    )
)

fixupFsrPhotons = cms.Sequence(
    ( fsrPhotonPFIsoChHad03pt02 + fsrPhotonPFIsoNHad03 + fsrPhotonPFIsoPhoton03 ) *
    boostedFsrPhotons
)

fsrPhotonSeq.replace(fsrPhotonsNoEleSC, fixupFsrPhotons * fsrPhotonsNoEleSC)
