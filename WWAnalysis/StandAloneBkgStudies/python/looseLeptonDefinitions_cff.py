import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_ISO as MUON_ISO_TIGHT_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_IP as MUON_IP_TIGHT_CUT




from WWAnalysis.AnalysisStep.wwMuons_cfi import SMURF_ISO as MU_REL_ISO


MUON_ISO_LOOSE_CUT = ("(" + MU_REL_ISO + "< 0.4)")


MUON_IP_LOOSE_CUT  = ("( abs(userFloat('tip')) < 0.1 && abs(userFloat('dzPV')) < 0.1)")


onlyIDMuons = cms.EDFilter("PATMuonSelector",
   src = cms.InputTag("boostedMuons"),
   filter = cms.bool(False),
   #cut = cms.string(MUON_ID_CUT + "&&" + MUON_ISO_LOOSE_CUT + "&&" + MUON_IP_LOOSE_CUT)
   cut = cms.string(MUON_ID_CUT)
)

IdAndIpMuons = cms.EDFilter("PATMuonSelector",
   src = cms.InputTag("boostedMuons"),
   filter = cms.bool(False),
   cut = cms.string(MUON_ID_CUT + "&&" + MUON_IP_LOOSE_CUT)
)

#tightMuons = cms.EDFilter("PATMuonSelector",
#   src = cms.InputTag("boostedMuons"),
#   filter = cms.bool(False),
#   cut = cms.string(MUON_ID_CUT + "&&" + MUON_ISO_TIGHT_CUT + "&&" + MUON_IP_TIGHT_CUT)
#)



