import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *

boostedElectronsEAPFIso.rho = 'rhoElFullEta'
boostedElectronsEAPFIso.effectiveAreaTarget = 'Data2012'

boostedMuonsEAPFIso.rho = 'rhoMuCentralNeutral'
boostedMuonsEAPFIso.effectiveAreaTarget = 'Data2012'
