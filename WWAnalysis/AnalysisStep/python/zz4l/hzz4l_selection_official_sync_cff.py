import FWCore.ParameterSet.Config as cms

## Include MELA cut

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import SEL_ZZ4L_STEP_4,SEL_ZZ4L_STEP_5,SEL_ZZ4L_STEP_6,SKIM_SEQUENCE
SEL_ZZ4L_STEP_4 = "mass > 70"
SEL_ZZ4L_STEP_5 = "mass > 100 && mz(1) > 12"
#SEL_ZZ4L_STEP_6 = "userFloat('melaSMH')/(userFloat('melaSMH')+userFloat('melaQQZZ')) > 0.1"
SEL_ZZ4L_STEP_6 = "userFloat('mela') > 0.1"
SKIM_SEQUENCE = ""
