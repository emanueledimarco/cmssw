import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.zz4l.hzz4l_selection_2012_cff import *

# No ISO in good lepton definition, is applied at Z step
MUID_GOOD = " && ".join(["pt > %f" % MU_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW])
ELID_GOOD = " && ".join(["pt > %f" % EL_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW])

DO_FSR_RECOVERY = True
FSR_MATCH = "goodLepIso"
