import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *

# make sure it's  7/5 for now
EL_PT_MIN=7; EL_PT_MIN_LOOSE=7
MU_PT_MIN=5; MU_PT_MIN_LOOSE=5

MUID_LOOSE = "pt > %f && %s" % (MU_PT_MIN_LOOSE, MUID_LOOSE_NO_PT_CUT)
ELID_LOOSE = "pt > %f && %s" % (EL_PT_MIN_LOOSE, ELID_LOOSE_NO_PT_CUT)

MUID_GOOD = " && ".join(["pt > %f" % MU_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_MVA_ISO]) ## MVA WP
ELID_GOOD = " && ".join(["pt > %f" % EL_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_MVA, SINGLE_MVA_ISO]) ## MVA WP

MUID_LOOSE_CR = " && ".join(["pt > %f" % MU_PT_MIN, MU_PRESELECTION, SINGLE_SIP_CUT])
ELID_LOOSE_CR = " && ".join(["pt > %f" % EL_PT_MIN, EL_PRESELECTION, SINGLE_SIP_CUT])
