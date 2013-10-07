import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.hzz4l_selection_cff import *

# make sure it's 7/5
MU_PT_MIN=7
EL_PT_MIN=5
# fix IDs
SINGLE_ID_PRL="userInt('prlID')"
MUID_GOOD = " && ".join(["pt > %f" % MU_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_PRL])
ELID_GOOD = " && ".join(["pt > %f" % EL_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_PRL])
MUID_LOOSE_CR = " && ".join(["pt > %f" % MU_PT_MIN, MU_PRESELECTION, SINGLE_SIP_CUT])
ELID_LOOSE_CR = " && ".join(["pt > %f" % EL_PT_MIN, EL_PRESELECTION, SINGLE_SIP_CUT])

SEL_ANY_Z="daughter(0).pdgId = - daughter(1).pdgId && mass > 4 && "+PAIR_DETISO_2D

#ARBITRATE_EARLY = True
SEL_BEST_Z1="40 < mass < 120"
SEL_ZZ4L_STEP_1="lByPt(0).pt > 20 && lByPt(1).pt > 10 && worsePairCombRelIsoBaseline < 0.35"
