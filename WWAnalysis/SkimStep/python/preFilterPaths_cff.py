import FWCore.ParameterSet.Config as cms

from WWAnalysis.SkimStep.preYieldFilter_cff import addPreYieldFilter
from WWAnalysis.SkimStep.preFakeFilter_cff import addPreFakeFilter

def addPreFilters(process,isMC,doBoris,isVV=False):
    addPreYieldFilter(process,isMC,doBoris,isVV)
    addPreFakeFilter(process,isMC)


