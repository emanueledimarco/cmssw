
import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.yieldProducer_cfi")
process.yieldParams = process.FWLiteParams.clone()

from WWAnalysis.AnalysisStep.cutPSets_cfi import *

swapOutPSetValues(process.yieldParams.selectionParams.wwelel,defaultWW)
process.yieldParams.selectionParams.wwmumu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu = process.yieldParams.selectionParams.wwelel.clone()
swapOutPSetValues(process.yieldParams.selectionParams.wwelmu,oppositeFlavor)

from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff         import wwSamples
isData = False
looper = wwSamples

if isData:
    lumis = LumiList(filename = '/nfs/bluearc/group/jsons/Skim_Mu_And_El_byGiovanni/Final.json')
    process.yieldParams.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.yieldParams.lumisToProcess = lumis.getCMSSWString().split(',')


#loop over the looper
for x in looper:
    #make new pset
    setattr(
        process.yieldParams.inputParams,
        looper[x][0],
        cms.PSet(
            files = cms.vstring(),
            scale = cms.double(looper[x][1]),
            printEvents = cms.bool(isData),
            isMC = cms.bool( not isData ),
            isSignal = cms.bool( looper[x][5]
        )
    )
    getattr(process.yieldParams.inputParams,looper[x][0]).files += looper[x][3]



