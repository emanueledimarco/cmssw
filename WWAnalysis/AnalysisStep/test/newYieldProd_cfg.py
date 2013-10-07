
import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.yieldProducer_cfi")
process.yieldParams = process.FWLiteParams.clone()

from WWAnalysis.AnalysisStep.cutPSets_cfi import *

swapOutPSetValues(process.yieldParams.selectionParams.wwelel0,hReOpt160)
process.yieldParams.selectionParams.wwmumu0 = process.yieldParams.selectionParams.wwelel0.clone()
process.yieldParams.selectionParams.wwelmu0 = process.yieldParams.selectionParams.wwelel0.clone()
process.yieldParams.selectionParams.wwmuel0 = process.yieldParams.selectionParams.wwelel0.clone()
swapOutPSetValues(process.yieldParams.selectionParams.wwelmu0,oppositeFlavor)
swapOutPSetValues(process.yieldParams.selectionParams.wwmuel0,oppositeFlavor)

# process.yieldParams.selectionParams.wwmumuID = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelmuID = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelelID = process.yieldParams.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.yieldParams.selectionParams.wwelmuID,oppositeFlavor)

# process.yieldParams.selectionParams.wwmumuISO = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelmuISO = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelelISO = process.yieldParams.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.yieldParams.selectionParams.wwelmuISO,oppositeFlavor)

# process.yieldParams.selectionParams.wwmumuNOCONV = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelmuNOCONV = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelelNOCONV = process.yieldParams.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.yieldParams.selectionParams.wwelmuNOCONV,oppositeFlavor)

# process.yieldParams.selectionParams.wwmumuIP = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelmuIP = process.yieldParams.selectionParams.wwelel0.clone()
# process.yieldParams.selectionParams.wwelelIP = process.yieldParams.selectionParams.wwelel0.clone()
# swapOutPSetValues(process.yieldParams.selectionParams.wwelmuIP,oppositeFlavor)

# delattr(process.yieldParams.selectionParams,"wwelel0")

process.yieldParams.inputParams = cms.PSet(
    RMMEFOLDER = cms.PSet(
        files = cms.vstring( ),
        scale = cms.double(1),
        printEvents = cms.bool(True),
        isMC = cms.bool(True),
        isSignal = cms.bool(True),
        maxEvents = cms.int32(-1),
    )
)
process.yieldParams.fileOutName = cms.string('data/RMMEFOLDER.root')

from glob import glob
process.yieldParams.inputParams.RMMEFOLDER.files = [ 'file:%s'%x for x in glob('/nfs/bluearc/group/trees/hww/WW_414_SKIM_V00/RMMEFOLDER/*.root') ]

# import commands
# myDir = "/castor/cern.ch/user/m/mangano/hww/skims/WW_413_SKIM_V00/id101160.GGtoH160toWWto2L2Nu/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 

