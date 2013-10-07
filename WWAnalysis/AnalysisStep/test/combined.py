import FWCore.ParameterSet.Config as cms
import os

from latinosYieldSkimZZ4L import process
from hzzTree2_cfg import process as process2

process.subProcess = cms.SubProcess(process2,
   # Optional SelectEvents parameter can go here.
   outputCommands = cms.untracked.vstring('keep *'),  #Optional parameter, defaulting to "keep *"
   SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('patPath'))                                    
)





### customization ###
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

###  uncomment the correct line depending if you want the output of the skim saved or not ###
#process.schedule = cms.Schedule( process.patPath, process.scrap, process.outpath)
process.schedule = cms.Schedule( process.patPath, process.scrap)



