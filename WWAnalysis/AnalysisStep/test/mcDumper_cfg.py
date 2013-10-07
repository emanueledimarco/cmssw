import FWCore.ParameterSet.Config as cms
import WWAnalysis.Misc.VarParsing as opts
options = opts.VarParsing('analysis')


#-------------------------------------------------------------------------------
# some basic cleanup
del options._register['filePrepend']
del options._register['totalSections']
del options._register['section']
del options._register['secondaryOutputFile']
del options._singletons['filePrepend']
del options._singletons['totalSections']
del options._singletons['section']
del options._singletons['secondaryOutputFile']
del options._register['secondaryInputFiles']
del options._lists['secondaryInputFiles']




options.register ('reportEvery',
                  1000,
                  opts.VarParsing.multiplicity.singleton, # singleton or list
                  opts.VarParsing.varType.int,          # string, int, or float
                  'Report every events')

options.register ('globalTag',
                  'START52_V9',
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.string,
                  'GlobalTag')


#-------------------------------------------------------------------------------
# defaults
options.outputFile = 'latinosYieldSkim.root'
options.maxEvents  = -1 #all events
#-------------------------------------------------------------------------------

options.parseArguments()

inputFiles       = options.inputFiles
outputFile       = options.outputFile
reportEvery      = options.reportEvery
globalTag        = options.globalTag + "::All"


# cmsRun mcDumper_cfg.py  inputFiles=file:/data2/amassiro/CMSSWRoot/Summer12/VBF_HToWWTo2LAndTau2Nu_M-600_8TeV-powheg-pythia6_AODSIM/F86C4A2D-9492-E111-847C-1CC1DE0540F0.root globalTag=START53_V7A  outputFile=testMC.root



# CMSSW Regular Stuff
processName = "MCDumper"
process = cms.Process(processName)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')


# Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery
process.MessageLogger.suppressWarning = cms.untracked.vstring('patTrigger')

process.GlobalTag.globaltag = globalTag

process.source = cms.Source('PoolSource',
      fileNames = cms.untracked.vstring( inputFiles )
      )

process.out    = cms.OutputModule("PoolOutputModule",
      outputCommands =  cms.untracked.vstring(),
      fileName = cms.untracked.string( "test.root" )
      )



## my pu dumper


#process.MyModule = cms.EDAnalyzer('MyRandomAnalyzer',
    #numBins = cms.untracked.int32(100),
    #minBin  = cms.untracked.double(0),
    #maxBin  = cms.untracked.double(100),
    #rootPlotFile = cms.untracked.string("myfile.root")
#)

process.load('WWAnalysis.SkimStep.vertexFiltering_cff')
process.goodPrimaryVertices.src = cms.InputTag("offlinePrimaryVertices")

process.latino = cms.EDAnalyzer('MCDumper',
              GenJetTag = cms.InputTag( "ak5GenJets" ),
              GenMetTag = cms.InputTag( "genMetTrue" )
)




#   _____      _              _       _
#  / ____|    | |            | |     | |
# | (___   ___| |__   ___  __| |_   _| | ___
#  \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \
#  ____) | (__| | | |  __/ (_| | |_| | |  __/
# |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|
#



# Defines which modules and sequences to run
process.mypath = cms.Path(process.goodPrimaryVertices*process.latino)


process.out.outputCommands =  cms.untracked.vstring(
    'drop *'
)


process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outputFile),
      closeFileFast = cms.untracked.bool(True)
  )




