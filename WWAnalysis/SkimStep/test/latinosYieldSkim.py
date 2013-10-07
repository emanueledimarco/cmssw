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
#-------------------------------------------------------------------------------

options.register ( 'summary',
                  True,
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.bool,
                  'Print run summary')

options.register ('reportEvery',
                  1000,
                  opts.VarParsing.multiplicity.singleton, # singleton or list
                  opts.VarParsing.varType.int,          # string, int, or float
                  'Report every events')

options.register ('eventsToProcess',
                  '',
                  opts.VarParsing.multiplicity.list,
                  opts.VarParsing.varType.string,
                  'Events to process')

options.register ('skipEvents',
                  0, # default value
                  opts.VarParsing.multiplicity.singleton, # singleton or list
                  opts.VarParsing.varType.int,          # string, int, or float
                  'Number of events to skip')

options.register ('isMC',
                  False, # default value
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.bool,
                  'Run in MonteCarlo mode')

options.register ('globalTag',
                  'START52_V9',
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.string,
                  'GlobalTag')

options.register ('doPF2PATAlso',
                  False, # default value
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.bool,
                  'Add pf2PAT leptons')

options.register ('doFakeRates',
                 'None', # default value
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Turn on Fake Rates mode (can be None, \'only\', \'also\'')

options.register ('doTauEmbed',
                 False, # default value
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Turn on DY embedding mode (can be \'True\' or \'False\'')

options.register ('correctMetPhi',
                 False, # default value
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Turn on MET phi corrections (can be \'True\' or \'False\'')

options.register ('doSusy',
                  False,                                    # default value
                  opts.VarParsing.multiplicity.singleton,   # singleton or list
                  opts.VarParsing.varType.bool,
                  'Turn on Susy MC dumper (can be \'True\' or \'False\'')

#-------------------------------------------------------------------------------
# defaults
options.outputFile = 'latinosYieldSkim.root'
options.maxEvents  = -1 #all events
#-------------------------------------------------------------------------------

options.parseArguments()

#   _____  __  __ __  __ ______     
#  |  __ \|  \/  |  \/  |  ____|    
#  | |__) | \  / | \  / | |__   ___ 
#  |  _  /| |\/| | |\/| |  __| / __|
#  | | \ \| |  | | |  | | |____\__ \
#  |_|  \_\_|  |_|_|  |_|______|___/
#                                   

isMC             = options.isMC
doPF2PATAlso     = options.doPF2PATAlso
doFakeRates      = options.doFakeRates
doTauEmbed       = options.doTauEmbed
inputFiles       = options.inputFiles
outputFile       = options.outputFile
skipEvents       = options.skipEvents
maxEvents        = options.maxEvents
summary          = options.summary
eventsToProcess  = options.eventsToProcess
reportEvery      = options.reportEvery
globalTag        = options.globalTag + "::All"
isVV             = False
doBorisGenFilter = False
correctMetPhi    = options.correctMetPhi
doSusy           = options.doSusy



labelJetRho = "RECO"
if doSusy :
    labelJetRho = "HLT"

# CMSSW Regular Stuff
process = cms.Process("Yield")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

#Options
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(summary))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (maxEvents) )
if eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange (eventsToProcess)

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery
process.MessageLogger.suppressWarning = cms.untracked.vstring('patTrigger')

process.GlobalTag.globaltag = globalTag
if not isMC:
         process.GlobalTag.toGet = cms.VPSet(
         cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
         tag = cms.string("TrackProbabilityCalibration_2D_2012DataTOT_v1_offline"),
         connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
         cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
         tag = cms.string("TrackProbabilityCalibration_3D_2012DataTOT_v1_offline"),
         connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
         )


process.source = cms.Source('PoolSource',fileNames=cms.untracked.vstring( inputFiles ), skipEvents=cms.untracked.uint32( skipEvents ) )
process.out    = cms.OutputModule("PoolOutputModule", outputCommands =  cms.untracked.vstring(), fileName = cms.untracked.string( outputFile ) )

# Gives us preFakeFilter and preYieldFilter
from WWAnalysis.SkimStep.preFilterPaths_cff import addPreFilters
addPreFilters(process,isMC,doBorisGenFilter,isVV)

# pat sequence
process.load("PhysicsTools.PatAlgos.patSequences_cff")
### this is necessary to avoid the conflict between PAT and RECO configurations
### see: https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1357.html
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.doAreaFastjet = True
process.kt6PFJets.voronoiRfact = 0.9
process.kt6PFJets.Rho_EtaMax   = cms.double( 4.4)
###


from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.helpers import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *

# Trigger filter
from WWAnalysis.SkimStep.triggerInformation_cff import addTriggerPaths
jetTrigMatches = addTriggerPaths(process)

process.preLeptonSequence = cms.Sequence()

process.load('WWAnalysis.SkimStep.vertexFiltering_cff')
process.preLeptonSequence += process.firstVertexIsGood
process.preLeptonSequence += process.goodPrimaryVertices




# Rho calculations
from WWAnalysis.SkimStep.rhoCalculations_cff import addRhoVariables

addRhoVariables(process,process.preLeptonSequence)
# added them above, so can remove them here
process.pfPileUp.PFCandidates = "particleFlow"
# this breaks the pfIso
# process.pfPileUp.checkClosestZVertex = cms.bool(False)
process.pfNoPileUp.bottomCollection = "particleFlow"
process.patDefaultSequence.remove( process.pfPileUp )
process.patDefaultSequence.remove( process.pfNoPileUp )

# fix for MC 512patch1 vs 522 production (begin)
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
# process.prePatSequence += process.PFTau ---> see process.PFTau
# fix for MC 512patch1 vs 522 production (end)


# electron ID stuff
from WWAnalysis.SkimStep.electronIDs_cff import addElectronIDs
eidModules = addElectronIDs(process,process.preLeptonSequence)

# generator stuff
from WWAnalysis.SkimStep.generatorInformation_cff import addGeneratorInfo
if isMC: 
     addGeneratorInfo(process,process.preLeptonSequence)

     if doSusy :
         from WWAnalysis.SkimStep.generatorInformationSUSY_cff import additionalselect
         process.prunedGen.select.extend (additionalselect)

#print process.prunedGen.select

#  ______ _           _                     _____      _   _      
# |  ____| |         | |                   |  __ \    | | | |     
# | |__  | | ___  ___| |_ _ __ ___  _ __   | |__) |_ _| |_| |__   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  |  ___/ _` | __| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | | |  | (_| | |_| | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| |_|   \__,_|\__|_| |_| 
#                                                                 

process.patElectrons.embedPFCandidate = False
process.patElectrons.embedSuperCluster = True
process.patElectrons.embedTrack = True
process.patElectrons.addElectronID = True
process.electronMatch.matched = "prunedGen"

process.patElectrons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("eleSmurfPF03"),
    cms.InputTag("eleSmurfPF04"),
    cms.InputTag("electronPFIsoChHad03"),
    cms.InputTag("electronPFIsoNHad03"),
    cms.InputTag("electronPFIsoPhoton03"),
    cms.InputTag("electronPFIsoChHad04"),
    cms.InputTag("electronPFIsoNHad04"),
    cms.InputTag("electronPFIsoPhoton04"),
    cms.InputTag("convValueMapProd","dist"),
    cms.InputTag("convValueMapProd","dcot"),
    cms.InputTag("convValueMapProd","passVtxConvert"),
    cms.InputTag("rhoEl"),
    cms.InputTag("rhoElNoPU"),
)


process.patElectrons.isolationValues = cms.PSet(
#     pfNeutralHadrons = cms.InputTag("isoValElectronWithNeutralIso"),
#     pfChargedHadrons = cms.InputTag("isoValElectronWithChargedIso"),
#     pfPhotons = cms.InputTag("isoValElectronWithPhotonIso")
)


#Set the Pat Electrons to use the eID
for module in eidModules:
    setattr(process.patElectrons.electronIDSources,module.label(),cms.InputTag(module.label()))

process.load("WWAnalysis.Tools.convValueMapProd_cfi")
process.convValueMapProd.conversionLabel = "allConversions"
process.load("WWAnalysis.Tools.electronPFIsoMapProd_cfi")
process.eleSmurfPF03 = process.electronPFIsoMapProd.clone()
process.eleSmurfPF03.deltaR = 0.3
process.eleSmurfPF04 = process.electronPFIsoMapProd.clone()
process.eleSmurfPF04.deltaR = 0.4
process.load("WWAnalysis.Tools.electronEGammaPFIsoProd_cfi")
process.preElectronSequence = cms.Sequence(process.convValueMapProd + process.eleSmurfPF03 + process.eleSmurfPF04 + process.pfEGammaIsolationSingleType)

#  __  __                     _____      _   _
# |  \/  |                   |  __ \    | | | |
# | \  / |_   _  ___  _ __   | |__) |_ _| |_| |__
# | |\/| | | | |/ _ \| '_ \  |  ___/ _` | __| '_ \
# | |  | | |_| | (_) | | | | | |  | (_| | |_| | | |
# |_|  |_|\__,_|\___/|_| |_| |_|   \__,_|\__|_| |_|
#


process.patMuons.embedPFCandidate = False
process.patMuons.embedTrack = True
process.patMuons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("muSmurfPF"),
    cms.InputTag("rhoMu"),
     cms.InputTag("rhoMuNoPU"),
)
process.patMuons.isolationValues = cms.PSet()
process.muonMatch.matched = "prunedGen"
process.load("WWAnalysis.Tools.muonPFIsoMapProd_cfi")
process.muSmurfPF = process.muonPFIsoMapProd.clone()
process.preMuonSequence = cms.Sequence(process.muSmurfPF)

#  _   _                             _     _   _____  ______ ___  _____     _______ 
# | \ | |                   /\      | |   | | |  __ \|  ____|__ \|  __ \ /\|__   __|
# |  \| | _____      __    /  \   __| | __| | | |__) | |__     ) | |__) /  \  | |   
# | . ` |/ _ \ \ /\ / /   / /\ \ / _` |/ _` | |  ___/|  __|   / /|  ___/ /\ \ | |   
# | |\  | (_) \ V  V /   / ____ \ (_| | (_| | | |    | |     / /_| |  / ____ \| |   
# |_| \_|\___/ \_/\_/   /_/    \_\__,_|\__,_| |_|    |_|    |____|_| /_/    \_\_|   
                                                                                     
def addFastJetCorrection(process,label,seq="patDefaultSequence",thisRho="kt6PFJets"):
    corrFact = getattr(process,"patJetCorrFactors"+label)
    setattr(process,"patJetCorrFactorsFastJet"+label,corrFact.clone())
    getattr(process,"patJetCorrFactorsFastJet"+label).levels[0] = 'L1FastJet'
    getattr(process,"patJetCorrFactorsFastJet"+label).rho = cms.InputTag(thisRho,"rho",labelJetRho)
    getattr(process,"patJetCorrFactorsFastJet"+label).useRho = cms.bool(True)
    getattr(process,seq).replace(
        getattr(process,"patJetCorrFactors"+label),
        getattr(process,"patJetCorrFactors"+label) +
        getattr(process,"patJetCorrFactorsFastJet"+label)
    )
    getattr(process,"patJets"+label).jetCorrFactorsSource = cms.VInputTag(
        cms.InputTag("patJetCorrFactorsFastJet"+label) ,
        cms.InputTag("patJetCorrFactors"+label) 
    )

from PhysicsTools.PatAlgos.tools.pfTools import *

if doPF2PATAlso:
    usePF2PAT(process,runPF2PAT=True, jetAlgo="AK5", runOnMC=isMC, postfix="PFlow") 


    ##################################
    #  MVA eID
    ##################################
    process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 )#+ process.mvaNonTrigV0 )
    #Electron ID
    process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    #process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")

    #pfLeptons Only
    process.pfLeptonsOnly=cms.Sequence(process.pfParticleSelectionSequencePFlow+process.pfPhotonSequencePFlow+process.pfMuonSequencePFlow+process.pfNoMuonPFlow+process.pfElectronSequencePFlow+process.eidMVASequence+process.patElectronsPFlow+process.patMuonsPFlow)
    #to use default sequence instead
    process.pfAllMuonsPFlow.src=cms.InputTag("pfNoPileUp")
    process.pfNoMuonPFlow.bottomCollection = cms.InputTag("pfNoPileUp")
    #add deltaBeta corr.
    process.pfIsolatedElectronsPFlow.doDeltaBetaCorrection = True
    process.pfIsolatedMuonsPFlow.doDeltaBetaCorrection = True

    ####################################
    #  changes for iso and deltaR
    ####################################
    #muons # special recipe 
    process.pfIsolatedMuonsPFlow.isolationValueMapsCharged = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03PFlow' ) )
    process.pfIsolatedMuonsPFlow.isolationValueMapsNeutral = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03PFlow' ), cms.InputTag( 'muPFIsoValueGamma03PFlow' ) )
    process.pfIsolatedMuonsPFlow.deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03PFlow' )
    process.patMuonsPFlow.isolationValues.pfNeutralHadrons = cms.InputTag( 'muPFIsoValueNeutral03PFlow' )
    process.patMuonsPFlow.isolationValues.pfPhotons = cms.InputTag( 'muPFIsoValueGamma03PFlow' )
    process.patMuonsPFlow.isolationValues.pfChargedHadrons = cms.InputTag( 'muPFIsoValueCharged03PFlow' )
    process.patMuonsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03PFlow' )
    process.patMuonsPFlow.isolationValues.pfChargedAll = cms.InputTag("muPFIsoValueChargedAll03PFlow")
    #all pfMuons considered as isolated
    process.pfIsolatedMuonsPFlow.combinedIsolationCut = cms.double(9999.)
    process.pfIsolatedMuonsPFlow.isolationCut = cms.double(9999.)
    #electrons
    process.pfIsolatedElectronsPFlow.isolationValueMapsCharged = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFIdPFlow' ) )
    process.pfIsolatedElectronsPFlow.isolationValueMapsNeutral = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFIdPFlow' ), cms.InputTag( 'elPFIsoValueGamma03PFIdPFlow' ) )
    process.pfIsolatedElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFIdPFlow' )
    process.patElectronsPFlow.isolationValues.pfNeutralHadrons = cms.InputTag( 'elPFIsoValueNeutral03PFIdPFlow' )
    process.patElectronsPFlow.isolationValues.pfPhotons = cms.InputTag( 'elPFIsoValueGamma03PFIdPFlow' )
    process.patElectronsPFlow.isolationValues.pfChargedHadrons = cms.InputTag( 'elPFIsoValueCharged03PFIdPFlow' )
    process.patElectronsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFIdPFlow' )
    process.patElectronsPFlow.isolationValues.pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFlow")
    #all pfElectrons considered as isolated
    process.pfIsolatedElectronsPFlow.combinedIsolationCut = cms.double(9999.)
    process.pfIsolatedElectronsPFlow.isolationCut = cms.double(9999.)

    if not isMC:
        process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
        removeMCMatchingPF2PAT( process, postfix="PFlow" )
        removeMCMatching( process)

else:
    if not isMC:
        removeMCMatching(process)




#      _      _      _____                                      
#     | |    | |    / ____|                                     
#     | | ___| |_  | (___   ___  __ _ _   _  ___ _ __   ___ ___ 
# _   | |/ _ \ __|  \___ \ / _ \/ _` | | | |/ _ \ '_ \ / __/ _ \
#| |__| |  __/ |_   ____) |  __/ (_| | |_| |  __/ | | | (_|  __/
# \____/ \___|\__| |_____/ \___|\__, |\__,_|\___|_| |_|\___\___|
#                                  | |                          
#                                  |_|                          

def addMETCorrections(process, isMC, correctMetPhi ):

    if not isMC:
         process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

    # for type 1/0
    process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
    process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
    process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMETcorrType0'),
        cms.InputTag('pfJetMETcorr', 'type1')
    )

    process.patMETCorrections.replace(process.producePFMETCorrections,process.type0PFMEtCorrection+process.producePFMETCorrections)
#     process.patMETCorrections.insert(0,process.type0PFMEtCorrection)
#     process.makePatMETs.insert(0,process.type0PFMEtCorrection)

    #for met xy shift
    if correctMetPhi:
        process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

        # use for 2012 Data
        if not isMC:
            process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data

        # use for Spring'12 MC
        if isMC:
            process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc

        process.patDefaultSequence.replace(process.type0PFMEtCorrection,process.type0PFMEtCorrection+process.pfMEtSysShiftCorrSequence)

        process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
            cms.InputTag('pfMETcorrType0'),
            cms.InputTag('pfJetMETcorr', 'type1') ,
            cms.InputTag('pfMEtSysShiftCorr'),
        )


if isMC:
    myCorrLabels = cms.vstring('L1Offset', 'L2Relative', 'L3Absolute')
else:
    myCorrLabels = cms.vstring('L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual')

addJetCollection(
    process,
    cms.InputTag("ak5PFJetsNoPU"),
    algoLabel    = "NoPU",
    typeLabel    = "",
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF',myCorrLabels),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)

#all the other jets:
switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF',myCorrLabels),
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True
)



#### MET corrections ####

addMETCorrections(process, isMC, correctMetPhi)


# add TCVHE
#### experimental configuration from Andrea Rizzi
process.trackCounting3D1st = cms.ESProducer("TrackCountingESProducer",
    impactParameterType = cms.int32(0), ## 0 = 3D, 1 = 2D

    maximumDistanceToJetAxis = cms.double(0.07),
    deltaR = cms.double(-1.0), ## use cut from JTA

    maximumDecayLength = cms.double(5.0),
    nthTrack = cms.int32(1),
    trackQualityClass = cms.string("any")
)
process.trackCountingVeryHighEffBJetTagsAOD = process.trackCountingHighEffBJetTagsAOD.clone( jetTagComputer = 'trackCounting3D1st' )
process.patDefaultSequence.replace(
    process.trackCountingHighEffBJetTagsAOD,
    process.trackCountingHighEffBJetTagsAOD +
    process.trackCountingVeryHighEffBJetTagsAOD
)
process.patJets.discriminatorSources.append(cms.InputTag("trackCountingVeryHighEffBJetTagsAOD"))


# Need to comment out for MVAMet
# Only keep 'em above 7 GeV as the f'in Smurfs
process.selectedPatJets.cut = "correctedJet('Uncorrected').pt > 7"


# add TCVHE
process.trackCountingVeryHighEffBJetTagsNoPU = process.trackCountingHighEffBJetTagsNoPU.clone( jetTagComputer = 'trackCounting3D1st' )
process.patDefaultSequence.replace(
    process.trackCountingHighEffBJetTagsNoPU,
    process.trackCountingHighEffBJetTagsNoPU +
    process.trackCountingVeryHighEffBJetTagsNoPU
)
process.patJetsNoPU.discriminatorSources.append(cms.InputTag("trackCountingVeryHighEffBJetTagsNoPU"))

# Some stuff to save space
process.patJets.embedCaloTowers = False
process.patJetsNoPU.embedCaloTowers = False
process.patJets.addTagInfos = False
process.patJetsNoPU.addTagInfos = False
process.patJets.embedPFCandidates = False
process.patJetsNoPU.embedPFCandidates = False
process.patJets.addAssociatedTracks = False
process.patJetsNoPU.addAssociatedTracks = False

# Not set up correctly by PAT:
process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("selectedPatJetsNoPU") )
process.patDefaultSequence.replace(
    process.cleanPatJets,
    process.cleanPatJets +
    process.cleanPatJetsNoPU 
)

for X in jetTrigMatches:
    oldmatch = getattr(process,X)
    newmatch = oldmatch.clone( src = cms.InputTag("cleanPatJetsNoPU") )
    setattr(process, X+'NoPU', newmatch)
    process.patDefaultSequence.replace(oldmatch, oldmatch+newmatch)

process.cleanPatJetsTriggerMatchNoPU = process.cleanPatJetsTriggerMatch.clone( matches = [X+"NoPU" for X in jetTrigMatches], src = cms.InputTag("cleanPatJetsNoPU") )
process.patDefaultSequence.replace(
    process.cleanPatJetsTriggerMatch,
    process.cleanPatJetsTriggerMatch +
    process.cleanPatJetsTriggerMatchNoPU
)


# Phil Jet ID:
from CMGTools.External.puJetIDAlgo_cff import PhilV1

process.JetIDcleanPatJetsTriggerMatch = cms.EDProducer('PileupJetIdProducer',
             produceJetIds = cms.bool(True),
             jetids = cms.InputTag(""),
                         runMvas = cms.bool(True),
                         jets = cms.InputTag("cleanPatJetsTriggerMatch"),
                         vertexes = cms.InputTag("goodPrimaryVertices"),
                         algos = cms.VPSet(PhilV1),
)
process.JetIDcleanPatJetsTriggerMatchNoPU = process.JetIDcleanPatJetsTriggerMatch.clone( jets ="cleanPatJetsTriggerMatchNoPU" )


process.boostedPatJetsTriggerMatch = cms.EDProducer("PatJetBooster",
    jetTag = cms.InputTag("cleanPatJetsTriggerMatch"),
    vertexTag = cms.InputTag("goodPrimaryVertices"),
    storeJetId = cms.untracked.bool(True),  
    jetIdTag   = cms.InputTag("JetIDcleanPatJetsTriggerMatch:philv1Id"),
    jetMvaTag  = cms.InputTag("JetIDcleanPatJetsTriggerMatch:philv1Discriminant")
)
process.boostedPatJetsTriggerMatchNoPU = process.boostedPatJetsTriggerMatch.clone( 
    jetTag    = "cleanPatJetsTriggerMatchNoPU" ,
    jetIdTag  = "JetIDcleanPatJetsTriggerMatchNoPU:philv1Id" ,
    jetMvaTag = "JetIDcleanPatJetsTriggerMatchNoPU:philv1Discriminant"
) 

process.slimPatJetsTriggerMatch = cms.EDProducer("PATJetSlimmer",
    src = cms.InputTag("boostedPatJetsTriggerMatch"),
    clearJetVars = cms.bool(True),
    clearDaughters = cms.bool(True),
    dropSpecific = cms.bool(False),
)
process.slimPatJetsTriggerMatchNoPU = process.slimPatJetsTriggerMatch.clone( src = "boostedPatJetsTriggerMatchNoPU" ) 

process.patDefaultSequence += (
    ( process.JetIDcleanPatJetsTriggerMatch +
      process.JetIDcleanPatJetsTriggerMatchNoPU ) *
    ( process.boostedPatJetsTriggerMatch +
      process.boostedPatJetsTriggerMatchNoPU ) * 
    ( process.slimPatJetsTriggerMatch     +
      process.slimPatJetsTriggerMatchNoPU )
)

# Add the fast jet correction:
addFastJetCorrection(process,"")
addFastJetCorrection(process,"NoPU","patDefaultSequence","kt6PFJets")

#               _               _____      _ _           _   _                 
#    /\        | |             / ____|    | | |         | | (_)                
#   /  \  _   _| |_ _ __ ___  | |     ___ | | | ___  ___| |_ _  ___  _ __  ___ 
#  / /\ \| | | | __| '__/ _ \ | |    / _ \| | |/ _ \/ __| __| |/ _ \| '_ \/ __|
# / ____ \ |_| | |_| | |  __/ | |___| (_) | | |  __/ (__| |_| | (_) | | | \__ \
#/_/    \_\__,_|\__|_|  \___|  \_____\___/|_|_|\___|\___|\__|_|\___/|_| |_|___/
#                                                                              

process.load("WWAnalysis.Tools.vertexSumPtMapProd_cfi")


process.autreSeq = cms.Sequence( 
    process.vertexMapProd 
#     process.mergedSuperClusters
)

process.load("WWAnalysis.SkimStep.hzz4lDetectorIsolation_cff")

# make the crazy sequence
process.hzzIsoSequence = cms.Sequence(process.hzz4lDetectorIsolationSequence)


#  _____ _                              _   __  __ ______ _______ 
# / ____| |                            | | |  \/  |  ____|__   __|
#| |    | |__   __ _ _ __ __ _  ___  __| | | \  / | |__     | |   
#| |    | '_ \ / _` | '__/ _` |/ _ \/ _` | | |\/| |  __|    | |   
#| |____| | | | (_| | | | (_| |  __/ (_| | | |  | | |____   | |   
# \_____|_| |_|\__,_|_|  \__, |\___|\__,_| |_|  |_|______|  |_|   
#                         __/ |                                   
#                        |___/                                    

# Needed for MVAMet
# process.allReducedPFCands = cms.EDProducer("ReducedCandidatesProducer",
#                                            srcCands    = cms.InputTag("particleFlow",""),
#                                            srcVertices = cms.InputTag("goodPrimaryVertices"),
#                                            dz          = cms.double(999999999.),
#                                            ptThresh    = cms.double(0.0)
#                                            )

process.reducedPFCands = cms.EDProducer("ReducedCandidatesProducer",
                                        srcCands = cms.InputTag("particleFlow",""),
                                        srcVertices = cms.InputTag("goodPrimaryVertices"),
                                        dz = cms.double(0.1),
                                        ptThresh = cms.double(0.5)
                                        )

# Needed for pfIso
process.reducedPFCandsPfNoPU = cms.EDProducer("ReducedCandidatesProducer",
                                              srcCands = cms.InputTag("pfNoPileUp",""),
                                              srcVertices = cms.InputTag("goodPrimaryVertices"),
                                              dz = cms.double(999999999.),
                                              ptThresh = cms.double(0.0)
                                              )


process.load("WWAnalysis.Tools.interestingVertexRefProducer_cfi")
process.load("WWAnalysis.Tools.chargedMetProducer_cfi")

process.patMuonsWithTriggerNoSA = cms.EDFilter("PATMuonRefSelector",
    cut = cms.string("type!=8"),
    src = cms.InputTag("boostedMuons"),
    filter = cms.bool(False)
)

process.lepsForMET = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("patMuonsWithTriggerNoSA"), cms.InputTag("boostedElectrons"))
)


process.lowPtLeps = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("lepsForMET"),
    cut = cms.string("pt>8"),
    filter = cms.bool(False),
)

process.interestingVertexRefProducer.leptonTags = [cms.InputTag("lowPtLeps")]

process.chargedMetProducer.collectionTag = "particleFlow"
process.chargedMetProducer.vertexTag = "interestingVertexRefProducer"
process.trackMetProducer = process.chargedMetProducer.clone(minNeutralPt = 99999., maxNeutralEta = 0)

process.chargedMetSeq = cms.Sequence( ( 
        process.patMuonsWithTriggerNoSA *
        process.lepsForMET * 
        process.lowPtLeps *
        process.interestingVertexRefProducer ) * 
    #process.chargedMetProducer +
    #process.trackMetProducer + 
    #process.allReducedPFCands *
    process.reducedPFCands +
    process.reducedPFCandsPfNoPU
)


# _______              
#|__   __|             
#   | | __ _ _   _ ___ 
#   | |/ _` | | | / __|
#   | | (_| | |_| \__ \
#   |_|\__,_|\__,_|___/
#                      

switchToPFTauHPS(
   process,
   pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
   pfTauLabelNew = cms.InputTag('hpsPFTauProducer'),
   postfix=""
)

# process.selectedPatTaus.cut = (
#     "pt > 8"
#    "pt > 8 && " +
#    "tauID('leadingTrackFinding') > 0.2 && tauID('byLooseIsolation') > 0.2"
# )




#  _____                        _                _              
# / ____|                      | |              | |             
#| (___   ___  _ __ ___   ___  | |     ___ _ __ | |_ ___  _ __  
# \___ \ / _ \| '_ ` _ \ / _ \ | |    / _ \ '_ \| __/ _ \| '_ \ 
# ____) | (_) | | | | | |  __/ | |___|  __/ |_) | || (_) | | | |
#|_____/ \___/|_| |_| |_|\___| |______\___| .__/ \__\___/|_| |_|
#                                         | |                   
#                                         |_|                   
# ____                  _   _             
#|  _ \                | | (_)            
#| |_) | ___   ___  ___| |_ _ _ __   __ _ 
#|  _ < / _ \ / _ \/ __| __| | '_ \ / _` |
#| |_) | (_) | (_) \__ \ |_| | | | | (_| |
#|____/ \___/ \___/|___/\__|_|_| |_|\__, |
#                                    __/ |
#                                   |___/ 

# First boost to get the IP values
# Then boost to add the PF isolation and the 

# needed rho for electron BDTISO
process.rhoElFullEta    = process.rhoEl.clone(rhoTag = cms.untracked.InputTag("kt6PFJets","rho",labelJetRho))
process.patElectrons.userData.userFloats.src  += [ cms.InputTag("rhoElFullEta") ]
process.preLeptonSequence.replace(process.rhoEl, process.rhoEl + process.rhoElFullEta)
process.rhoMuFullEta    = process.rhoMu.clone(rhoTag = cms.untracked.InputTag("kt6PFJets","rho",labelJetRho))
process.patElectrons.userData.userFloats.src  += [ cms.InputTag("rhoMuFullEta") ]
process.preLeptonSequence.replace(process.rhoMu, process.rhoMu + process.rhoMuFullEta)

# add track IP information?
process.load("WWAnalysis.AnalysisStep.leptonBoosting_cff")
if doTauEmbed == False:
    process.preBoostedElectrons = process.boostedElectrons.clone( electronTag = cms.InputTag("cleanPatElectronsTriggerMatch") )
else :
    process.preBoostedElectrons = process.fakePreBoostedElectrons.clone (src = cms.InputTag("cleanPatElectronsTriggerMatch"))

# to correct for btag
if doTauEmbed == True :
    process.jetTracksAssociatorAtVertex.tracks = cms.InputTag("tmfTracks")


# process.preBoostedElectrons = process.boostedElectrons.clone( electronTag = cms.InputTag("cleanPatElectronsTriggerMatch") )
process.preBoostedMuons     = process.boostedMuons.clone( muonTag = cms.InputTag("cleanPatMuonsTriggerMatch") )
process.patDefaultSequence += process.preBoostedElectrons
process.patDefaultSequence += process.preBoostedMuons


# add Iso deposits 
process.load("WWAnalysis.AnalysisStep.isoAdding_cff")
process.boostedElectronsIso = process.isoAddedElectrons.clone( electronTag = "preBoostedElectrons" )
process.boostedMuonsIso = process.isoAddedMuons.clone( muonTag = "preBoostedMuons" )
# from WWAnalysis.SkimStep.hzz4lDetectorIsolation_cff import muIsoFromDepsZZ4L, eleIsoFromDepsZZ4L
# process.boostedMuonsIso.deposits     += muIsoFromDepsZZ4L
# process.boostedElectronsIso.deposits += eleIsoFromDepsZZ4L

# process.patDefaultSequence += process.hzzIsoSequence
process.patDefaultSequence += process.boostedElectronsIso
process.patDefaultSequence += process.boostedMuonsIso

# add MVA Id and MVA Iso
if doTauEmbed == True:
    process.boostedElectronsBDTID = cms.EDProducer("PatElectronUserFloatAdder",
        src = cms.InputTag("boostedElectronsIso"),
        variables = cms.PSet( 
            bdttrig = cms.string("1000"),
            bdtnontrig = cms.string("1000")
        )
    )
    process.boostedElectrons = cms.EDProducer("PatElectronUserFloatAdder",
        src = cms.InputTag("boostedElectronsBDTID"),
        variables = cms.PSet( 
            bdtisonontrig = cms.string("0")
        )
    )
else:
    process.boostedElectronsBDTID = cms.EDProducer("PatElectronBoosterBDTID", src = cms.InputTag("boostedElectronsIso"))
    process.boostedElectrons = cms.EDProducer("PatElectronBoosterBDTIso", src = cms.InputTag("boostedElectronsBDTID"), effectiveAreaTarget = cms.string("Data2011"),  rho = cms.string("rhoElFullEta"))

#process.boostedElectronsBDTID = cms.EDProducer("PatElectronBoosterBDTID", src = cms.InputTag("boostedElectronsIso"))
#process.boostedElectrons = cms.EDProducer("PatElectronBoosterBDTIso", src = cms.InputTag("boostedElectronsBDTID"), effectiveAreaTarget = cms.string("Data2011"),  rho = cms.string("rhoElFullEta"))



process.boostedMuonsBDTID = cms.EDProducer("PatMuonBoosterBDTID", 
                                           src = cms.InputTag("boostedMuonsIso"), 
                                           vertexs = cms.InputTag("goodPrimaryVertices"),
                                           pfCands = cms.InputTag("particleFlow"),
                                           rho = cms.InputTag("kt6PFJets","rho",labelJetRho),
                                           dzCut = cms.double(0.2),
                                           outputName = cms.string("bdtidnontrigDZ"))

process.boostedMuonsBDTIso = cms.EDProducer("PatMuonBoosterBDTIso", 
                                            src = cms.InputTag("boostedMuonsBDTID"),
                                            vertexs = cms.InputTag("goodPrimaryVertices"),
                                            pfCands = cms.InputTag("particleFlow"),
                                            rho = cms.InputTag("kt6PFJets","rho",labelJetRho),
                                            effectiveAreaTarget = cms.string("Fall11MC"),
                                            dzCut = cms.double(0.2),
                                            outputName = cms.string("bdtisonontrigDZ"))

process.boostedMuons = cms.EDProducer("PatMuonBoosterBDTIso", 
                                      src = cms.InputTag("boostedMuonsBDTIso"),
                                      vertexs = cms.InputTag("goodPrimaryVertices"),
                                      pfCands = cms.InputTag("pfNoPileUp"),
                                      rho = cms.InputTag("kt6PFJets","rho",labelJetRho),
                                      effectiveAreaTarget = cms.string("Fall11MC"),
                                      dzCut = cms.double(999999.),
                                      outputName = cms.string("bdtisonontrigPFNOPU"))


process.patDefaultSequence += process.boostedElectronsBDTID
process.patDefaultSequence += process.boostedElectrons
process.patDefaultSequence += process.boostedMuonsBDTID  
process.patDefaultSequence += process.boostedMuonsBDTIso  
process.patDefaultSequence += process.boostedMuons


#   _____      _              _       _      
#  / ____|    | |            | |     | |     
# | (___   ___| |__   ___  __| |_   _| | ___ 
#  \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \
#  ____) | (__| | | |  __/ (_| | |_| | |  __/
# |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|
#                                            

process.out.outputCommands =  cms.untracked.vstring(
    'drop *',
    # Leptons
    'keep *_boostedElectrons_*_*',
    'keep *_boostedMuons_*_*',
    'keep *_cleanPatTausTriggerMatch*_*_*',
    # Jets
    'keep patJets_slimPatJetsTriggerMatch_*_*',
    'keep patJets_slimPatJetsTriggerMatchPFlow_*_*',
    'keep patJets_slimPatJetsTriggerMatchNoPU_*_*',
    #'keep recoGenJets_patJets_genJets_*',
    #'keep recoGenJets_patJetsPFlow_genJets_*',
    #'keep recoGenJets_patJetsNoPU_genJets_*',
    'keep recoGenJets_selectedPatJets_genJets_*',
    'keep recoGenJets_selectedPatJetsPFlow_genJets_*',
    'keep recoGenJets_selectedPatJetsNoPU_genJets_*',
#     'keep patJets_slimPatJetsTriggerMatchCalo_*_*',
#     'keep patJets_slimPatJetsTriggerMatchJPT_*_*',
    # Tracking
    'keep *_goodPrimaryVertices_*_*',
    #'keep *_offlineBeamSpot_*_*',
    # MET
    #'keep *_tcMet_*_*',
    #'keep *_met_*_*',
    'keep *_pfMet_*_*',
    # MC
    'keep *_prunedGen_*_*',
    'keep *_genMetTrue_*_*',
    'keep GenEventInfoProduct_generator_*_*',
    'keep LHEEventProduct_source_*_*',
    'keep *_genParticles_*_*',
    # Trigger
    'keep *_TriggerResults_*_*',
    #'keep *_vertexMapProd_*_*',
    # Misc
    'keep *_addPileupInfo_*_*',
    #'keep *_chargedMetProducer_*_*',
    #'keep *_trackMetProducer_*_*',
# Needed for MVAMet
#    'keep *_allReducedPFCands_*_*',
    'keep *_reducedPFCands_*_*',
#to be checked if replaces above collection
    'keep *_reducedPFCandsPfNoPU_*_*',
#     'keep *_mergedSuperClusters_*_'+process.name_(),
    'keep *_kt6PF*_rho_'+process.name_(),
    'keep double_kt6PFJetsCentralNeutral_rho_RECO',
    # Debug info, usually commented out
    #'keep *_patTrigger_*_*',  
    #'keep *_l1extraParticles_*_*',  
#if doPF2PATAlso...
    'keep *_patMuonsPFlow_*_Yield',
    'keep *_patElectronsPFlow_*_Yield',
    'keep *_pfType1CorrectedMet_*_Yield',
 # met xy shift correction
    'keep *_pfMEtSysShiftCorr_*_Yield',
)

process.prePatSequence  = cms.Sequence( process.preLeptonSequence + process.preElectronSequence + process.preMuonSequence + process.PFTau)
process.postPatSequence = cms.Sequence( process.autreSeq + process.chargedMetSeq )


# In order to use the good primary vertices everywhere (It would be nicer to set the proper inputTags in the first place)
massSearchReplaceAnyInputTag(process.preYieldFilter,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.prePatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.patDefaultSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.postPatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)

#Vertex definition
if doPF2PATAlso:
    massSearchReplaceAnyInputTag(process.patPF2PATSequencePFlow,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"))
    removeMCMatching( process,['Muons','Electrons'], postfix="PFlow" )

process.firstVertexIsGood.vertices = cms.InputTag("offlinePrimaryVertices")
process.goodPrimaryVertices.src = cms.InputTag("offlinePrimaryVertices")

process.scrap      = cms.Path( process.noscraping ) 
process.outpath    = cms.EndPath(process.out)

if  doPF2PATAlso:
    process.patPath = cms.Path( process.preYieldFilter + process.prePatSequence * process.patDefaultSequence * process.pfLeptonsOnly * process.postPatSequence )
    process.fakPath = cms.Path( process.preFakeFilter + process.prePatSequence * process.patDefaultSequence * process.pfLeptonsOnly * process.postPatSequence )
else:
    process.patPath = cms.Path( process.preYieldFilter + process.prePatSequence * process.patDefaultSequence * process.postPatSequence)
    process.fakPath = cms.Path( process.preFakeFilter + process.prePatSequence * process.patDefaultSequence * process.postPatSequence )

process.out.SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('patPath'))

# from WWAnalysis.SkimStep.skimTools import addIsolationInformation
# addIsolationInformation(process)

if doFakeRates == 'None':
    process.schedule = cms.Schedule( process.patPath, process.scrap, process.outpath)
elif doFakeRates == 'also':
    process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'patPath', 'fakPath' ))
    process.countOverlaps = cms.Path(process.preYieldFilter + process.preFakeFilter)
    process.schedule = cms.Schedule( process.patPath, process.fakPath, process.countOverlaps, process.scrap, process.outpath)
elif doFakeRates == 'only':
    process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'fakPath' ))
    process.schedule = cms.Schedule( process.fakPath, process.scrap, process.outpath)
    

if doTauEmbed == True:
    process.out.outputCommands.extend(
   [
       'keep *_generator_*_*' , 
       'keep *_goldenZmumuCandidatesGe2IsoMuons_*_*'
   ]
  )
