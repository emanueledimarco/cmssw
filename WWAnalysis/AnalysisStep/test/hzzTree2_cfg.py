import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Tree")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")

#### in order to use standard regression module ####
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
############



process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    'file:file:hzz4lSkim.root'
    #'root://pcmssd12//data/mangano/MC/8TeV/hzz/step1/step1_id201_42X_S1_V07.root'
    #'root://pcmssd12//data/mangano/MC/8TeV/hzz/step1/step1_id1125_53X_S1_V10.root'
    #'root://pcmssd12//data/mangano/DATA/DoubleMu_HZZ_53X_S1_V10_step1_id010.root'
    #'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/step1/sync/S1_V11/GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A.root'
    #'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/step1/sync/S1_V11/VBF_HToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A_4228BABE-70FA-E111-941B-001A92971B26.S1_V11.root'
    #'file:/data/gpetrucc/7TeV/hzz/step1/sync/S1_V11/HToZZTo4L_M-120_Fall11S6.00215E21D5C4.root',
    #'file:/data/gpetrucc/7TeV/hzz/step1/sync/S1_V11/HToZZTo4L_M-120_Fall11S6.00215E21D5C4.root',
]

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

###### HERE IS THE PART THAT YOU WANT TO CONFIGURE #######
isMC = True
doUseCaliforniaElectronModule = False
doEleRegression = True
EleRegressionType = 2
doEleCalibration = True
doDummyEcalCalib = False
doMuonScaleCorrection = True
NONBLIND = ""
addLeptonPath = False
addZPath = False
doMITBDT = False
doVBF = True
isMuEG = False # used only for data - to turn on/off 3e trigger for 2e2mu events
doZ1Refit = False
doKDAfterZ1Refit = True # will be enabled only if doZ1Refit is True
E_LHC  = 8 # will be set to 7 automatically on 42X, see below
doSyncPaths = False
SKIM_SEQUENCE = ""
###########################################################

cmsswVer=os.environ["CMSSW_VERSION"]
releaseVer="53X" #default
if "CMSSW_5_3_" in cmsswVer:
    releaseVer="53X"

if "CMSSW_5_2_" in cmsswVer:
    releaseVer="52X"

if "CMSSW_4_2_" in cmsswVer:
    releaseVer="42X"
    E_LHC = 7

if "CMSSW_4_4_" in cmsswVer:
    releaseVer="44X"
    E_LHC = 7

muonCalibString = "Summer12_DR53X"

if releaseVer == "42X" : 
    if isMC:
        process.GlobalTag.globaltag = 'START42_V14B::All'   #for 42X MC
        muonCalibString = "Fall11_START42"
    else:
        process.GlobalTag.globaltag = 'GR_R_42_V25::All'  #for 42X DATA
        muonCalibString = "Data2011_42X"
elif releaseVer == "44X" : 
    if isMC:
        process.GlobalTag.globaltag = 'START44_V13::All'   #to read 42X DATA in 44X
        muonCalibString = "Fall11_START42"
    else:
        process.GlobalTag.globaltag = 'GR_R_44_V15C::All'   #to read 42X DATA in 44X
        muonCalibString = "Data2011_42X"
elif releaseVer == "52X" : 
    if isMC:
        process.GlobalTag.globaltag = 'START52_V5::All'   #for 52X MC
    else:
        process.GlobalTag.globaltag = 'GR_R_52_V7::All'   #for 52X DATA
        muonCalibString = "Data2012_53X"
elif releaseVer == "53X" : 
    if isMC:
        process.GlobalTag.globaltag = 'START53_V23::All'   #for 53X MC  
        muonCalibString = "Summer12_DR53X_smearReReco"
    else:
        #process.GlobalTag.globaltag = 'FT_53_V6_AN3::All'  #for 53X DATA July13 ReReco  
        #process.GlobalTag.globaltag = 'FT_53_V6C_AN3::All' #for Aug06 recover  
        #process.GlobalTag.globaltag = 'GR_P_V41_AN3::All'  #for 2012C prompt-data and >=533 release
        #process.GlobalTag.globaltag = 'FT_53_V10_AN3::All' #for 2012C v1 August 24 ReReco    
        #process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'  #for all 2012 data run-ranges ???
        process.GlobalTag.globaltag = 'FT_53_V21_AN4::All'   #22Jan 2013 rereco
        muonCalibString = "Data2012_53X_ReReco"


### =========== Selection ==============
process.load("WWAnalysis.AnalysisStep.hzz4l_selection_cff")

process.load("WWAnalysis.AnalysisStep.zz4l.reSkim_cff")
process.load("WWAnalysis.AnalysisStep.zz4l.mcSequences_cff")
process.load("WWAnalysis.AnalysisStep.zz4l.recoFinalStateClassifiers_cff")
process.load("WWAnalysis.AnalysisStep.zz4lTree_cfi")

## Overrides for synch exercise (note: leave also the other pieces above uncommented as necessary)
#from WWAnalysis.AnalysisStep.zz4l.hzz4l_selection_official_sync_cff import *  

if (releaseVer == "42X" or releaseVer == "44X") :
    print "Use 2011 Selection"
    from WWAnalysis.AnalysisStep.zz4l.hzz4l_selection_2011_fsr_cff import *  
else:
    print "Use 2012 Selection"
    from WWAnalysis.AnalysisStep.zz4l.hzz4l_selection_2012_fsr_cff import *  


### =========== Trigger ==============
if (releaseVer == "42X" or releaseVer == "44X") :
    TRIGGER_FILTER = 'triggerFilter7TeV_MC' if isMC else 'triggerFilter7TeV_DATA'
    doMITBDT = False # Incompatible with this version of TMVA
    if ((not isMC) and isMuEG):
        process.triggerFilter7TeV_DATA.tripleEl = "none"
else:
    TRIGGER_FILTER = 'triggerFilter8TeV'
    if ((not isMC) and isMuEG):
        TRIGGER_FILTER = 'triggerFilter8TeV_No3e'




### =========== BEGIN COMMON PART ==============

### 0a) Do electron energy regression

if doUseCaliforniaElectronModule:     # TO DO: get rid of this once we have the new skims
    process.load("WWAnalysis.AnalysisStep.regressionEnergyPatElectrons_cfi")  
    process.boostedRegressionElectrons = process.regressionEnergyPatElectrons.clone()
else:
    process.load("EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi")
    process.boostedRegressionElectrons = process.eleRegressionEnergy.clone()

process.boostedRegressionElectrons.energyRegressionType = cms.uint32(EleRegressionType)

if (releaseVer == "42X" or releaseVer == "44X"):
    process.boostedRegressionElectrons.rhoCollection = cms.InputTag("kt6PFJetsForIsoActiveArea","rho")
    if (EleRegressionType == 1):
        process.boostedRegressionElectrons.regressionInputFile = cms.string("EGamma/EGammaAnalysisTools/data/eleEnergyRegWeights_V1.root")
        #process.boostedRegressionElectrons.regressionInputFile = cms.string("EGamma/EGammaAnalysisTools/data/weightFile_NoTrkVar_42X.root")
    if (EleRegressionType == 2):
        process.boostedRegressionElectrons.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2011Weights_WithSubClusters_VMay28.root")
else:
    process.boostedRegressionElectrons.rhoCollection = cms.InputTag("kt6PFJets","rho","RECO")
    if (EleRegressionType == 1):
        process.boostedRegressionElectrons.regressionInputFile = cms.string("EGamma/EGammaAnalysisTools/data/eleEnergyRegWeights_V1.root")
        #process.boostedRegressionElectrons.regressionInputFile = cms.string("EGamma/EGammaAnalysisTools/data/weightFile_NoTrkVar_53X.root")
    if (EleRegressionType == 2):
        process.boostedRegressionElectrons.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root")


process.boostedRegressionElectrons.vertexCollection = cms.InputTag('goodPrimaryVertices')
process.boostedRegressionElectrons.debug = cms.untracked.bool(False)


### 0b) Do electron scale calibration
if not hasattr(process, 'RandomNumberGeneratorService'):
    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")



if doUseCaliforniaElectronModule:     # TO DO: get rid of this once we have the new skims
    process.load("WWAnalysis.AnalysisStep.calibratedPatElectrons_cfi")
    process.electronCalibrationAndCombine = process.calibratedPatElectrons.clone()
    process.electronCalibrationAndCombine.isMC = isMC

    #set energy measurement type
    process.electronCalibrationAndCombine.energyMeasurementType = cms.uint32(0)
    if doEleRegression:
        process.electronCalibrationAndCombine.energyMeasurementType = cms.uint32(EleRegressionType)
        
    #set dummy or real corrections
    process.electronCalibrationAndCombine.debug = cms.bool(doDummyEcalCalib)
    if doDummyEcalCalib: process.electronCalibrationAndCombine.smearingRatio = cms.double(1.0)
    
    if isMC : 
        if (releaseVer == "42X" or releaseVer == "44X") : process.electronCalibrationAndCombine.inputDataset = 'Fall11_ICHEP2012'
        else     : process.electronCalibrationAndCombine.inputDataset = 'Summer12_Moriond2013'
    else    : 
        if (releaseVer == "42X" or releaseVer == "44X") : process.electronCalibrationAndCombine.inputDataset = 'Jan16ReReco'
        else     : process.electronCalibrationAndCombine.inputDataset = 'Moriond2013'

    process.electronCalibrationAndCombine.updateEnergyError = cms.bool(True)
    process.electronCalibrationAndCombine.isAOD = cms.bool(True)

else:
    process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
    process.electronCalibrationAndCombine = process.calibratedPatElectrons.clone()
    process.electronCalibrationAndCombine.isMC = cms.bool(isMC)

    #set energy measurement type
    process.electronCalibrationAndCombine.correctionsType = cms.int32(3) #for default correction
    process.electronCalibrationAndCombine.combinationType = cms.int32(1) #std E-P combination by default
    if doEleRegression :        
        if (EleRegressionType == 1):
            process.electronCalibrationAndCombine.correctionsType = cms.int32(1) #corr for old regression
            process.electronCalibrationAndCombine.combinationType = cms.int32(1) #std E-P combination
        if (EleRegressionType == 2):
            process.electronCalibrationAndCombine.correctionsType = cms.int32(2) #corr for new regression
            process.electronCalibrationAndCombine.combinationType = cms.int32(3) #regression E-P combination

    #set dummy or real corrections
    if doDummyEcalCalib:
        process.electronCalibrationAndCombine.verbose = cms.bool(False)
        process.electronCalibrationAndCombine.synchronization = cms.bool(doDummyEcalCalib)
        process.electronCalibrationAndCombine.lumiRatio = cms.double(1.0)
    else:
        process.electronCalibrationAndCombine.lumiRatio = cms.double(0.607)
        
    if isMC : 
        if (releaseVer == "42X" or releaseVer == "44X") : process.electronCalibrationAndCombine.inputDataset = cms.string("Fall11")
        else     : process.electronCalibrationAndCombine.inputDataset = cms.string("Summer12_LegacyPaper")
    else    : 
        if (releaseVer == "42X" or releaseVer == "44X") : process.electronCalibrationAndCombine.inputDataset = cms.string("Jan16ReReco")
        else     : process.electronCalibrationAndCombine.inputDataset = cms.string("22Jan2013ReReco")


    process.electronCalibrationAndCombine.updateEnergyError = cms.bool(True)
#######


process.RandomNumberGeneratorService.electronCalibrationAndCombine = cms.PSet(
    initialSeed = cms.untracked.uint32(1),
    engineName = cms.untracked.string('TRandom3')
)



process.postreboosting = cms.Sequence(
    process.boostedRegressionElectrons * process.electronCalibrationAndCombine * boostedElectronsEAPFIso * boostedElectronsStep2 
)

if doEleRegression:
    if doUseCaliforniaElectronModule:     # TO DO: get rid of this once we have the new skims
        process.boostedRegressionElectrons.inputPatElectronsTag = "boostedElectronsID"
    else:
        process.boostedRegressionElectrons.inputElectronsTag = "boostedElectronsID"
    process.electronCalibrationAndCombine.inputPatElectronsTag = "boostedRegressionElectrons"
    process.boostedElectronsEAPFIso.src = "electronCalibrationAndCombine"   
    if doUseCaliforniaElectronModule and (not doEleCalibration):    # TO DO: get rid of this once we have the new skims
        process.boostedElectronsEAPFIso.src = "boostedRegressionElectrons"   
        process.postreboosting.remove(process.electronCalibrationAndCombine)
 
if (not doEleRegression) and doEleCalibration:
    process.electronCalibrationAndCombine.inputPatElectronsTag = "boostedElectronsID"
    process.boostedElectronsEAPFIso.src = "electronCalibrationAndCombine"   
    process.postreboosting.remove(process.boostedRegressionElectrons)

if (not doEleRegression) and (not doEleCalibration):
    process.postreboosting.remove(process.electronCalibrationAndCombine)
    process.postreboosting.remove(process.boostedRegressionElectrons)



##  0c) Do muon scale calibration 

if doMuonScaleCorrection:
    #process.scaledMuons = cms.EDProducer("RochesterPATMuonCorrector", src = cms.InputTag("boostedMuons"))
    process.scaledMuons = cms.EDProducer("MuScleFitPATMuonCorrector", 
        src = cms.InputTag("boostedMuons"),
        debug = cms.bool(False),
        identifier = cms.string(muonCalibString),
        applySmearing = cms.bool(isMC),
        fakeSmearing = cms.bool(doDummyEcalCalib)
    )
    process.boostedMuonsEAPFIso.src = "scaledMuons"
    process.reboosting.replace(process.boostedMuonsEAPFIso, process.scaledMuons + process.boostedMuonsEAPFIso)



## 1) DEFINE LOOSE LEPTONS 
process.looseMuNoClean = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuonsStep2"),
    cut = cms.string(MUID_LOOSE),
)

process.looseMu = cms.EDProducer(("PATMuonCleanerBySegments" if releaseVer >= "52X" else "PATMuonCleanerBySegmentsUser"),
    src = cms.InputTag("looseMuNoClean"),
    preselection = cms.string("track.isNonnull"),
    passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
    fractionOfSharedSegments = cms.double(0.499),
)

process.looseElNoClean = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectronsStep2"),
    cut = cms.string(ELID_LOOSE),
)

## 2) LEPTONS CLEANING
process.looseEl = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("looseElNoClean"),
    preselection = cms.string(""),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src = cms.InputTag("boostedMuonsStep2"),
           algorithm = cms.string("byDeltaR"),
           preselection = cms.string("userInt('pfMuId')>0 || isGlobalMuon"),
           deltaR  = cms.double(0.05),
           checkRecoComponents = cms.bool(False),
           pairCut  = cms.string(""),
           requireNoOverlaps  = cms.bool(True),
        )
    ),
    finalCut = cms.string(''),
)

process.looseLep = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("looseEl", "looseMu")
)


## 3) DEFINE GOOD LEPTONS 

process.goodMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("looseMu"),
    cut = cms.string(MUID_GOOD),
)

process.goodEl = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("looseEl"),
    cut = cms.string(ELID_GOOD),
)

process.goodLep = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("goodMu", "goodEl")
)

process.countGoodLep = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodLep"),
    minNumber = cms.uint32(4),
)
process.countGoodLep2 = process.countGoodLep.clone(minNumber = 2)
process.countGoodLep3 = process.countGoodLep.clone(minNumber = 3)
process.countLooseLep = process.countGoodLep.clone(src = "looseLep")
process.countLooseLep3 = process.countLooseLep.clone(minNumber = 3)
process.countLooseLep2 = process.countLooseLep.clone(minNumber = 2)

process.countSequenceLLG = cms.Sequence(
    process.countLooseLep3 + 
    process.countLooseLep  +
    process.countGoodLep3  +
    process.countGoodLep
)
process.countSequenceLGG = cms.Sequence(
    process.countLooseLep3 +
    process.countGoodLep3  +
    process.countLooseLep  +
    process.countGoodLep
)

process.slimPatJetsReCorr = cms.EDProducer("PatJetReCorrector",
    jets = cms.InputTag("slimPatJets"),
    payload = cms.string('AK5PF'),
    rho = cms.InputTag('kt6PFJets', 'rho'),
)

if isMC:
    process.slimPatJetsReCorr.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute')
else:
    process.slimPatJetsReCorr.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')

process.goodJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("slimPatJetsReCorr"),
    preselection = cms.string(JETID_GOOD),
    checkOverlaps = cms.PSet(),
    finalCut = cms.string(""),
)

# no need to do the recorrection for 44X
if (releaseVer == "42X" or releaseVer == "44X") : process.goodJets.src = cms.InputTag("slimPatJets")

## 4) MAKE Z CANDIDATES

process.goodLL = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goodLep@+ goodLep@-"),
    cut = cms.string("")
)
process.zllAnyNoFSR = cms.EDProducer("SkimEvent2LProducer",
    src = cms.InputTag("goodLL"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    doMassRes = cms.bool(False),
    isMC = cms.bool(isMC),
)

process.zllAny = cms.EDProducer("SkimEvent2LFsrCollector",
    zll = cms.InputTag("zllAnyNoFSR"),
    photons = cms.InputTag("fsrPhotonsStep2"),
    photonSelection = cms.string("!hasOverlaps('eleVeto') && hasOverlaps('%s')" % FSR_MATCH),
    photonMatch     = cms.string(FSR_MATCH),
    isolationLabel = cms.string("pfCombIso04EACorr"),
    isolationCut   = cms.double(0.4),
    doFsrRecovery = cms.bool(DO_FSR_RECOVERY),
    sortPhotonsByPt = cms.bool(False), ## Use Mike's sorting (by pt if pt > 4, by deltaR otherwise)
)

process.zll = cms.EDFilter("SkimEvent2LSelector",
    src = cms.InputTag("zllAny"),
    cut = cms.string(SEL_ANY_Z),
)

process.oneZ = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zll"),
    minNumber = cms.uint32(1),
)

process.bestZ = cms.EDFilter("SkimEvent2LSorter",
    src = cms.InputTag("zll"),
    sortBy = cms.string("abs(mass - 91.188)"),
    sortAscending = cms.bool(True), 
    maxNumber = cms.uint32(1),
)

process.selectedZ1 = cms.EDFilter("SkimEvent2LSelector",
    src = cms.InputTag("bestZ"),
    cut = cms.string(SEL_BEST_Z1),
    filter = cms.bool(ARBITRATE_EARLY)
)


### =========== BEGIN SIGNAL PART ==============
# Here starts the ZZ step
process.zz = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedZ1 zll" if ARBITRATE_EARLY else "zll zll"),
    cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&"+
                     "deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02 &&"+
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&"+
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02"),
)
process.oneZZ = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zz"),
    minNumber = cms.uint32(1),
)

process.skimEvent4LNoArb = cms.EDProducer("SkimEvent4LProducer",
    src = cms.InputTag("zz"),
    reducedPFCands = cms.InputTag("reducedPFCands"),
    jets = cms.InputTag("goodJets"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    isMC = cms.bool(isMC),
    mcMatch = cms.InputTag(""),
    doswap = cms.bool(False if ARBITRATE_EARLY else True), # sort the two Z's to ensure that Z1 is closest to the nominal mass (unless fixed already)
    doMELA = cms.bool(True),
    energyForMELA = cms.double(E_LHC),                                      
    doMassRes = cms.bool(True),
    doBDT = cms.bool(doMITBDT),
    weightfile_ScalarVsBkgBDT = cms.string("WWAnalysis/AnalysisStep/data/BDTWeights/ScalarVsBkg/hzz4l_mH125_BDTG.weights.xml"),
    gensTag = cms.InputTag("prunedGen"),
    higgsMassWeightFile = cms.string("WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs700_8TeV_Lineshape+Interference.txt"),
    doExtendedMassRes = cms.bool(True),
    doZ1Refit = cms.bool(doZ1Refit),
    doKDAfterZ1Refit = cms.bool(doZ1Refit and doKDAfterZ1Refit),
)

process.zz4lTreeNoArb = process.zz4lTree.clone(src = cms.InputTag("skimEvent4LNoArb"))


process.selectedZZs1 = cms.EDFilter("SkimEvent4LSelector",
    src = cms.InputTag("best4L" if ARBITRATE_EARLY else "skimEvent4LNoArb"),
    cut = cms.string(SEL_ZZ4L_STEP_1),
    filter = cms.bool(True),
)
process.selectedZZs2 = process.selectedZZs1.clone(
    src = "selectedZZs1",
    cut = SEL_ZZ4L_STEP_2,
)
process.selectedZZs3 = process.selectedZZs1.clone(
    src = "selectedZZs2",
    cut = SEL_ZZ4L_STEP_3,
)
process.selectedZZs4 = process.selectedZZs1.clone(
    src = "selectedZZs3",
    cut = SEL_ZZ4L_STEP_4,
)
process.selectedZZs5 = process.selectedZZs1.clone(
    src = "selectedZZs4",
    cut = SEL_ZZ4L_STEP_5,
)
process.selectedZZs6 = process.selectedZZs1.clone(
    src = "selectedZZs5",
    cut = SEL_ZZ4L_STEP_6,
)


process.best4Lpass1 = cms.EDProducer("SkimEvent4LSorterWithTies",
    src = cms.InputTag("skimEvent4LNoArb" if ARBITRATE_EARLY else "selectedZZs6"),
    sortBy = cms.string(SEL_ZZ4L_ARBITRATION_1),
    sortAscending = cms.bool(False),
    maxNumber = cms.uint32(1),
    keepTies = cms.bool(True),
)
process.best4L  = cms.EDProducer("SkimEvent4LSorterWithTies",
    src = cms.InputTag("best4Lpass1"),
    sortBy = cms.string(SEL_ZZ4L_ARBITRATION_2),
    sortAscending = cms.bool(False),
    maxNumber = cms.uint32(1),
    keepTies = cms.bool(False),
)

process.zz4lTree.src = cms.InputTag("selectedZZs6" if ARBITRATE_EARLY else "best4L")

if ARBITRATE_EARLY:
    process.zzCombinatoric = cms.Sequence(
        process.best4Lpass1 + process.best4L +
        process.selectedZZs1 + process.selectedZZs2 + process.selectedZZs3 + process.selectedZZs4 + process.selectedZZs5 + process.selectedZZs6
    )
else:
    process.zzCombinatoric = cms.Sequence(
        process.selectedZZs1 + process.selectedZZs2 + process.selectedZZs3 + process.selectedZZs4 + process.selectedZZs5 + process.selectedZZs6 +
        process.best4Lpass1 + process.best4L 
    )

### =========== CONTROL SELECTION ==========================

### ========= INCLUSIVE MONITORING OF ALL LEPTONS  =============
process.muonTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("boostedMuonsStep2"),
    sortDescendingBy = cms.string("pt"),
    cut = cms.string(""),
    variables   = cms.PSet(
       pt     = cms.string("pt"),
       eta    = cms.string("eta"),
       phi    = cms.string("phi"),
       pfIsoChHad04      = cms.string("userFloat('pfChHadIso04')"),
       pfIsoNHad04_NoEA  = cms.string("userFloat('pfNHadIso04')"),
       pfIsoPhot04_NoEA  = cms.string("userFloat('pfPhotonIso04')"),
       pfIsoEAtot        = cms.string("userFloat('pfCombIso04EACorrEAtot')"),
       pfIsoComb04EACorr = cms.string("userFloat('pfCombIso04EACorr')"),
       bdtIso = cms.string("userFloat('bdtIso')"),
       sip  = cms.string("userFloat('sip')"),
       rho  = cms.string("userFloat('rhoMu')"),
       rhoAA  = cms.string("userFloat('rhoMuActiveArea')"),
    ),
    flags = cms.PSet(
       looseID = cms.string(MUID_LOOSE_NO_PT_CUT),
       mvaISO  = cms.string("userInt('mvaIso')"),
       newID  = cms.string("userInt('newID')"),
       prlID = cms.string("userInt('prlID')"),
    ),
    addRunLumiInfo = cms.bool(True),
)

process.electronTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("boostedElectronsStep2"),
    sortDescendingBy = cms.string("pt"),
    cut = cms.string(""),
    variables   = cms.PSet(
       sceta     = cms.string("superCluster.eta"),
       scphi     = cms.string("superCluster.phi"),
       r9        = cms.string("userFloat('e3x3')/userFloat('rawEnergy')"),
       etawidth  = cms.string("userFloat('etaWidth')"),
       phiwidth  = cms.string("userFloat('phiWidth')"),
       sietaieta = cms.string("sigmaIetaIeta"),
       hOverE    = cms.string("hcalOverEcal"),
       esOverSC  = cms.string("userFloat('esEnergy')/userFloat('rawEnergy')"),
       gsfp      = cms.string("gsfTrack.p"),
       gsfpmode  = cms.string("gsfTrack.momentumMode.R"),
       genp      = cms.string("? genParticleById(11, 1, 1).isNonnull ? genParticleById(11, 1, 1).p() : -1"),
       p         = cms.string("p"),
       pt     = cms.string("pt"),
       eta    = cms.string("eta"),
       phi    = cms.string("phi"),
       scE    = cms.string("superCluster.energy"),
       scrawE = cms.string("userFloat('rawEnergy')"),
       ecalE  = cms.string("ecalEnergy"),
       classification = cms.string("classification"),
       pfIsoChHad04      = cms.string("userFloat('pfChHadIso04')"),
       pfIsoNHad04_NoEA  = cms.string("userFloat('pfNHadIso04')"),
       pfIsoPhot04_NoEA  = cms.string("userFloat('pfPhotonIso04')"),
       pfIsoEAtot        = cms.string("userFloat('pfCombIso04EACorrEAtot')"),
       pfIsoComb04EACorr = cms.string("userFloat('pfCombIso04EACorr')"),
       bdtID  = cms.string("userFloat('bdtID')"),
       bdtIso = cms.string("userFloat('bdtIso')"),
       sip  = cms.string("userFloat('sip')"),
       rho  = cms.string("userFloat('rhoEl')"),
       rhoAA  = cms.string("userFloat('rhoElActiveArea')"),
       missingHits = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits"),
    ),
    flags = cms.PSet(
       looseID = cms.string(ELID_LOOSE_NO_PT_CUT),
       mvaISO = cms.string("userInt('mvaIso')"),
       mvaID  = cms.string("userInt('mvaID')"),
       newID = cms.string("userInt('newID')"),
       prlID = cms.string("userInt('prlID')"),
       isMatched = cms.string("genParticleById(11, 1, 1).isNonnull"),
       isEB      = cms.string("isEB"),
       isEE      = cms.string("isEE"),
       isGold    = cms.string("userFloat('e3x3')/userFloat('rawEnergy') >= 0.94")
    ),
    addRunLumiInfo = cms.bool(True),
)

process.photonTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("fsrPhotonsStep2"),
    sortDescendingBy = cms.string("pt"),
    cut = cms.string("pt > 2"),
    variables   = cms.PSet(
       pt     = cms.string("pt"),
       eta    = cms.string("eta"),
       phi    = cms.string("phi"),
       #pfIsoChHad04 = cms.string("userFloat('fsrPhotonPFIsoChHad04')"),
       #pfIsoChHad04pt02 = cms.string("userFloat('fsrPhotonPFIsoChHad04pt02')"),
       #pfIsoNHad04 = cms.string("userFloat('fsrPhotonPFIsoNHad04')"),  
       #pfIsoPhot04 = cms.string("userFloat('fsrPhotonPFIsoPhoton04')"),  
       #pfIsoChHadPU04 = cms.string("userFloat('fsrPhotonPFIsoChHadPU04')"),
       #pfIsoChHadPU04pt02 = cms.string("userFloat('fsrPhotonPFIsoChHadPU04pt02')"),
       #pfIsoChHad03 = cms.string("userFloat('fsrPhotonPFIsoChHad03')"),
       pfIsoChHad03pt02 = cms.string("userFloat('fsrPhotonPFIsoChHad03pt02')"),
       pfIsoNHad03 = cms.string("userFloat('fsrPhotonPFIsoNHad03')"),  
       pfIsoPhot03 = cms.string("userFloat('fsrPhotonPFIsoPhoton03')"),  
       #pfIsoChHadPU03 = cms.string("userFloat('fsrPhotonPFIsoChHadPU03')"),
       pfIsoChHadPU03pt02 = cms.string("userFloat('fsrPhotonPFIsoChHadPU03pt02')"),
       rho  = cms.string("userFloat('rhoFsrPhoton')"),
       rhoCN  = cms.string("userFloat('rhoFsrPhotonCentralNeutral')"),
       status = cms.string("status"),
       deltaR = cms.string("? hasOverlaps('goodLepNoIso') ? deltaR(overlaps('goodLepNoIso')[0].eta, overlaps('goodLepNoIso')[0].phi, eta, phi) : 999"),
    ),
    flags = cms.PSet(
        eleMatch = cms.string("hasOverlaps('eleVeto')"),
        lepMatch = cms.string("hasOverlaps('%s')" % FSR_MATCH),
        lepMatchNoIso = cms.string("hasOverlaps('goodLepNoIso')"),
    ),
    addRunLumiInfo = cms.bool(True),
)

process.jetTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("slimPatJetsReCorr"),
    sortDescendingBy = cms.string("pt"),
    cut = cms.string("pt > 10"),
    variables   = cms.PSet(
       pt     = cms.string("pt"),
       eta    = cms.string("eta"),
       phi    = cms.string("phi"),
       puJetIDMask = cms.string('userInt(\'jetId\')'),
       jetIDMVA = cms.string('userFloat(\'jetMva\')'),
    ),
    flags = cms.PSet(
        passID = cms.string(JET_ID),
        passPUID = cms.string(JET_PUID),
    ),
    addRunLumiInfo = cms.bool(True),
)
for X in "neutralMultiplicity chargedMultiplicity neutralEmEnergyFraction neutralHadronEnergyFraction chargedEmEnergyFraction chargedHadronEnergyFraction chargedMultiplicity".split():
    setattr(process.jetTree.variables,X,cms.string(X))


### ========= INCLUSIVE MONITORING OF ALL Z'S =============
process.anyZllTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("bestZ"),
    sortDescendingBy = cms.string("-abs(mass - 91.188)"),
    cut = cms.string(""),
    variables   = cms.PSet(
       intimeSimVertices = cms.string("numsimvertices"),
       numTrueInteractions = cms.string("numTrueInteractions"),
       rhoMuCentral = cms.string("luserFloat(0, 'rhoMuCentralNeutral')"),
       l1pt     = cms.string("daughter(0).pt"),
       l1eta    = cms.string("daughter(0).eta"),
       l1phi    = cms.string("daughter(0).phi"),
       l1pdgId  = cms.string("daughter(0).pdgId"),
       l2pt     = cms.string("daughter(1).pt"),
       l2eta    = cms.string("daughter(1).eta"),
       l2phi    = cms.string("daughter(1).phi"),
       l2pdgId  = cms.string("daughter(1).pdgId"),
       zmass    = cms.string("mass"),
       zmll     = cms.string("mll"), # can differ in case of FSR
       zpt      = cms.string("pt"),
       zeta     = cms.string("eta"),
       zphi     = cms.string("phi"),
       l1tkIso  = cms.string("luserFloat(0,'tkZZ4L')"),
       l1ecalIso= cms.string("luserFloat(0,'ecalZZ4L')"),
       l1hcalIso= cms.string("luserFloat(0,'hcalZZ4L')"),
       l2tkIso  = cms.string("luserFloat(0,'tkZZ4L')"),
       l2ecalIso= cms.string("luserFloat(1,'ecalZZ4L')"),
       l2hcalIso= cms.string("luserFloat(1,'hcalZZ4L')"),
       l1combDetIso = cms.string("combinedIso(0)"),
       l2combDetIso = cms.string("combinedIso(1)"),
       pairCombDetIso = cms.string("combinedPairRelativeIso"),
       l1pfIsoComb04EACorr  = cms.string("luserFloat(0, 'pfCombIso04EACorr')"),
       l2pfIsoComb04EACorr  = cms.string("luserFloat(1, 'pfCombIso04EACorr')"),
       l1pfIsoChHad04  = cms.string("luserFloat(0, 'pfChHadIso04')"),
       l2pfIsoChHad04  = cms.string("luserFloat(1, 'pfChHadIso04')"),
       l1pfIsoNHad04  = cms.string("luserFloat(0, 'pfNHadIso04')"),
       l2pfIsoNHad04  = cms.string("luserFloat(1, 'pfNHadIso04')"),
       l1pfIsoPhoton04  = cms.string("luserFloat(0, 'pfPhotonIso04')"),
       l2pfIsoPhoton04  = cms.string("luserFloat(1, 'pfPhotonIso04')"),
       l1pfIsoEAtot     = cms.string("luserFloat(0, 'pfCombIso04EACorrEAtot')"),
       l2pfIsoEAtot     = cms.string("luserFloat(1, 'pfCombIso04EACorrEAtot')"),
       l1bdtID   = cms.string("luserFloat(0,'bdtID')"),
       l2bdtID   = cms.string("luserFloat(1,'bdtID')"),
       l1bdtIso  = cms.string("luserFloat(0,'bdtIso')"),
       l2bdtIso  = cms.string("luserFloat(1,'bdtIso')"),
       l1sip  = cms.string("luserFloat(0,'sip')"),
       l2sip  = cms.string("luserFloat(1,'sip')"),
       rho      = cms.string("luserFloat(1,'rhoEl') + luserFloat(1,'rhoMu')"),
       rhoAA    = cms.string("luserFloat(1,'rhoElActiveArea') + luserFloat(1,'rhoMuActiveArea')"),
       nvtx     = cms.string("numvertices"),
       massErr  = cms.string("userFloat('massErr')"),
       phopt    = cms.string("? numberOfDaughters > 2 ? daughter(2).pt  : -1"),
       phoeta   = cms.string("? numberOfDaughters > 2 ? daughter(2).eta : -1"),
       phophi   = cms.string("? numberOfDaughters > 2 ? daughter(2).phi : -1"),
    ),
    flags = cms.PSet(
       l1idNew  = cms.string("luserInt(0,'newID')"),
       l2idNew  = cms.string("luserInt(1,'newID')"),
       l1idPRL  = cms.string("luserInt(0,'prlID')"),
       l2idPRL  = cms.string("luserInt(1,'prlID')"),
       l1idMVA  = cms.string("luserInt(0,'mvaID')"),
       l2idMVA  = cms.string("luserInt(1,'mvaID')"),
       l1isoMVA  = cms.string("luserInt(0,'mvaIso')"),
       l2isoMVA  = cms.string("luserInt(1,'mvaIso')"),
       fsr       = cms.string("numberOfDaughters > 2"),
    ),
    addRunLumiInfo = cms.bool(True),
)

process.zeetree = process.anyZllTree.clone(cut = cms.string("abs(daughter(1).pdgId) == 11"))
process.zeetree.variables.massErr  = cms.string("userFloat('massErr')")
process.zeetree.flags.l1ConvR  = cms.string("l(0).gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")
process.zeetree.flags.l2ConvR  = cms.string("l(1).gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")

process.zeetree.variables.l1sceta     = cms.string("l(0).superCluster.eta")
process.zeetree.variables.l1scphi     = cms.string("l(0).superCluster.phi")
process.zeetree.variables.l1r9        = cms.string("luserFloat(0, 'e3x3')/luserFloat(0, 'rawEnergy')")
process.zeetree.variables.l1etawidth  = cms.string("luserFloat(0, 'etaWidth')")
process.zeetree.variables.l1phiwidth  = cms.string("luserFloat(0, 'phiWidth')")
process.zeetree.variables.l1sietaieta = cms.string("l(0).sigmaIetaIeta")
process.zeetree.variables.l1hOverE    = cms.string("l(0).hcalOverEcal")
process.zeetree.variables.l1esOverSC  = cms.string("luserFloat(0, 'esEnergy')/luserFloat(0, 'rawEnergy')")
process.zeetree.variables.l1gsfp      = cms.string("l(0).gsfTrack.p")
process.zeetree.variables.l1gsfpmode  = cms.string("l(0).gsfTrack.momentumMode.R")
process.zeetree.variables.l1p         = cms.string("l(0).p")
process.zeetree.variables.l1classification = cms.string("l(0).classification")
process.zeetree.variables.l2sceta     = cms.string("l(1).superCluster.eta")
process.zeetree.variables.l2scphi     = cms.string("l(1).superCluster.phi")
process.zeetree.variables.l2r9        = cms.string("luserFloat(1,'e3x3')/luserFloat(1,'rawEnergy')")
process.zeetree.variables.l2etawidth  = cms.string("luserFloat(1,'etaWidth')")
process.zeetree.variables.l2phiwidth  = cms.string("luserFloat(1,'phiWidth')")
process.zeetree.variables.l2sietaieta = cms.string("l(1).sigmaIetaIeta")
process.zeetree.variables.l2hOverE    = cms.string("l(1).hcalOverEcal")
process.zeetree.variables.l2esOverSC  = cms.string("luserFloat(1, 'esEnergy')/luserFloat(1, 'rawEnergy')")
process.zeetree.variables.l2gsfp      = cms.string("l(1).gsfTrack.p")
process.zeetree.variables.l2gsfpmode  = cms.string("l(1).gsfTrack.momentumMode.R")
process.zeetree.variables.l2p         = cms.string("l(1).p")
process.zeetree.variables.l2classification = cms.string("l(1).classification")



process.zmmtree = process.anyZllTree.clone(cut = cms.string("abs(daughter(1).pdgId) == 13"))
process.zmmtree.variables.l1tkLayers =  cms.string("l(0).track.hitPattern.trackerLayersWithMeasurement")
process.zmmtree.variables.l2tkLayers =  cms.string("l(1).track.hitPattern.trackerLayersWithMeasurement")
process.zmmtree.flags.l1idGlb  = cms.string("l(0).isGlobalMuon")
process.zmmtree.flags.l2idGlb  = cms.string("l(1).isGlobalMuon")

# ========= FAKE RATE MEASUREMENT PART (Z + 1L) =================
process.looseMuCR = process.goodMu.clone(cut = MUID_LOOSE_CR)
process.looseElCR = process.goodEl.clone(cut = ELID_LOOSE_CR)
process.looseLepCR = process.goodLep.clone(src = [ 'looseMuCR', 'looseElCR' ])
process.looseLepCRwithFSRIso = cms.EDProducer("SingleLeptonFsrCollector",
    leptons = cms.InputTag("looseLepCR"),
    label   = cms.string("pfCombRelIso04EACorr_WithFSR"),
    photons = cms.InputTag("fsrPhotonsCR"),
    photonSelection = process.zllAny.photonSelection,
    photonMatch     = process.zllAny.photonMatch,
    isolationLabel  = process.zllAny.isolationLabel,
    isolationCut    = process.zllAny.isolationCut,
    doFsrRecovery   = process.zllAny.doFsrRecovery,
    sortPhotonsByPt = process.zllAny.sortPhotonsByPt,
) 

process.metVeto = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("pfMet"),
    cut = cms.string("et < 25"),
    filter = cms.bool(True)
)
process.lepMinFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("looseLepCR"),
    minNumber = cms.uint32(3),
)

process.lepMaxFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("looseLepCR"), 
    minNumber = cms.uint32(4),
)

Z_PLUS_LEP_MIJ=("sqrt(pow(daughter(0).daughter({0}).energy+daughter(1).energy, 2) - " +
                "     pow(daughter(0).daughter({0}).px    +daughter(1).px, 2) -" +  
                "     pow(daughter(0).daughter({0}).py    +daughter(1).py, 2) -" +  
                "     pow(daughter(0).daughter({0}).pz    +daughter(1).pz, 2))")
process.zPlusLep = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedZ1 looseLepCRwithFSRIso"),
    cut = cms.string( "deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 && "+
                      "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)>0.02 && "+
                     ("(daughter(0).daughter(0).charge == daughter(1).charge || %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(0))) +
                     ("(daughter(0).daughter(1).charge == daughter(1).charge || %s > 4)    " % ( Z_PLUS_LEP_MIJ.format(1)))),
    checkCharge = cms.bool(False)
)


process.zllmtree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("zPlusLep"),
    cut = cms.string("abs(daughter(1).pdgId) == 13"),
    addRunLumiInfo = cms.bool(True),
    variables   = cms.PSet(
       zmass    = cms.string("daughter(0).mass"),
       met      = cms.string("daughter(0).masterClone.met"),
       l1pt     = cms.string("daughter(0).masterClone.daughter(0).pt"),
       l1eta    = cms.string("daughter(0).masterClone.daughter(0).eta"),
       l1phi    = cms.string("daughter(0).masterClone.daughter(0).phi"),
       l1pdgId  = cms.string("daughter(0).masterClone.daughter(0).pdgId"),
       l2pt     = cms.string("daughter(0).masterClone.daughter(1).pt"),
       l2eta    = cms.string("daughter(0).masterClone.daughter(1).eta"),
       l2phi    = cms.string("daughter(0).masterClone.daughter(1).phi"),
       l2pdgId  = cms.string("daughter(0).masterClone.daughter(1).pdgId"),
       pt       = cms.string("daughter(1).pt"),
       eta      = cms.string("daughter(1).eta"),
       phi      = cms.string("daughter(1).phi"),
       tkIso    = cms.string("daughter(1).masterClone.userFloat('tkZZ4L')"),
       ecalIso  = cms.string("daughter(1).masterClone.userFloat('ecalZZ4L')"),
       hcalIso  = cms.string("daughter(1).masterClone.userFloat('hcalZZ4L')"),
       rho      = cms.string("daughter(1).masterClone.userFloat('rhoMu')"),
       rhoAA    = cms.string("daughter(1).masterClone.userFloat('rhoMuActiveArea')"),
       bdtIdDz = cms.string("daughter(1).masterClone.userFloat('bdtidnontrigDZ')"),
       bdtIso  = cms.string("daughter(1).masterClone.userFloat('bdtIso')"),
       bdtIsoDz = cms.string("daughter(1).masterClone.userFloat('bdtisonontrigDZ')"), # deprecated
       bdtIsoPnp = cms.string("daughter(1).masterClone.userFloat('bdtisonontrigPFNOPU')"), # deprecated
       pfIsoChHad04       = cms.string("daughter(1).masterClone.userFloat('pfChHadIso04')"),
       pfIsoNHad04_NoEA   = cms.string("daughter(1).masterClone.userFloat('pfNHadIso04')"),
       pfIsoPhoton04_NoEA = cms.string("daughter(1).masterClone.userFloat('pfPhotonIso04')"),
       pfCombRelIso04EACorr = cms.string("daughter(1).masterClone.userFloat('pfCombRelIso04EACorr')"),
       pfCombRelIso04EACorr_WithFSR = cms.string("daughter(1).masterClone.userFloat('pfCombRelIso04EACorr_WithFSR')"),
       sip    = cms.string("daughter(1).masterClone.userFloat('sip')"),
       sip3d = cms.string("daughter(1).masterClone.userFloat('ip')/daughter(1).masterClone.userFloat('ipErr')"),
       numvertices = cms.string("daughter(0).masterClone.numvertices"),
       intimeSimVertices = cms.string("daughter(0).masterClone.numsimvertices"),
       numTrueInteractions = cms.string("daughter(0).masterClone.numTrueInteractions"),
       dr13 = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)"),
       dr23 = cms.string("deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)"),
       q13  = cms.string("daughter(0).daughter(0).charge + daughter(1).charge"),
       q23  = cms.string("daughter(0).daughter(1).charge + daughter(1).charge"),
       m13  = cms.string("sqrt(pow(daughter(0).daughter(0).energy+daughter(1).energy, 2) - " +
                         "     pow(daughter(0).daughter(0).px    +daughter(1).px, 2) -" +  
                         "     pow(daughter(0).daughter(0).py    +daughter(1).py, 2) -" +  
                         "     pow(daughter(0).daughter(0).pz    +daughter(1).pz, 2))"),   
       m23  = cms.string("sqrt(pow(daughter(0).daughter(1).energy+daughter(1).energy, 2) - "+
                         "     pow(daughter(0).daughter(1).px    +daughter(1).px, 2) -"+   
                         "     pow(daughter(0).daughter(1).py    +daughter(1).py, 2) -"+   
                         "     pow(daughter(0).daughter(1).pz    +daughter(1).pz, 2))"),   
    ),
    flags = cms.PSet(
       globalmu = cms.string("daughter(1).masterClone.isGlobalMuon"),
       trackermu = cms.string("daughter(1).masterClone.isTrackerMuon"),
       newID  = cms.string("daughter(1).masterClone.userInt('newID')"), 
       prlID  = cms.string("daughter(1).masterClone.userInt('prlID')"), 
       mvaIso = cms.string("daughter(1).masterClone.userInt('mvaIso')"), 
       pfid   = cms.string("daughter(1).masterClone.userInt('pfMuId')"), 
       id     = cms.string("daughter(1).masterClone.isGlobalMuon && daughter(1).masterClone.track.numberOfValidHits() > 10"), 
       id2012 = cms.string("daughter(1).masterClone.isGlobalMuon && daughter(1).masterClone.track.hitPattern.trackerLayersWithMeasurement > 5"), 
       l1trig   = cms.string("daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*', 1, 0).size() > 0 || " + 
                             "daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Mu13_Mu8_v*', 1, 0).size() > 0"),
       l2trig   = cms.string("daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*', 1, 0).size() > 0 || " + 
                             "daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Mu13_Mu8_v*', 1, 0).size() > 0"),
       l3trig   = cms.string("daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*', 1, 0).size() > 0 || " + 
                             "daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Mu13_Mu8_v*', 1, 0).size() > 0"),
       l1trig3e7 = cms.string("daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_TripleEle10_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l2trig3e7 = cms.string("daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_TripleEle10_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l3trig3e7 = cms.string("daughter(1).masterClone.triggerObjectMatchesByPath('HLT_TripleEle10_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l1trig3e8 = cms.string("daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l2trig3e8 = cms.string("daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l3trig3e8 = cms.string("daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
    )
)

process.zlletree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("zPlusLep"),
    cut = cms.string("abs(daughter(1).pdgId) == 11"),
    addRunLumiInfo = cms.bool(True),
    variables   = cms.PSet(
       zmass    = cms.string("daughter(0).mass"),
       met      = cms.string("daughter(0).masterClone.met"),
       pt       = cms.string("daughter(1).pt"),
       l1pt     = cms.string("daughter(0).masterClone.daughter(0).pt"),
       l1eta    = cms.string("daughter(0).masterClone.daughter(0).eta"),
       l1phi    = cms.string("daughter(0).masterClone.daughter(0).phi"),
       l1pdgId  = cms.string("daughter(0).masterClone.daughter(0).pdgId"),
       l2pt     = cms.string("daughter(0).masterClone.daughter(1).pt"),
       l2eta    = cms.string("daughter(0).masterClone.daughter(1).eta"),
       l2phi    = cms.string("daughter(0).masterClone.daughter(1).phi"),
       l2pdgId  = cms.string("daughter(0).masterClone.daughter(1).pdgId"),
       eta      = cms.string("daughter(1).eta"),
       phi      = cms.string("daughter(1).phi"),
       tkIso    = cms.string("daughter(1).masterClone.userFloat('tkZZ4L')"),
       ecalIso  = cms.string("daughter(1).masterClone.userFloat('ecalZZ4L')"),
       hcalIso  = cms.string("daughter(1).masterClone.userFloat('hcalZZ4L')"),
       rho      = cms.string("daughter(1).masterClone.userFloat('rhoEl')"),
       rhoAA    = cms.string("daughter(1).masterClone.userFloat('rhoElActiveArea')"),
       bdtIdYtDz = cms.string("daughter(1).masterClone.userFloat('bdtIDTrig')"),
       bdtIdNtDz = cms.string("daughter(1).masterClone.userFloat('bdtID')"),
       bdtIsoNtDz = cms.string("daughter(1).masterClone.userFloat('bdtIso')"),
       bdtID     = cms.string("daughter(1).masterClone.userFloat('bdtID')"),
       bdtIso    = cms.string("daughter(1).masterClone.userFloat('bdtIso')"),
       pfIsoChHad04       = cms.string("daughter(1).masterClone.userFloat('pfChHadIso04')"),
       pfIsoNHad04_NoEA   = cms.string("daughter(1).masterClone.userFloat('pfNHadIso04')"),
       pfIsoPhoton04_NoEA = cms.string("daughter(1).masterClone.userFloat('pfPhotonIso04')"),
       pfCombRelIso04EACorr = cms.string("daughter(1).masterClone.userFloat('pfCombRelIso04EACorr')"),
       pfCombRelIso04EACorr_WithFSR = cms.string("daughter(1).masterClone.userFloat('pfCombRelIso04EACorr_WithFSR')"),
       sip    = cms.string("daughter(1).masterClone.userFloat('sip')"),
       nmisshits = cms.string('daughter(1).masterClone.gsfTrack.trackerExpectedHitsInner.numberOfHits'), 
       numvertices = cms.string("daughter(0).masterClone.numvertices"),
       intimeSimVertices = cms.string("daughter(0).masterClone.numsimvertices"),
       numTrueInteractions = cms.string("daughter(0).masterClone.numTrueInteractions"),
       dr13 = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)"),
       dr23 = cms.string("deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)"),
       q13  = cms.string("daughter(0).daughter(0).charge + daughter(1).charge"),
       q23  = cms.string("daughter(0).daughter(1).charge + daughter(1).charge"),
       m13  = cms.string("sqrt(pow(daughter(0).daughter(0).energy+daughter(1).energy, 2) - " +
                         "     pow(daughter(0).daughter(0).px    +daughter(1).px, 2) -" +  
                         "     pow(daughter(0).daughter(0).py    +daughter(1).py, 2) -" +  
                         "     pow(daughter(0).daughter(0).pz    +daughter(1).pz, 2))"),   
       m23  = cms.string("sqrt(pow(daughter(0).daughter(1).energy+daughter(1).energy, 2) - "+
                         "     pow(daughter(0).daughter(1).px    +daughter(1).px, 2) -"+   
                         "     pow(daughter(0).daughter(1).py    +daughter(1).py, 2) -"+   
                         "     pow(daughter(0).daughter(1).pz    +daughter(1).pz, 2))"),   
    ),
    flags = cms.PSet(
       newID  = cms.string("daughter(1).masterClone.userInt('newID')"), 
       prlID  = cms.string("daughter(1).masterClone.userInt('prlID')"), 
       mvaIso = cms.string("daughter(1).masterClone.userInt('mvaIso')"), 
       mvaID  = cms.string("daughter(1).masterClone.userInt('mvaID')"), 
       id     = cms.string("test_bit(daughter(1).masterClone.electronID('cicTight'),0) == 1 && daughter(1).masterClone.gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1"),   
       l1trig   = cms.string("daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 1, 0).size() > 0 || " + 
                             "daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*', 1, 0).size() > 0"),
       l2trig   = cms.string("daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',1, 0).size() > 0 || " +
                             "daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*').size() > 0"),
       l3trig   = cms.string("daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 1, 0).size() > 0 || " +
                             "daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*', 1, 0).size() > 0"),
       l1trig3e7 = cms.string("daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_TripleEle10_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l2trig3e7 = cms.string("daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_TripleEle10_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l3trig3e7 = cms.string("daughter(1).masterClone.triggerObjectMatchesByPath('HLT_TripleEle10_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l1trig3e8 = cms.string("daughter(0).daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l2trig3e8 = cms.string("daughter(0).daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
       l3trig3e8 = cms.string("daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*',1,0).size() > 0"),
    )
)

# ========= FAKE RATE APPLICATION PART (Z + 2L) =================
process.diLepCRbare = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("looseLepCR looseLepCR"),
    cut = cms.string(SEL_ANY_CR_Z),
    checkCharge = cms.bool(False)
)

process.diLepCRnoFSR = cms.EDProducer("SkimEvent2LProducer",
    src = cms.InputTag("diLepCRbare"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    isMC = cms.bool(isMC),
)

process.diLepCR = process.zllAny.clone(
    zll = 'diLepCRnoFSR',
    photons = 'fsrPhotonsCR',
    isolationCut = 9e9, # don't reject non-isolated Z's
    label = cms.string("pfCombRelIso04EACorr_WithFSR"),
)

process.zx = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedZ1 diLepCR"),
    cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&"+
                     "deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02 &&"+
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&"+
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02"),
    checkCharge = cms.bool(False)
)


process.skimEventZX = cms.EDProducer("SkimEvent4LProducer",
    src = cms.InputTag("zx"),
    reducedPFCands = cms.InputTag("reducedPFCands"),
    jets = cms.InputTag("goodJets"),
    pfMet = cms.InputTag("pfMet"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    isMC = cms.bool(isMC),
    mcMatch = cms.InputTag(""),
    doMELA = cms.bool(True),
    doMassRes = cms.bool(True),
    doExtendedMassRes = cms.bool(True),
    energyForMELA = cms.double(E_LHC),                                      
    doswap = cms.bool(False) ## Leave the Z1 as is
)

process.anyZxTree = process.zz4lTree.clone( src = "skimEventZX")

process.skimEventZXsort1 = process.best4Lpass1.clone( src = "skimEventZX" )
process.bestZX = process.best4L.clone( src = "skimEventZXsort1")
process.skimEventZXcut1 = process.selectedZZs1.clone( src = "bestZX" )
process.skimEventZXcut2 = process.selectedZZs2.clone( src = "skimEventZXcut1" )
process.skimEventZXcut3 = process.selectedZZs3.clone( src = "skimEventZXcut2" )
process.skimEventZXcut4 = process.selectedZZs4.clone( src = "skimEventZXcut3" )
process.zxTree = process.zz4lTree.clone( src = "skimEventZXcut4")
process.zxTree.variables.l3pfIsoComb04EACorr_WithFSR = cms.string("z(1).masterClone.userFloat('pfCombRelIso04EACorr_WithFSR[0]')*lpt(1,0)")
process.zxTree.variables.l4pfIsoComb04EACorr_WithFSR = cms.string("z(1).masterClone.userFloat('pfCombRelIso04EACorr_WithFSR[1]')*lpt(1,1)")
process.skimEventZXSS = cms.EDFilter("SkimEvent4LSelector", 
        src = cms.InputTag("skimEventZX"), 
        cut = cms.string("lpdgId(1,0) == lpdgId(1,1)"), 
)
process.skimEventZXOS = process.skimEventZXSS.clone(cut = "lpdgId(1,0) == -lpdgId(1,1)")


# Setting up paths
skimseq = process.reskim
if SKIM_SEQUENCE not in ['', 'reskim']:
    skimseq = getattr(process, SKIM_SEQUENCE)
elif SKIM_SEQUENCE == "":
    process.dummyFilter = cms.EDFilter("HLTBool", result = cms.bool(True))
    skimseq = process.dummyFilter
if TRIGGER_FILTER:
    skimseq += getattr(process, TRIGGER_FILTER)

process.common = cms.Sequence(
    process.reboosting +
    process.postreboosting + 
    skimseq + 
    process.looseMuNoClean + process.looseMu +
    process.looseElNoClean + process.looseEl +
    process.looseLep +
    process.countLooseLep2 + 
    process.goodMu +
    process.goodEl +
    process.goodLep +
    process.slimPatJetsReCorr +
    process.goodJets +
    process.goodLL +
    process.zllAnyNoFSR +
    process.zllAny + 
    process.zll +
    process.bestZ
)

if DO_FSR_RECOVERY: process.common.replace(process.zllAnyNoFSR, process.zllAnyNoFSR + process.fsrPhotonSeq)

if False:
    process.load("WWAnalysis.AnalysisStep.genRecoMatcher_cfi")
    process.common.replace(process.zll, process.genRecoMatcher + process.zll)
    process.skimEvent4LNoArb.isSignal = cms.bool(True)
    process.skimEvent4LNoArb.mcMatch  = cms.InputTag("genRecoMatcher")

process.zzPathSeq = cms.Sequence( # make as sequence, so I can include in other sequences/paths
    process.common +
    process.oneZ +
    process.selectedZ1 +
    process.zz     +
    process.oneZZ  +
    process.skimEvent4LNoArb +  process.zz4lTreeNoArb +
    process.zzCombinatoric + 
    process.zz4lTree  
)
process.zzPath = cms.Path( process.zzPathSeq )
if FOUR_LEPTON_FILTER_PRE_Z:
    process.zzPath.replace(process.oneZ, process.countSequenceLLG + process.oneZ)
if FOUR_LEPTON_FILTER_POST_Z:
    process.zzPath.replace(process.zz, process.countSequenceLLG + process.zz)

process.count4lPath  = cms.Path(
    process.common +
    process.countGoodLep
)
process.countZZ1FSR = process.selectedZZs1.clone(
    src = process.zz4lTree.src,
    cut = cms.string("z(0).numberOfDaughters + z(1).numberOfDaughters == 5"),
)
process.countZZ2FSR = process.countZZ1FSR.clone(
    cut = cms.string("z(0).numberOfDaughters + z(1).numberOfDaughters == 6"),
)
process.zzPath_1FSR = cms.Path(process.zzPathSeq + process.countZZ1FSR)
process.zzPath_2FSR = cms.Path(process.zzPathSeq + process.countZZ2FSR)

if doVBF:
    from PhysicsTools.PatAlgos.tools.helpers import listModules
    modules = listModules(process.zzPath)
    process.vbfPath = cms.Path(sum(modules[1:], modules[0]))
    if ARBITRATE_EARLY:
        process.selectedZZsVBF1 = process.selectedZZs6.clone(cut = SEL_ZZ4L_VBF_STEP_1)
        process.selectedZZsVBF2 = process.selectedZZs6.clone(cut = SEL_ZZ4L_VBF_STEP_2, src ="selectedZZsVBF1")
        process.selectedZZsVBF3 = process.selectedZZs6.clone(cut = SEL_ZZ4L_VBF_STEP_3, src ="selectedZZsVBF2")
        process.selectedZZsVBF4 = process.selectedZZs6.clone(cut = SEL_ZZ4L_VBF_STEP_4, src ="selectedZZsVBF3")
        process.zz4lVBFTree = process.zz4lTree.clone(src = "selectedZZsVBF4")
        process.vbfPath.replace(process.selectedZZs6, process.selectedZZsVBF1+process.selectedZZsVBF2+process.selectedZZsVBF3+
                                                      process.selectedZZsVBF4+process.zz4lVBFTree)
    else:
        raise RuntimeError, "Not yet implemented for ARBITRATE_EARLY = False"

process.leptonPath = cms.Path(
    process.common +   
    process.muonTree + process.electronTree 
)
if DO_FSR_RECOVERY: process.leptonPath._seq += process.photonTree
if doVBF: process.leptonPath._seq += process.jetTree

process.zPath = cms.Path(
    process.common +   
    process.oneZ +
    process.zmmtree + process.zeetree 
)
process.countZ1FSR = process.selectedZ1.clone(
    src = cms.InputTag("selectedZ1"),
    cut = cms.string("numberOfDaughters > 2"),
)
process.countZ1eeFSR = process.countZ1FSR.clone(cut = "numberOfDaughters > 2 && abs(daughter(0).pdgId) == 11")
process.countZ1mmFSR = process.countZ1FSR.clone(cut = "numberOfDaughters > 2 && abs(daughter(0).pdgId) == 13")

process.countZ1FSRPath   = cms.Path( process.common + process.oneZ + process.selectedZ1 + process.countZ1FSR)
process.countZ1eeFSRPath = cms.Path( process.common + process.oneZ + process.selectedZ1 + process.countZ1eeFSR)
process.countZ1mmFSRPath = cms.Path( process.common + process.oneZ + process.selectedZ1 + process.countZ1mmFSR)

process.zlPath = cms.Path(
    process.common +
    process.bestZ +
    process.selectedZ1 +
    process.looseMuCR + 
    process.looseElCR + 
    process.looseLepCR +
    process.looseLepCRwithFSRIso +
    process.metVeto +
    process.lepMinFilter +
    ~process.lepMaxFilter +    
    process.zPlusLep     + process.zllmtree     + process.zlletree 
)

process.zllPathSeq = cms.Sequence(
    process.common +
    process.bestZ +
    process.selectedZ1 +
    process.looseMuCR + 
    process.looseElCR + 
    process.looseLepCR +
    process.diLepCRbare +
    process.diLepCRnoFSR +
    process.diLepCR +
    process.zx    +
    process.skimEventZX     +  process.anyZxTree
)
process.zxTreeSeq = cms.Sequence(
    process.skimEventZXsort1 + process.bestZX +  
    process.skimEventZXcut1 +
    process.skimEventZXcut2 +
    process.skimEventZXcut3 +
    process.skimEventZXcut4 +
    process.zxTree 
)
from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
process.zxTreeSeqOS = cms.Sequence(process.skimEventZXOS + cloneProcessingSnippet(process, process.zxTreeSeq, "OS"))
process.zxTreeSeqSS = cms.Sequence(process.skimEventZXSS + cloneProcessingSnippet(process, process.zxTreeSeq, "SS"))
process.skimEventZXsort1OS.src = "skimEventZXOS"
process.skimEventZXsort1SS.src = "skimEventZXSS"

process.zllPath   = cms.Path( process.zllPathSeq + process.zxTreeSeq )
process.zllSSPath = cms.Path( process.zllPathSeq + process.zxTreeSeqSS )
process.zllOSPath = cms.Path( process.zllPathSeq + process.zxTreeSeqOS )

if DO_FSR_RECOVERY:
    process.zlPath.replace( process.looseLepCR,  process.looseLepCR + process.fsrPhotonsCR)
    process.zllPath.replace(process.looseLepCR,  process.looseLepCR + process.fsrPhotonsCR)

### Paths with MC matching sequence
from WWAnalysis.AnalysisStep.zz4l.ptEtaFilters import adaEtaFilter, adaPtMinFilter
process.gen3FilterEta254PtMin = process.gen3FilterEta254PtMin5.clone(cut = cms.string(adaEtaFilter(2.5,2.4)+" && "+adaPtMinFilter(EL_PT_MIN,MU_PT_MIN)))
process.genSkim      = cms.Sequence(process.gen3RecoSeq + process.gen3FilterAny)
process.genSkimPtEta = cms.Sequence(process.genSkim + process.gen3FilterEta254PtMin)
process.ZZ_Any   = cms.Path(  process.genSkim                      + process.zzPath._seq )
process.ZZ_4Mu   = cms.Path(  process.genSkim + process.gen3ZZ4M   + process.zzPath._seq )
process.ZZ_4E    = cms.Path(  process.genSkim + process.gen3ZZ4E   + process.zzPath._seq )
process.ZZ_2E2Mu = cms.Path(  process.genSkim + process.gen3ZZ2E2M + process.zzPath._seq )
process.ZZ_GenPtEta_Any   = cms.Path(  process.genSkimPtEta                      + process.zzPath._seq  )
process.ZZ_GenPtEta_4Mu   = cms.Path(  process.genSkimPtEta + process.gen3ZZ4M   + process.zzPath._seq  )
process.ZZ_GenPtEta_4E    = cms.Path(  process.genSkimPtEta + process.gen3ZZ4E   + process.zzPath._seq  )
process.ZZ_GenPtEta_2E2Mu = cms.Path(  process.genSkimPtEta + process.gen3ZZ2E2M + process.zzPath._seq  )
process.ZZ_LepMonitor = cms.Path( process.gen3RecoSeq + process.gen3FilterAny + process.leptonPath._seq )

## Schedule without MC matching 
process.schedule = cms.Schedule(process.zzPath)
if addLeptonPath:
    process.schedule.extend([process.leptonPath, process.count4lPath])
if addZPath:
    process.schedule.extend([process.zPath])

process.schedule.extend([process.zlPath])
process.schedule.extend([process.zllPath, process.zllSSPath, process.zllOSPath])
 
if doSyncPaths:
    ### make paths with reco classification
    from WWAnalysis.AnalysisStep.zz4l.recoFinalStateClassifiers_cff import makeSplittedPaths4L
    makeSplittedPaths4L(process, 'zzPath', TRIGGER_FILTER)
    ### add them to schedule
    process.schedule.extend([ process.zzPath_4E, process.zzPath_4M, process.zzPath_2E2M ])
    ##process.schedule.extend([ process.zzPath_4E_3Path, process.zzPath_4M_3Path ]) # not commissioned yet with FSR
if doSyncPaths:
    from WWAnalysis.AnalysisStep.zz4l.recoFinalStateClassifiers_cff import makeSplittedPaths4L
    makeSplittedPaths4L(process, 'zzPath_1FSR', TRIGGER_FILTER, doThreePathLogic=False)
    makeSplittedPaths4L(process, 'zzPath_2FSR', TRIGGER_FILTER, doThreePathLogic=False)
    process.schedule.extend([ process.countZ1FSRPath, process.countZ1eeFSRPath, process.countZ1mmFSRPath] )
    process.schedule.extend([ process.zzPath_1FSR, process.zzPath_2FSR ])
    process.schedule.extend([ process.zzPath_1FSR_4E, process.zzPath_1FSR_4M, process.zzPath_1FSR_2E2M ])
    process.schedule.extend([ process.zzPath_2FSR_4E, process.zzPath_2FSR_4M, process.zzPath_2FSR_2E2M ])
if doVBF and doSyncPaths:
    process.schedule.extend([process.vbfPath])
    from WWAnalysis.AnalysisStep.zz4l.recoFinalStateClassifiers_cff import makeSplittedPaths4L
    makeSplittedPaths4L(process, 'vbfPath', TRIGGER_FILTER, doThreePathLogic=False)
    process.schedule.extend([ process.vbfPath_4E, process.vbfPath_4M, process.vbfPath_2E2M ])

## Add to schedules paths with MC matching
if False and isMC:
    process.schedule.extend([ process.ZZ_Any, process.ZZ_4Mu, process.ZZ_4E, process.ZZ_2E2Mu ])
    process.schedule.extend([ process.ZZ_GenPtEta_Any, process.ZZ_GenPtEta_4Mu, process.ZZ_GenPtEta_4E, process.ZZ_GenPtEta_2E2Mu ])

process.TFileService = cms.Service("TFileService", fileName = cms.string("hzzTree.root"))

if False:
    process.out = cms.OutputModule("PoolOutputModule", 
                fileName = cms.untracked.string("hzzEventsMissing.root"),
                #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
                outputCommands = cms.untracked.vstring("keep *", 
                        "drop *_*_*_Tree", 
                        "keep *_boostedElectronsStep2_*_Tree", "keep *_boostedMuonsStep2_*_Tree", 
                        "keep *_looseMu*_*_*", "keep *_looseEl*_*_*", 
                        "keep *_goodMu*_*_*", "keep *_goodEl*_*_*", 
                        "keep *_middleMu*_*_*", "keep *_middleEl*_*_*", 
                        "keep *_zllAny_*_*", "keep *_bestZ_*_*", 
                        "keep *_goodJets_*_*",
                        "keep *_TriggerResults_*_Tree"),
                )
    process.end = cms.EndPath(process.out)
    process.schedule.extend([process.end])
#process.selectedZZs6.cut = ""
#process.TFileService.fileName = "hzzTree.NoKDCut.root"
#process.source.eventsToProcess = cms.untracked.VEventRange('1:45744','1:35912','1:24264','1:26251')
