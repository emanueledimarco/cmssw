import FWCore.ParameterSet.Config as cms
from FWCore.PythonUtilities.LumiList import LumiList
import os

from WWAnalysis.AnalysisStep.yieldProducer_cfi import *
from WWAnalysis.AnalysisStep.cutPSets_cfi import *
from WWAnalysis.AnalysisStep.cutMassDependPSets_cfi import *
from WWAnalysis.AnalysisStep.pileupReweighting_cfi import reWeightVector

process = cms.Process("Test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

dataIsLocal=RMMEISLOCAL
# dataIsLocal = False
# dataIsLocal = True

isGluGlu = RMMEGLUGLU
# isGluGlu = False
# isGluGlu = True

isMC = RMMEMC
# isMC = True
# isMC = False

process.source = cms.Source("PoolSource")
if dataIsLocal:
    process.source.fileNames = cms.untracked.vstring(RMMEINPUTFILES)
else:
    process.source.fileNames = cms.untracked.vstring(
        'file:RMMEFN'
#         'file:hypoEvents.root'                                                                                             
#         'file:/nfs/bluearc/group/trees/hww/R414_S1_V06_S2_V03/ggToH160toWWto2L2Nu.step2.1.root'                            
#         'file:/nfs/bluearc/group/trees/hww/R414_S1_V06_S2_V03_Reduced/101160/ggToH160toWWto2L2Nu_1_1_ZAS.root'
#         'file:/data/mwlebour/WW_414_SKIM_V04_STEP2_V00/101160/ggToH160toWWto2L2Nu_1_1_jGU.root'                            
    )

 
# from glob import glob
# process.source.fileNames = [ 'file:%s'%x for x in glob('/nfs/bluearc/group/trees/hww/R414_S1_V06_S2_V03/*.root') ]

# import commands
# myDir = "/castor/cern.ch/user/m/mangano/hww/skims/WW_413_SKIM_V00/id101160.GGtoH160toWWto2L2Nu/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) #RMME
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RMMEFN')
#     fileName = cms.string('hists.root')
)



process.eventHists = cms.EDAnalyzer("CreateEventHistsEDMWrapped",
    FWLiteParams.clone(),
)

if isMC:
    # PU re-weighting:
    process.eventHists.puWeights = reWeightVector[:]
    process.eventHists.puLabel = cms.InputTag("addPileupInfo")
else:
    lumis = LumiList(filename = 'RMMEJSON')
#     lumis = LumiList(filename = os.getenv('CMSSW_BASE')+'/src/WWAnalysis/Misc/Jsons/certifiedUCSD.json')
    process.eventHists.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.eventHists.lumisToProcess = lumis.getCMSSWString().split(',')

process.eventHists.sampleName = cms.string("RMMENUM_RMMENAME")
# process.eventHists.sampleName = cms.string("101160_ggToH160toWWto2L2Nu")
process.eventHists.doNMinus1 = False
# process.eventHists.doByCuts  = False
# process.eventHists.printSummary  = True


# RMME
# do this here or in the step 2 files?
if isGluGlu:
    process.higgsPt = cms.EDProducer("HWWKFactorProducer",
        genParticlesTag = cms.InputTag("onlyHiggsGen"),
        inputFilename = cms.untracked.string("HiggsAnalysis/HiggsToWW2Leptons/data/kfactors_Std/kfactors_mhRMMEMASS_renRMMEMASS_facRMMEMASS.dat")
#         inputFilename = cms.untracked.string("HiggsAnalysis/HiggsToWW2Leptons/data/kfactors_Std/kfactors_mh160_ren160_fac160.dat")
        ProcessID = cms.untracked.int32(10010),
        Debug =cms.untracked.bool(False)
    )
    process.eventHists.ptWeight = cms.InputTag("higgsPt")
    process.p = cms.Path(process.higgsPt * process.eventHists)
else:
    process.p = cms.Path(process.eventHists)
# process.p = cms.Path(process.eventHists)

addMassDependentCuts(process.eventHists.hypotheses.wwelel0.cuts,hGammaMRRMMEMASS)
# addMassDependentCuts(process.eventHists.hypotheses.wwelel0.cuts,hGammaMR160)
# addMassDependentCuts(process.eventHists.hypotheses.wwelel0.cuts,hThursdayRMMEMASS)
# addMassDependentCuts(process.eventHists.hypotheses.wwelel0.cuts,hThursday160)

#MonteCarlo     SingleMuon     DoubleMuon     MuEG           DoubleElectron 
switchTrigger(process.eventHists.hypotheses.wwelel0.cuts,RMMESAMPLE)
# switchTrigger(process.eventHists.hypotheses.wwelel0.cuts,MC)
# switchTrigger(process.eventHists.hypotheses.wwelel0.cuts,SingleMuon)
# switchTrigger(process.eventHists.hypotheses.wwelel0.cuts,DoubleMuon)
# switchTrigger(process.eventHists.hypotheses.wwelel0.cuts,MuEG)
#switchTrigger(process.eventHists.hypotheses.wwelel0.cuts,DoubleElectron)

process.eventHists.hypotheses.wwmumu0 = cms.PSet(src = cms.InputTag("wwmumu0"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmu0 = cms.PSet(src = cms.InputTag("wwelmu0"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuel0 = cms.PSet(src = cms.InputTag("wwmuel0"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmu0.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuel0.cuts)
# switchToTrailingElectron(process.eventHists.hypotheses.wwmuel0.cuts)
# switchToTrailingElectron(process.eventHists.hypotheses.wwelel0.cuts)
# injectTightIsolationCut(process.eventHists.hypotheses.wwmumu0.cuts)
# injectTightIsolationCut(process.eventHists.hypotheses.wwelmu0.cuts)
# injectTightIsolationCut(process.eventHists.hypotheses.wwmuel0.cuts)
# injectTightIsolationCut(process.eventHists.hypotheses.wwelel0.cuts)
# injectIPCut(process.eventHists.hypotheses.wwmumu0.cuts)
# injectIPCut(process.eventHists.hypotheses.wwelmu0.cuts)
# injectIPCut(process.eventHists.hypotheses.wwmuel0.cuts)
# injectIPCut(process.eventHists.hypotheses.wwelel0.cuts)

process.eventHists.hypotheses.wwmumuIPLHT = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuIPLHT = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelIPLHT = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelIPLHT = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

process.eventHists.hypotheses.wwmumuIDLHT = cms.PSet(src = cms.InputTag("wwmumuIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuIDLHT = cms.PSet(src = cms.InputTag("wwelmuIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelIDLHT = cms.PSet(src = cms.InputTag("wwmuelIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelIDLHT = cms.PSet(src = cms.InputTag("wwelelIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

process.eventHists.hypotheses.wwmumuCONVLHT = cms.PSet(src = cms.InputTag("wwmumuCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuCONVLHT = cms.PSet(src = cms.InputTag("wwelmuCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelCONVLHT = cms.PSet(src = cms.InputTag("wwmuelCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelCONVLHT = cms.PSet(src = cms.InputTag("wwelelCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

process.eventHists.hypotheses.wwmumuISOLHT = cms.PSet(src = cms.InputTag("wwmumuISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuISOLHT = cms.PSet(src = cms.InputTag("wwelmuISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelISOLHT = cms.PSet(src = cms.InputTag("wwmuelISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelISOLHT = cms.PSet(src = cms.InputTag("wwelelISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

process.eventHists.hypotheses.wwmumuIPMerge = cms.PSet(src = cms.InputTag("wwmumuIPMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuIPMerge = cms.PSet(src = cms.InputTag("wwelmuIPMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelIPMerge = cms.PSet(src = cms.InputTag("wwmuelIPMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelIPMerge = cms.PSet(src = cms.InputTag("wwelelIPMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

process.eventHists.hypotheses.wwmumuIDMerge = cms.PSet(src = cms.InputTag("wwmumuIDMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuIDMerge = cms.PSet(src = cms.InputTag("wwelmuIDMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelIDMerge = cms.PSet(src = cms.InputTag("wwmuelIDMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelIDMerge = cms.PSet(src = cms.InputTag("wwelelIDMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

process.eventHists.hypotheses.wwmumuCONVMerge = cms.PSet(src = cms.InputTag("wwmumuCONVMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuCONVMerge = cms.PSet(src = cms.InputTag("wwelmuCONVMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelCONVMerge = cms.PSet(src = cms.InputTag("wwmuelCONVMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelCONVMerge = cms.PSet(src = cms.InputTag("wwelelCONVMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

process.eventHists.hypotheses.wwmumuISOMerge = cms.PSet(src = cms.InputTag("wwmumuISOMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
process.eventHists.hypotheses.wwelmuISOMerge = cms.PSet(src = cms.InputTag("wwelmuISOMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
process.eventHists.hypotheses.wwmuelISOMerge = cms.PSet(src = cms.InputTag("wwmuelISOMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
process.eventHists.hypotheses.wwelelISOMerge = cms.PSet(src = cms.InputTag("wwelelISOMerge"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

# process.eventHists.hypotheses.wwmumuCONVPFLHT = cms.PSet(src = cms.InputTag("wwmumuCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
# process.eventHists.hypotheses.wwelmuCONVPFLHT = cms.PSet(src = cms.InputTag("wwelmuCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
# process.eventHists.hypotheses.wwmuelCONVPFLHT = cms.PSet(src = cms.InputTag("wwmuelCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
# process.eventHists.hypotheses.wwelelCONVPFLHT = cms.PSet(src = cms.InputTag("wwelelCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

# process.eventHists.hypotheses.wwmumuISOPFLHT = cms.PSet(src = cms.InputTag("wwmumuISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
# process.eventHists.hypotheses.wwelmuISOPFLHT = cms.PSet(src = cms.InputTag("wwelmuISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
# process.eventHists.hypotheses.wwmuelISOPFLHT = cms.PSet(src = cms.InputTag("wwmuelISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
# process.eventHists.hypotheses.wwelelISOPFLHT = cms.PSet(src = cms.InputTag("wwelelISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

# process.eventHists.hypotheses.wwmumuIPPFLHT = cms.PSet(src = cms.InputTag("wwmumuIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmumu0.cuts))
# process.eventHists.hypotheses.wwelmuIPPFLHT = cms.PSet(src = cms.InputTag("wwelmuIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelmu0.cuts))
# process.eventHists.hypotheses.wwmuelIPPFLHT = cms.PSet(src = cms.InputTag("wwmuelIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwmuel0.cuts))
# process.eventHists.hypotheses.wwelelIPPFLHT = cms.PSet(src = cms.InputTag("wwelelIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))

# ttBarEstimate = RMMETTBAR
ttBarEstimate = False
if ttBarEstimate:
    switchTrigger(ttBar,RMMESAMPLE)
#     switchTrigger(ttBar,MC)
#     switchTrigger(ttBar,SingleMuon)
#     switchTrigger(ttBar,DoubleMuon)
#     switchTrigger(ttBar,MuEG)
#     switchTrigger(ttBar,DoubleElectron)
    
    mumuTopTemp = cloneVPSet(ttBar)
    
    elmuTopTemp = cloneVPSet(ttBar)
    switchToOppoFlavor(elmuTopTemp)
    
    muelTopTemp = cloneVPSet(ttBar)
    switchToOppoFlavor(muelTopTemp)
    switchToTrailingElectron(muelTopTemp)
    
    elelTopTemp = cloneVPSet(ttBar)
    switchToTrailingElectron(elelTopTemp)

    # if RMMESAMPLE == "MC" :
    if "MC" == "MCasdf" :
        process.eventHists.hypotheses.wwmumuTopBMC = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(mumuTopTemp))
        process.eventHists.hypotheses.wwelmuTopBMC = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(elmuTopTemp))
        process.eventHists.hypotheses.wwmuelTopBMC = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(muelTopTemp))
        process.eventHists.hypotheses.wwelelTopBMC = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(elelTopTemp))
        process.eventHists.hypotheses.wwmumuTopBMC.cuts.append(bTag.clone())
        process.eventHists.hypotheses.wwelmuTopBMC.cuts.append(bTag.clone())
        process.eventHists.hypotheses.wwmuelTopBMC.cuts.append(bTag.clone())
        process.eventHists.hypotheses.wwelelTopBMC.cuts.append(bTag.clone())
        
        process.eventHists.hypotheses.wwmumuTopMuMC = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(mumuTopTemp))
        process.eventHists.hypotheses.wwelmuTopMuMC = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(elmuTopTemp))
        process.eventHists.hypotheses.wwmuelTopMuMC = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(muelTopTemp))
        process.eventHists.hypotheses.wwelelTopMuMC = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(elelTopTemp))
        process.eventHists.hypotheses.wwmumuTopMuMC.cuts.append(softMu.clone())
        process.eventHists.hypotheses.wwelmuTopMuMC.cuts.append(softMu.clone())
        process.eventHists.hypotheses.wwmuelTopMuMC.cuts.append(softMu.clone())
        process.eventHists.hypotheses.wwelelTopMuMC.cuts.append(softMu.clone())
        
        process.eventHists.hypotheses.wwmumuTopTagMC = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(mumuTopTemp))
        process.eventHists.hypotheses.wwelmuTopTagMC = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(elmuTopTemp))
        process.eventHists.hypotheses.wwmuelTopTagMC = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(muelTopTemp))
        process.eventHists.hypotheses.wwelelTopTagMC = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(elelTopTemp))
        process.eventHists.hypotheses.wwmumuTopTagMC.cuts.append(topTag.clone())
        process.eventHists.hypotheses.wwelmuTopTagMC.cuts.append(topTag.clone())
        process.eventHists.hypotheses.wwmuelTopTagMC.cuts.append(topTag.clone())
        process.eventHists.hypotheses.wwelelTopTagMC.cuts.append(topTag.clone())
    
        process.eventHists.hypotheses.wwmumuTopTag = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(mumuTopTemp))
        process.eventHists.hypotheses.wwelmuTopTag = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(elmuTopTemp))
        process.eventHists.hypotheses.wwmuelTopTag = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(muelTopTemp))
        process.eventHists.hypotheses.wwelelTopTag = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(elelTopTemp))
        process.eventHists.hypotheses.wwmumuTopTag.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwelmuTopTag.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwmuelTopTag.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwelelTopTag.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwmumuTopTag.cuts.append(topTag.clone())
        process.eventHists.hypotheses.wwelmuTopTag.cuts.append(topTag.clone())
        process.eventHists.hypotheses.wwmuelTopTag.cuts.append(topTag.clone())
        process.eventHists.hypotheses.wwelelTopTag.cuts.append(topTag.clone())
     
        process.eventHists.hypotheses.wwmumuTopB = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(mumuTopTemp))
        process.eventHists.hypotheses.wwelmuTopB = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(elmuTopTemp))
        process.eventHists.hypotheses.wwmuelTopB = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(muelTopTemp))
        process.eventHists.hypotheses.wwelelTopB = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(elelTopTemp))
        process.eventHists.hypotheses.wwmumuTopB.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwelmuTopB.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwmuelTopB.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwelelTopB.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwmumuTopB.cuts.append(bTag.clone())
        process.eventHists.hypotheses.wwelmuTopB.cuts.append(bTag.clone())
        process.eventHists.hypotheses.wwmuelTopB.cuts.append(bTag.clone())
        process.eventHists.hypotheses.wwelelTopB.cuts.append(bTag.clone())
     
        process.eventHists.hypotheses.wwmumuTopMu = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(mumuTopTemp))
        process.eventHists.hypotheses.wwelmuTopMu = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(elmuTopTemp))
        process.eventHists.hypotheses.wwmuelTopMu = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(muelTopTemp))
        process.eventHists.hypotheses.wwelelTopMu = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(elelTopTemp))
        process.eventHists.hypotheses.wwmumuTopMu.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwelmuTopMu.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwmuelTopMu.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwelelTopMu.cuts.append(highPtB.clone())
        process.eventHists.hypotheses.wwmumuTopMu.cuts.append(NoBtag.clone())
        process.eventHists.hypotheses.wwelmuTopMu.cuts.append(NoBtag.clone())
        process.eventHists.hypotheses.wwmuelTopMu.cuts.append(NoBtag.clone())
        process.eventHists.hypotheses.wwelelTopMu.cuts.append(NoBtag.clone())
        process.eventHists.hypotheses.wwmumuTopMu.cuts.append(softMu.clone())
        process.eventHists.hypotheses.wwelmuTopMu.cuts.append(softMu.clone())
        process.eventHists.hypotheses.wwmuelTopMu.cuts.append(softMu.clone())
        process.eventHists.hypotheses.wwelelTopMu.cuts.append(softMu.clone())
     
        process.eventHists.histParams.nBJetsOver = cms.PSet(
            variable = cms.string("bTaggedJetsOver(30.0,2.1)"),
            nbins = cms.uint32(10),
            low = cms.double(-0.5),
            high = cms.double(9.5),
            xtitle = cms.string("N_{b-Jets}^{over}")
            )
    
        process.eventHists.histParams.nBJetsUnder = cms.PSet(
            variable = cms.string("bTaggedJetsUnder(30.0,2.1)"),
            nbins = cms.uint32(10),
            low = cms.double(-0.5),
            high = cms.double(9.5),
            xtitle = cms.string("N_{b-Jets}^{under}")
            )
    
        process.eventHists.histParams.nBJets = cms.PSet(
            variable = cms.string("bTaggedJetsOver(30.0,2.1) + bTaggedJetsUnder(30.0,2.1)"),
            nbins = cms.uint32(10),
            low = cms.double(-0.5),
            high = cms.double(9.5),
            xtitle = cms.string("N_{b-Jets}^{over}")
            )
    
# delattr(process.eventHists.hypotheses,"wwelel0")
# delattr(process.eventHists.hypotheses,"wwmuel0")
# delattr(process.eventHists.hypotheses,"wwelmu0")
# delattr(process.eventHists.hypotheses,"wwmumu0")

