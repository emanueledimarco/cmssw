import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.yieldProducer_cfi import *
from WWAnalysis.AnalysisStep.cutPSets_cfi import *

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        'file:hypoEvents.root'
#         'file:/data/mwlebour/WW_414_SKIM_V04_STEP2_V00/101160/ggToH160toWWto2L2Nu_1_1_jGU.root'
#         'file:/nfs/bluearc/group/trees/hww/R414_S1_V06pre2_S2_V02_S3_V00/023/TTJetsMad_1_1_aDa.root'
    ), 
    maxEvents   = cms.int32(-1),
    outputEvery = cms.uint32(1000),
)
# from glob import glob
# process.fwliteInput.fileNames = [ 'file:%s'%x for x in glob('/nfs/bluearc/group/trees/hww/R414_S1_V06pre2_S2_V02_S3_V00/023/*.root') ]
# process.fwliteInput.fileNames = [ 'file:%s'%x for x in glob('/nfs/bluearc/group/trees/hww/R414_S1_V06pre2_S2_V02_S3_V00/101160/*.root') ]
# process.fwliteInput.fileNames = [ 'file:%s'%x for x in glob('/nfs/bluearc/group/trees/hww/WW_414_SKIM_V00/RMMEFOLDER/*.root') ]

# import commands
# myDir = "/castor/cern.ch/user/m/mangano/hww/skims/WW_413_SKIM_V00/id101160.GGtoH160toWWto2L2Nu/"
# process.fwliteInput.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 

    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('hists.root'),
)

process.eventHists = FWLiteParams.clone()
process.eventHists.sampleName = cms.string("101160.ggToH160toWWto2L2Nu")
# process.eventHists.doNMinus1 = False
# process.eventHists.doByCuts  = False

# Setup
addMassDependentCuts(process.eventHists.hypotheses.wwelel0.cuts,hReOpt160)

#MonteCarlo     SingleMuon     DoubleMuon     MuEG           DoubleElectron 
switchTrigger(process.eventHists.hypotheses.wwelel0.cuts,MC)
# switchTrigger(process.eventHists.hypotheses.wwelel0,SingleMuon)
# switchTrigger(process.eventHists.hypotheses.wwelel0,DoubleMuon)
# switchTrigger(process.eventHists.hypotheses.wwelel0,MuEG)
# switchTrigger(process.eventHists.hypotheses.wwelel0,DoubleElectron)

process.eventHists.hypotheses.wwmumuIDLHT = cms.PSet(src = cms.InputTag("wwmumuIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmuIDLHT = cms.PSet(src = cms.InputTag("wwelmuIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuelIDLHT = cms.PSet(src = cms.InputTag("wwmuelIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelelIDLHT = cms.PSet(src = cms.InputTag("wwelelIDLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuIDLHT.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelIDLHT.cuts)

process.eventHists.hypotheses.wwmumuCONVLHT = cms.PSet(src = cms.InputTag("wwmumuCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmuCONVLHT = cms.PSet(src = cms.InputTag("wwelmuCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuelCONVLHT = cms.PSet(src = cms.InputTag("wwmuelCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelelCONVLHT = cms.PSet(src = cms.InputTag("wwelelCONVLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuCONVLHT.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelCONVLHT.cuts)

process.eventHists.hypotheses.wwmumuISOLHT = cms.PSet(src = cms.InputTag("wwmumuISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmuISOLHT = cms.PSet(src = cms.InputTag("wwelmuISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuelISOLHT = cms.PSet(src = cms.InputTag("wwmuelISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelelISOLHT = cms.PSet(src = cms.InputTag("wwelelISOLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuISOLHT.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelISOLHT.cuts)

process.eventHists.hypotheses.wwmumuIPLHT = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmuIPLHT = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuelIPLHT = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelelIPLHT = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuIPLHT.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelIPLHT.cuts)

process.eventHists.hypotheses.wwmumuCONVPFLHT = cms.PSet(src = cms.InputTag("wwmumuCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmuCONVPFLHT = cms.PSet(src = cms.InputTag("wwelmuCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuelCONVPFLHT = cms.PSet(src = cms.InputTag("wwmuelCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelelCONVPFLHT = cms.PSet(src = cms.InputTag("wwelelCONVPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuCONVPFLHT.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelCONVPFLHT.cuts)

process.eventHists.hypotheses.wwmumuISOPFLHT = cms.PSet(src = cms.InputTag("wwmumuISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmuISOPFLHT = cms.PSet(src = cms.InputTag("wwelmuISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuelISOPFLHT = cms.PSet(src = cms.InputTag("wwmuelISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelelISOPFLHT = cms.PSet(src = cms.InputTag("wwelelISOPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuISOPFLHT.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelISOPFLHT.cuts)

process.eventHists.hypotheses.wwmumuIPPFLHT = cms.PSet(src = cms.InputTag("wwmumuIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelmuIPPFLHT = cms.PSet(src = cms.InputTag("wwelmuIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwmuelIPPFLHT = cms.PSet(src = cms.InputTag("wwmuelIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
process.eventHists.hypotheses.wwelelIPPFLHT = cms.PSet(src = cms.InputTag("wwelelIPPFLHT"), cuts = cloneVPSet(process.eventHists.hypotheses.wwelel0.cuts))
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuIPPFLHT.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelIPPFLHT.cuts)

process.eventHists.hypotheses.wwmumuTopBMC = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelmuTopBMC = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmuelTopBMC = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelelTopBMC = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmumuTopBMC.cuts.append(bTag.clone())
process.eventHists.hypotheses.wwelmuTopBMC.cuts.append(bTag.clone())
process.eventHists.hypotheses.wwmuelTopBMC.cuts.append(bTag.clone())
process.eventHists.hypotheses.wwelelTopBMC.cuts.append(bTag.clone())
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuTopBMC.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelTopBMC.cuts)

process.eventHists.hypotheses.wwmumuTopMuMC = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelmuTopMuMC = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmuelTopMuMC = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelelTopMuMC = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmumuTopMuMC.cuts.append(softMu.clone())
process.eventHists.hypotheses.wwelmuTopMuMC.cuts.append(softMu.clone())
process.eventHists.hypotheses.wwmuelTopMuMC.cuts.append(softMu.clone())
process.eventHists.hypotheses.wwelelTopMuMC.cuts.append(softMu.clone())
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuTopMuMC.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelTopMuMC.cuts)

process.eventHists.hypotheses.wwmumuTopTagMC = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelmuTopTagMC = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmuelTopTagMC = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelelTopTagMC = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmumuTopTagMC.cuts.append(topTag.clone())
process.eventHists.hypotheses.wwelmuTopTagMC.cuts.append(topTag.clone())
process.eventHists.hypotheses.wwmuelTopTagMC.cuts.append(topTag.clone())
process.eventHists.hypotheses.wwelelTopTagMC.cuts.append(topTag.clone())
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuTopTagMC.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelTopTagMC.cuts)

process.eventHists.hypotheses.wwmumuTopB = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelmuTopB = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmuelTopB = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelelTopB = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmumuTopB.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwelmuTopB.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwmuelTopB.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwelelTopB.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwmumuTopB.cuts.append(bTag.clone())
process.eventHists.hypotheses.wwelmuTopB.cuts.append(bTag.clone())
process.eventHists.hypotheses.wwmuelTopB.cuts.append(bTag.clone())
process.eventHists.hypotheses.wwelelTopB.cuts.append(bTag.clone())
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuTopB.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelTopB.cuts)

process.eventHists.hypotheses.wwmumuTopMu = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelmuTopMu = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmuelTopMu = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelelTopMu = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmumuTopMu.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwelmuTopMu.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwmuelTopMu.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwelelTopMu.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwmumuTopMu.cuts.append(softMu.clone())
process.eventHists.hypotheses.wwelmuTopMu.cuts.append(softMu.clone())
process.eventHists.hypotheses.wwmuelTopMu.cuts.append(softMu.clone())
process.eventHists.hypotheses.wwelelTopMu.cuts.append(softMu.clone())
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuTopMu.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelTopMu.cuts)

process.eventHists.hypotheses.wwmumuTopTag = cms.PSet(src = cms.InputTag("wwmumuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelmuTopTag = cms.PSet(src = cms.InputTag("wwelmuIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmuelTopTag = cms.PSet(src = cms.InputTag("wwmuelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwelelTopTag = cms.PSet(src = cms.InputTag("wwelelIPLHT"), cuts = cloneVPSet(ttBar))
process.eventHists.hypotheses.wwmumuTopTag.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwelmuTopTag.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwmuelTopTag.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwelelTopTag.cuts.append(highPtB.clone())
process.eventHists.hypotheses.wwmumuTopTag.cuts.append(topTag.clone())
process.eventHists.hypotheses.wwelmuTopTag.cuts.append(topTag.clone())
process.eventHists.hypotheses.wwmuelTopTag.cuts.append(topTag.clone())
process.eventHists.hypotheses.wwelelTopTag.cuts.append(topTag.clone())
switchToOppoFlavor(process.eventHists.hypotheses.wwelmuTopTag.cuts)
switchToOppoFlavor(process.eventHists.hypotheses.wwmuelTopTag.cuts)

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

delattr(process.eventHists.hypotheses,"wwelel0")

