import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from glob import glob
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/DoubleElectron2012_Z/hzz4lSkim_137_1_8MW.root',
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/DoubleElectron2012_Z/hzz4lSkim_142_1_cqT.root',
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/DoubleElectron2012_Z/hzz4lSkim_37_1_fBh.root',
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/DoubleElectron2012_Z/hzz4lSkim_42_1_8Yd.root',
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/SingleMu2012_Z/hzz4lSkim_137_0_MS7.root',
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/SingleMu2012_Z/hzz4lSkim_142_0_c3x.root',
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/SingleMu2012_Z/hzz4lSkim_37_1_538.root',
	'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/hzz/skims/52X/2012_05_01/SingleMu2012_Z/hzz4lSkim_42_0_C7u.root',
    )
)
#process.source.fileNames = [ 'file:/afs/cern.ch/work/g/gpetrucc/HZZ/CMSSW_5_2_4_patch4/src/WWAnalysis/SkimStep/test/hzz4lSkim.root' ]
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("WWAnalysis.AnalysisStep.hzz4l_selection_cff")
process.load("WWAnalysis.AnalysisStep.tp_lepton_effs_cff")

from WWAnalysis.AnalysisStep.tp_lepton_effs_cff import *

process.tagMuons.cut = cms.string(
    "pt > 20 && abs(eta) < 2.1 && " + 
    "userInt('pfMuId') && "+
    "muonID('GlobalMuonPromptTight') && numberOfMatchedStations > 1 && "+
    "innerTrack.hitPattern.numberOfValidPixelHits > 0 && "+
    "track.hitPattern.trackerLayersWithMeasurement > 5 && "
    "userFloat('sip') < 4 && "+
    "userFloat('pfChHadIso04')/pt < 0.1 && "+
    "("+HLT_Any1M+")"
)
process.tagElectrons.cut = cms.string(
    "pt > 20 && abs(eta) < 2.5 && " + 
    "test_bit(electronID('cicSuperTight'),0) == 1 && "+
    "gsfTrack.trackerExpectedHitsInner.numberOfHits == 0 && "+
    "userFloat('sip') < 4 && "+
    "userFloat('pfChHadIso04')/pt < 0.1 && "
    "("+"||".join([HLT_Any1E,EleTriggers.HLT_Ele32_SC17_TnP_Ele32Leg.value(),EleTriggers.HLT_Ele17_Ele8_TnP_Ele17Leg.value()])+")"
)


process.tpTreeMuMu.variables.bdtIdDz  = cms.string("userFloat('bdtidnontrigDZ')")
process.tpTreeMuMu.variables.bdtIsoDz  = cms.string("userFloat('bdtisonontrigDZ')")
process.tpTreeMuMu.variables.bdtIsoPnP = cms.string("userFloat('bdtisonontrigPFNOPU')")
process.tpTreeMuMu.variables.tkIso    = cms.string("userFloat('tkZZ4L')")
process.tpTreeMuMu.variables.ecalIso  = cms.string("userFloat('ecalZZ4L')")
process.tpTreeMuMu.variables.hcalIso  = cms.string("userFloat('hcalZZ4L')")
process.tpTreeMuMu.variables.rho      = cms.string("userFloat('rhoMu')")
process.tpTreeMuMu.variables.pfIsoChHad04       = cms.string("userFloat('muonPFIsoChHad04')")
process.tpTreeMuMu.variables.pfIsoNHad04_NoEA   = cms.string("userFloat('muonPFIsoNHad04')")
process.tpTreeMuMu.variables.pfIsoPhoton04_NoEA = cms.string("userFloat('muonPFIsoPhoton04')")
process.tpTreeMuMu.variables.pfIsoComb04EACorr  = cms.string("userFloat('pfCombIso04EACorr')")
process.tpTreeMuMu.variables.pfIsoCombRel04EACorr  = cms.string("userFloat('pfCombIso04EACorr')/pt")
process.tpTreeMuMu.variables.sip = cms.string("userFloat('sip')")
process.tpTreeMuMu.flags.newID = cms.string("userInt('newID')")
process.tpTreeMuMu.flags.pfID  = cms.string("userInt('pfMuId')")
process.tpTreeMuMu.flags.prlID = cms.string("userInt('prlID')")
process.tpTreeMuMu.flags.prlID2012 = cms.string("isGlobalMuon && track.hitPattern.trackerLayersWithMeasurement > 5")

process.tpTreeElEl.variables.bdtIdYtDz  = cms.string("userFloat('bdttrig')")
process.tpTreeElEl.variables.bdtIdNtDz  = cms.string("userFloat('bdtnontrig')")
process.tpTreeElEl.variables.bdtIsoNtDz  = cms.string("userFloat('bdtisonontrig')")
process.tpTreeElEl.variables.tkIso    = cms.string("userFloat('tkZZ4L')")
process.tpTreeElEl.variables.ecalIso  = cms.string("userFloat('ecalZZ4L')")
process.tpTreeElEl.variables.hcalIso  = cms.string("userFloat('hcalZZ4L')")
process.tpTreeElEl.variables.rho      = cms.string("userFloat('rhoEl')")
process.tpTreeElEl.variables.pfIsoChHad04       = cms.string("userFloat('electronPFIsoChHad04')")
process.tpTreeElEl.variables.pfIsoNHad04_NoEA   = cms.string("userFloat('electronPFIsoNHad04')")
process.tpTreeElEl.variables.pfIsoPhoton04_NoEA = cms.string("userFloat('electronPFIsoPhoton04')")
process.tpTreeElEl.variables.pfIsoComb04EACorr  = cms.string("userFloat('pfCombIso04EACorr')")
process.tpTreeElEl.variables.pfIsoCombRel04EACorr  = cms.string("userFloat('pfCombIso04EACorr')/pt")
process.tpTreeElEl.variables.sip = cms.string("userFloat('sip')")
process.tpTreeElEl.flags.newID = cms.string("userInt('newID')")
process.tpTreeElEl.flags.bdtID = cms.string("userInt('bdtID')")
process.tpTreeElEl.flags.prlID = cms.string("userInt('prlID')")
process.tpTreeElEl.flags.cicID   = cms.string("test_bit(electronID('cicTight'),0) == 1")
process.tpTreeElEl.flags.prlConv = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")

process.nJetsMuons.objects     = "slimPatJets"
process.nJetsElectrons.objects = "slimPatJets"

process.runMuMu = cms.Path(process.reboosting + process.tnpSimpleSequenceMuMu)
process.runElEl = cms.Path(process.reboosting + process.tnpSimpleSequenceElEl) 

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpLeptonID_HZZ.root"))

