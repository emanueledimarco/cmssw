import FWCore.ParameterSet.Config as cms

process = cms.Process("Iso")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")
process.GlobalTag.globaltag = 'START38_V12::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
# process.source.fileNames = ['file:/home/mwlebour/isoTest/setup/Hww160.skim.root']
# process.source.fileNames = ['file:/home/mwlebour/isoTest/setup/Hww160.skim.root']
from glob import glob
process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/hww/WW_39X_ISO_V01/id101160.Flat/*.root') ]

process.load("RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi")
process.load("RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi")
process.muIsoDepositTk.IOPSet.inputMuonCollection = 'goodMuons'

process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequencePAT_cff")
process.eleIsoDepositTk.src = 'goodElectrons'

process.dummy = cms.EDProducer("asdfadf")

dZValues   = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00, 99.9 ]
d0Values   = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00, 99.9 ]
tkStrips   = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025]
tkDeltaRs  = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025]
outerDRs   = [0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]
eneThresh  = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
              0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 ] 
etThresh   = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
              0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 ] 


#   ____  _   _                  _____ _          __  __ 
#  / __ \| | | |                / ____| |        / _|/ _|
# | |  | | |_| |__   ___ _ __  | (___ | |_ _   _| |_| |_ 
# | |  | | __| '_ \ / _ \ '__|  \___ \| __| | | |  _|  _|
# | |__| | |_| | | |  __/ |     ____) | |_| |_| | | | |  
#  \____/ \__|_| |_|\___|_|    |_____/ \__|\__,_|_| |_|  
#                                                        

# process.vetoAbleJets = cms.EDFilter("CandViewSelector",
#     src = cms.InputTag("cleanPatJets"),
#     filter = cms.bool(False),
#     cut = cms.string("pt > 25 && abs(eta) < 5.0")
# )

process.vetoAbleJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("cleanPatJetsPF"),
    preselection = cms.string('pt>25 && abs(eta)<5.0'),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("goodElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
        ),
        mu = cms.PSet(
           src       = cms.InputTag("goodMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
        ),
    ),
    finalCut = cms.string(''),
)

process.genJets25 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("highPtGenJets"),
    filter = cms.bool(False),
    cut = cms.string("pt > 25. && abs(eta) < 5.")
)

process.ele80 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("goodElectrons"),
    filter = cms.bool(False),
    cut = cms.string("electronID('eidVBTFRel80')==1 || electronID('eidVBTFRel80')==5 || electronID('eidVBTFRel80')==7 || electronID('eidVBTFRel80')==3")
)

process.ele95 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("goodElectrons"),
    filter = cms.bool(False),
    cut = cms.string("electronID('eidVBTFRel95')==1 || electronID('eidVBTFRel95')==5 || electronID('eidVBTFRel95')==7 || electronID('eidVBTFRel95')==3")
)

from WWAnalysis.AnalysisStep.betterMuonTupleProducer_cff import isGood
from WWAnalysis.AnalysisStep.betterMuonTupleProducer_cff import isPrompt
process.muID = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("goodMuons"),
    filter = cms.bool(False),
    cut = cms.string(isGood+'&&'+isPrompt)
)

process.goodVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag("goodPrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
    filter = cms.bool(False),   # otherwise it won't filter the events, just produce an empty vertex collection.
)


process.otherStuff = cms.Sequence( process.genJets25 + process.goodVertices + (process.ele80 + process.ele95 + process.muID) * process.vetoAbleJets )


#  ______ _           _                    
# |  ____| |         | |                   
# | |__  | | ___  ___| |_ _ __ ___  _ __   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| 
#                                                                 
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.boostedElectrons.electronTag = 'goodElectrons'

process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff")
process.boostedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string('tkDefault')
process.boostedElectrons.deposits[-1].deltaR = 0.3
process.boostedElectrons.deposits.append( process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string('ecalDefault')
process.boostedElectrons.deposits[-1].deltaR = 0.3
process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string('hcalDefault')
process.boostedElectrons.deposits[-1].deltaR = 0.3


process.load("WWAnalysis.AnalysisStep.betterElectronTupleProducer_cff")
process.eleTuple = process.betterElectronTupleProducer.clone()
process.eleTuple.vtxLabel = "goodVertices"
process.eleTuple.candCounts.nJets = cms.untracked.InputTag("vetoAbleJets")
process.eleTuple.candCounts.nEle80 = cms.untracked.InputTag("ele80")
process.eleTuple.candCounts.nEle95 = cms.untracked.InputTag("ele95")
process.eleTuple.candCounts.nMuID  = cms.untracked.InputTag("muID")
process.eleTuple.candCounts.nGenJet  = cms.untracked.InputTag("genJets25")

process.eleSeq = cms.Sequence( 
    process.boostedElectrons * 
    process.eleTuple
)

#  ______ _        _____                 
# |  ____| |      |  __ \                
# | |__  | | ___  | |  | | ___ _ __  ___ 
# |  __| | |/ _ \ | |  | |/ _ \ '_ \/ __|
# | |____| |  __/ | |__| |  __/ |_) \__ \
# |______|_|\___| |_____/ \___| .__/|___/
#                             | |        
#                             |_|        

process.eleIsoSeq = cms.Sequence(
    process.eleIsoDepositTk 
)


for deposit in process.boostedElectrons.deposits[:]:
    tempDep = deposit.clone()
    tempDep.label = tempDep.label.value()+'95'
    for veto in tempDep.vetos[:]:
        tempDep.vetos.append('ele95:'+veto)
    process.boostedElectrons.deposits.append(tempDep)
    process.eleTuple.variables.append(
        cms.PSet( 
            tag = cms.untracked.string(tempDep.label.value()),        
            quantity = cms.untracked.string("userFloat('"+tempDep.label.value()+"')")
        )
    )
    tempDep = deposit.clone()
    tempDep.label = tempDep.label.value()+'80'
    for veto in tempDep.vetos[:]:
        tempDep.vetos.append('ele80:'+veto)
    process.boostedElectrons.deposits.append(tempDep)
    process.eleTuple.variables.append(
        cms.PSet( 
            tag = cms.untracked.string(tempDep.label.value()),        
            quantity = cms.untracked.string("userFloat('"+tempDep.label.value()+"')")
        )
    )


#  __  __                   
# |  \/  |                  
# | \  / |_   _  ___  _ __  
# | |\/| | | | |/ _ \| '_ \ 
# | |  | | |_| | (_) | | | |
# |_|  |_|\__,_|\___/|_| |_|
#                           
#                                                     

process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")
process.boostedMuons.muonTag = 'goodMuons'


tempPSet =  cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag("muIsoDepositTk"),
        weight = cms.string('1'),
        deltaR = cms.double(0.3),
        vetos = cms.vstring('0.01','Threshold(0.0)'),
        skipDefaultVeto = cms.bool(True),
        label = cms.string('isoTk')
)


process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('tkDep')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositTk')
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('ecalDep')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','ecal')
process.boostedMuons.deposits[-1].vetos = [
    'EcalBarrel:ThresholdFromTransverse(0.12)',
    'EcalEndcaps:ThresholdFromTransverse(0.45)',
    'Threshold(0.2)',
    '0.07'
]
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('hcalDep')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','hcal')
process.boostedMuons.deposits[-1].vetos = [
    'ThresholdFromTransverse(0.6)',
    'Threshold(0.5)',
    '0.1'
]

process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('tkDefault')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositTk')
process.boostedMuons.deposits[-1].vetos = []
process.boostedMuons.deposits[-1].skipDefaultVeto = False
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('ecalDefault')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','ecal')
process.boostedMuons.deposits[-1].vetos = []
process.boostedMuons.deposits[-1].skipDefaultVeto = False
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('hcalDefault')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','hcal')
process.boostedMuons.deposits[-1].vetos = []
process.boostedMuons.deposits[-1].skipDefaultVeto = False

process.load("WWAnalysis.AnalysisStep.betterMuonTupleProducer_cff")
process.muTuple = process.betterMuonTupleProducer.clone()
process.muTuple.vtxLabel = "goodVertices"
process.muTuple.candCounts.nJet = cms.untracked.InputTag("vetoAbleJets")
process.muTuple.candCounts.nEle80 = cms.untracked.InputTag("ele80")
process.muTuple.candCounts.nEle95 = cms.untracked.InputTag("ele95")
process.muTuple.candCounts.nMuID  = cms.untracked.InputTag("muID")
process.muTuple.candCounts.nGenJet  = cms.untracked.InputTag("genJets25")


process.muSeq = cms.Sequence(
    process.boostedMuons * 
    process.muTuple
)


#  __  __         _______             _               _____                 
# |  \/  |       |__   __|           | |             |  __ \                
# | \  / |_   _     | |_ __ __ _  ___| | _____ _ __  | |  | | ___ _ __  ___ 
# | |\/| | | | |    | | '__/ _` |/ __| |/ / _ \ '__| | |  | |/ _ \ '_ \/ __|
# | |  | | |_| |    | | | | (_| | (__|   <  __/ |    | |__| |  __/ |_) \__ \
# |_|  |_|\__,_|    |_|_|  \__,_|\___|_|\_\___|_|    |_____/ \___| .__/|___/
#                                                                 | |        
#                                                                 |_|        

for deposit in process.boostedMuons.deposits[:]:
    if deposit.label.value().find('Dep') != -1:
        tempDep = deposit.clone()
        tempDep.label = tempDep.label.value()+'ID'
        for veto in tempDep.vetos[:]:
            tempDep.vetos.append('muID:'+veto)
        process.boostedMuons.deposits.append(tempDep)
        process.muTuple.variables.append(
            cms.PSet( 
                tag = cms.untracked.string(tempDep.label.value()),        
                quantity = cms.untracked.string("userFloat('"+tempDep.label.value()+"')")
            )
        )


process.muIsoSeq = cms.Sequence(
    process.muIsoDepositTk
)



#
#  _____  ______    _____ _          __  __ 
# |  __ \|  ____|  / ____| |        / _|/ _|
# | |__) | |__    | (___ | |_ _   _| |_| |_ 
# |  ___/|  __|    \___ \| __| | | |  _|  _|
# | |    | |       ____) | |_| |_| | | | |  
# |_|    |_|      |_____/ \__|\__,_|_| |_|  
#                                           
#                                                                                     

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.isoDepMuonWithNeutralHadrons = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("goodMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)

process.pfAllPhotonsNoPU = process.pfAllPhotons.clone()
process.pfAllNeutralHadronsNoPU = process.pfAllNeutralHadrons.clone()
process.pfAllChargedHadronsNoPU = process.pfAllChargedHadrons.clone()
process.pfAllPhotons.src = "particleFlow"
process.pfAllNeutralHadrons.src = "particleFlow"
process.pfAllChargedHadrons.src = "particleFlow"

process.isoDepMuonWithPhotons = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithPhotons.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotons")
process.isoDepMuonWithChargedPFlow = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithChargedPFlow.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedHadrons")

process.isoDepMuonWithNeutralHadronsNoPU = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithNeutralHadronsNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllNeutralHadronsNoPU")
process.isoDepMuonWithPhotonsNoPU = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithPhotonsNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotonsNoPU")
process.isoDepMuonWithChargedPFlowNoPU = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithChargedPFlowNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedHadronsNoPU")

process.pfSeq = cms.Sequence( (
        process.pfAllPhotons + 
        process.pfAllNeutralHadrons + 
        process.pfAllChargedHadrons  
    ) * (
        process.isoDepMuonWithNeutralHadrons +
        process.isoDepMuonWithPhotons + 
        process.isoDepMuonWithChargedPFlow
    ) + 
    process.pfPileUp * 
    process.pfNoPileUp * ( 
        process.pfAllPhotonsNoPU + 
        process.pfAllNeutralHadronsNoPU + 
        process.pfAllChargedHadronsNoPU 
    ) * (
        process.isoDepMuonWithNeutralHadronsNoPU +
        process.isoDepMuonWithPhotonsNoPU + 
        process.isoDepMuonWithChargedPFlowNoPU
    )
)


process.elFromW = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==11 && abs(mother.mother.pdgId)==24")
)

process.elFromTau = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==11 && abs(mother.pdgId)==15")
)

process.muFromW = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==13 && abs(mother.mother.pdgId)==24")
)

process.muFromTau = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("prunedGen"),
    filter = cms.bool(True),
    cut = cms.string("abs(pdgId)==13 && abs(mother.pdgId)==15")
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("isoIpStudy.root"))

isBkg = False
if isBkg:
    process.elePath = cms.Path(~process.elFromW * ~process.elFromTau * process.eleIsoSeq * process.eleSeq)
    process.muPath  = cms.Path(~process.muFromW * ~process.muFromTau * process.muIsoSeq  * process.muSeq)
else:
    process.elePath = cms.Path(process.otherStuff * process.eleIsoSeq * process.eleSeq)
    process.muPath  = cms.Path(process.otherStuff * process.muIsoSeq  * process.muSeq)

