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
process.GlobalTag.globaltag = 'RMMEGlobalTag'
# process.GlobalTag.globaltag = 'START38_V12::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['RMMEFN']
# process.source.fileNames = ['file:/home/mwlebour/isoTest/setup/Hww160.skim.root']
# process.source.fileNames = ['file:/home/mwlebour/isoTest/setup/Hww160.skim.root']
# process.source.fileNames = ['file:latinosYieldSkim.root']
# process.source.fileNames = ['file:/data/mwlebour/MC/Spring11/R414_S1_V06pre2_S2_V02_S3_V00_ID101200/101160/ggToH160toWWto2L2Nu_9_1_Rdc.root']
# from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/hww/WW_39X_ISO_V01/id101160.Flat/*.root') ]
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi")
process.load("RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi")
process.muIsoDepositTk.IOPSet.inputMuonCollection = 'boostedMuons'

process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequencePAT_cff")
process.eleIsoDepositTk.src = 'boostedElectrons'

process.dummy = cms.EDProducer("asdfadf")

ptValues    = [0.00, 0.2, 0.4, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5 ]
outerDRs    = [0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10]
dZValues    = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00, 99.9 ]
tkDeltaRs   = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
tkStrips    = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]




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
    src = cms.InputTag("slimPatJetsTriggerMatch"),
    preselection = cms.string('pt>30 && abs(eta)<5.0'),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("boostedElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
        ),
        mu = cms.PSet(
           src       = cms.InputTag("boostedMuons"),
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

from WWAnalysis.AnalysisStep.electronIDs_cff import ELE_NOCONV,ELE_IP,ELE_ID_LH_90_2011
process.eleLH90 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("boostedElectrons"),
    filter = cms.bool(False),
    cut = cms.string(ELE_NOCONV + "&&" + ELE_IP + "&&" + ELE_ID_LH_90_2011)
)


from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT,MUON_IP_CUT
process.muID = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("boostedMuons"),
    filter = cms.bool(False),
    cut = cms.string(MUON_ID_CUT + "&&" + MUON_IP_CUT)
)

process.goodVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag("goodPrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
    filter = cms.bool(False),   # otherwise it won't filter the events, just produce an empty vertex collection.
)



process.otherStuff = cms.Sequence( process.goodVertices + (process.eleLH90 + process.muID) * process.vetoAbleJets )


#  ______ _           _                    
# |  ____| |         | |                   
# | |__  | | ___  ___| |_ _ __ ___  _ __   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| 
#                                                                 
process.load("WWAnalysis.AnalysisStep.isoAdding_cff")
process.isoAddedElectrons.electronTag = 'boostedElectrons'

process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff")
process.eleIsoFromDepsTk.deposits[0].deltaR = 0.3
process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].deltaR = 0.3
process.eleIsoFromDepsHcalFromTowers.deposits[0].deltaR = 0.3

process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('tkDefault')
process.isoAddedElectrons.deposits[-1].deltaR = 0.3
process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('ecalDefault')
process.isoAddedElectrons.deposits[-1].src = cms.InputTag('myEleIsoDepositEcalFromHits')
process.isoAddedElectrons.deposits[-1].deltaR = 0.3
process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('hcalDefault')
process.isoAddedElectrons.deposits[-1].src = cms.InputTag('myEleIsoDepositHcalFromTowers')
process.isoAddedElectrons.deposits[-1].deltaR = 0.3


process.load("WWAnalysis.AnalysisStep.betterElectronTupleProducer_cff")
process.eleTuple = process.betterElectronTupleProducer.clone()
process.eleTuple.src = "isoAddedElectrons"
process.eleTuple.vtxLabel = "goodVertices"
process.eleTuple.candCounts.nJets = cms.untracked.InputTag("vetoAbleJets")
process.eleTuple.candCounts.nEle80 = cms.untracked.InputTag("eleLH90")
process.eleTuple.candCounts.nMuID  = cms.untracked.InputTag("muID")
# process.eleTuple.candCounts.nGenJet  = cms.untracked.InputTag("genJets25")

process.eleSeq = cms.Sequence( 
    process.isoAddedElectrons * 
    process.eleTuple
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
process.isoAddedMuons.muonTag = 'boostedMuons'


muIsoFromDepsTk =  cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag("muIsoDepositTk"),
        weight = cms.string('1'),
        deltaR = cms.double(0.3),
        vetos = cms.vstring('0.01','Threshold(0.0)'),
        skipDefaultVeto = cms.bool(True),
        label = cms.string('tkDep')
)
process.isoAddedMuons.deposits.append( muIsoFromDepsTk.clone() )

muIsoFromDepsEcal =  cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag('muIsoDepositCalByAssociatorTowers','ecal'),
        weight = cms.string('1'),
        deltaR = cms.double(0.3),
        vetos = cms.vstring(
            'EcalBarrel:ThresholdFromTransverse(0.12)',
            'EcalEndcaps:ThresholdFromTransverse(0.45)',
            'Threshold(0.2)',
            '0.07'
        ),
        skipDefaultVeto = cms.bool(True),
        label = cms.string('ecalDep')
)
process.isoAddedMuons.deposits.append( muIsoFromDepsEcal.clone() )

muIsoFromDepsHcal =  cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag('muIsoDepositCalByAssociatorTowers','hcal'),
        weight = cms.string('1'),
        deltaR = cms.double(0.3),
        vetos = cms.vstring(
            'ThresholdFromTransverse(0.6)',
            'Threshold(0.5)',
            '0.1'
        ),
        skipDefaultVeto = cms.bool(True),
        label = cms.string('hcalDep')
)
process.isoAddedMuons.deposits.append( muIsoFromDepsHcal.clone() )

process.isoAddedMuons.deposits.append( muIsoFromDepsTk.clone() )
process.isoAddedMuons.deposits[-1].label = cms.string('tkDefault')
process.isoAddedMuons.deposits[-1].src = cms.InputTag('muIsoDepositTk')
process.isoAddedMuons.deposits[-1].vetos = []
process.isoAddedMuons.deposits[-1].skipDefaultVeto = False
process.isoAddedMuons.deposits.append( muIsoFromDepsEcal.clone() )
process.isoAddedMuons.deposits[-1].label = cms.string('ecalDefault')
process.isoAddedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','ecal')
process.isoAddedMuons.deposits[-1].vetos = []
process.isoAddedMuons.deposits[-1].skipDefaultVeto = False
process.isoAddedMuons.deposits.append( muIsoFromDepsHcal.clone() )
process.isoAddedMuons.deposits[-1].label = cms.string('hcalDefault')
process.isoAddedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','hcal')
process.isoAddedMuons.deposits[-1].vetos = []
process.isoAddedMuons.deposits[-1].skipDefaultVeto = False

process.load("WWAnalysis.AnalysisStep.betterMuonTupleProducer_cff")
process.muTuple = process.betterMuonTupleProducer.clone()
process.muTuple.src = "isoAddedMuons"
process.muTuple.vtxLabel = "goodVertices"
process.muTuple.candCounts.nJet = cms.untracked.InputTag("vetoAbleJets")
process.muTuple.candCounts.nEle80 = cms.untracked.InputTag("eleLH90")
process.muTuple.candCounts.nMuID  = cms.untracked.InputTag("muID")
# process.muTuple.candCounts.nGenJet  = cms.untracked.InputTag("genJets25")


process.muSeq = cms.Sequence(
    process.isoAddedMuons * 
    process.muTuple
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
process.pfAllNeutral = process.pfAllNeutralHadrons.clone()
process.pfAllCharged = process.pfAllChargedHadrons.clone()
process.pfAllPhotonsNoPU = process.pfAllPhotons.clone()
process.pfAllNeutralNoPU = process.pfAllNeutralHadrons.clone()
process.pfAllChargedNoPU = process.pfAllChargedHadrons.clone()
process.pfAllPhotons.src = "particleFlow"
process.pfAllNeutral.src = "particleFlow"
process.pfAllCharged.src = "particleFlow"
process.pfAllPhotonsNoPU.src = "pfNoPileUp"
process.pfAllNeutralNoPU.src = "pfNoPileUp"
process.pfAllChargedNoPU.src = "pfNoPileUp"

process.isoDepMuWithNeutral = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("boostedMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutral"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)
process.isoDepMuWithPhotons = process.isoDepMuWithNeutral.clone()
process.isoDepMuWithPhotons.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotons")
process.isoDepMuWithCharged = process.isoDepMuWithNeutral.clone()
process.isoDepMuWithCharged.ExtractorPSet.inputCandView = cms.InputTag("pfAllCharged")

process.isoDepMuWithNeutralNoPU = process.isoDepMuWithNeutral.clone()
process.isoDepMuWithNeutralNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllNeutralNoPU")
process.isoDepMuWithPhotonsNoPU = process.isoDepMuWithNeutral.clone()
process.isoDepMuWithPhotonsNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotonsNoPU")
process.isoDepMuWithChargedNoPU = process.isoDepMuWithNeutral.clone()
process.isoDepMuWithChargedNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedNoPU")

process.isoDepElWithNeutral = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("boostedElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutral"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)
process.isoDepElWithPhotons = process.isoDepElWithNeutral.clone()
process.isoDepElWithPhotons.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotons")
process.isoDepElWithCharged = process.isoDepElWithNeutral.clone()
process.isoDepElWithCharged.ExtractorPSet.inputCandView = cms.InputTag("pfAllCharged")

process.isoDepElWithNeutralNoPU = process.isoDepElWithNeutral.clone()
process.isoDepElWithNeutralNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllNeutralNoPU")
process.isoDepElWithPhotonsNoPU = process.isoDepElWithNeutral.clone()
process.isoDepElWithPhotonsNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotonsNoPU")
process.isoDepElWithChargedNoPU = process.isoDepElWithNeutral.clone()
process.isoDepElWithChargedNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedNoPU")

process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)
process.kt6PFJetsForIso.Ghost_EtaMax = cms.double(2.5)
process.kt6PFJetsForIsoNoPU = process.kt6PFJetsForIso.clone( src = "pfNoPileUp" )


process.pfSeq = cms.Sequence( (
        process.pfAllPhotons + 
        process.pfAllNeutral + 
        process.pfAllCharged  
    ) * (
        process.isoDepMuWithNeutral +
        process.isoDepMuWithPhotons + 
        process.isoDepMuWithCharged +
        process.isoDepElWithNeutral +
        process.isoDepElWithPhotons + 
        process.isoDepElWithCharged 
    ) + 
    process.pfPileUp * 
    process.pfNoPileUp * (
        process.kt6PFJetsForIso +
        process.kt6PFJetsForIsoNoPU
    ) * ( 
        process.pfAllChargedNoPU 
    ) * (
        process.isoDepMuWithChargedNoPU +
        process.isoDepElWithChargedNoPU
    )
)

process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('pfPhoton')
process.isoAddedElectrons.deposits[-1].src = cms.InputTag('isoDepElWithPhotons')
process.isoAddedElectrons.deposits[-1].vetos = ['Threshold(0.5)','RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)']
process.isoAddedElectrons.deposits[-1].deltaR = 0.3
process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('pfCharged')
process.isoAddedElectrons.deposits[-1].src = cms.InputTag('isoDepElWithCharged')
process.isoAddedElectrons.deposits[-1].vetos = ['Threshold(0.5)','0.01']
process.isoAddedElectrons.deposits[-1].deltaR = 0.3
process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('pfNeutral')
process.isoAddedElectrons.deposits[-1].src = cms.InputTag('isoDepElWithNeutral')
process.isoAddedElectrons.deposits[-1].vetos = ['Threshold(0.5)']
process.isoAddedElectrons.deposits[-1].deltaR = 0.3

process.isoAddedMuons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedMuons.deposits[-1].label = cms.string('pfPhoton')
process.isoAddedMuons.deposits[-1].src = cms.InputTag('isoDepMuWithPhotons')
process.isoAddedMuons.deposits[-1].vetos = ['Threshold(0.5)']
process.isoAddedMuons.deposits[-1].deltaR = 0.3
process.isoAddedMuons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedMuons.deposits[-1].label = cms.string('pfCharged')
process.isoAddedMuons.deposits[-1].src = cms.InputTag('isoDepMuWithCharged')
process.isoAddedMuons.deposits[-1].vetos = ['Threshold(0.5)','0.01']
process.isoAddedMuons.deposits[-1].deltaR = 0.3
process.isoAddedMuons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedMuons.deposits[-1].label = cms.string('pfNeutral')
process.isoAddedMuons.deposits[-1].src = cms.InputTag('isoDepMuWithNeutral')
process.isoAddedMuons.deposits[-1].vetos = ['Threshold(0.5)']
process.isoAddedMuons.deposits[-1].deltaR = 0.3

process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('pfChargedNoPU')
process.isoAddedElectrons.deposits[-1].src = cms.InputTag('isoDepElWithChargedNoPU')
process.isoAddedElectrons.deposits[-1].vetos = ['Threshold(0.5)','0.01']
process.isoAddedElectrons.deposits[-1].deltaR = 0.3

process.isoAddedMuons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.isoAddedMuons.deposits[-1].label = cms.string('pfChargedNoPU')
process.isoAddedMuons.deposits[-1].src = cms.InputTag('isoDepMuWithChargedNoPU')
process.isoAddedMuons.deposits[-1].vetos = ['Threshold(0.5)','0.01']
process.isoAddedMuons.deposits[-1].deltaR = 0.3


process.eleTuple.variables.append(cms.PSet(tag=cms.untracked.string("tkDefault"),quantity=cms.untracked.string("userFloat('tkDefault')")))
process.eleTuple.variables.append(cms.PSet(tag=cms.untracked.string("ecalDefault"),quantity=cms.untracked.string("userFloat('ecalDefault')")))
process.eleTuple.variables.append(cms.PSet(tag=cms.untracked.string("hcalDefault"),quantity=cms.untracked.string("userFloat('hcalDefault')")))
process.eleTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfPhoton"),quantity=cms.untracked.string("userFloat('pfPhoton')")))
process.eleTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfCharged"),quantity=cms.untracked.string("userFloat('pfCharged')")))
process.eleTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfNeutral"),quantity=cms.untracked.string("userFloat('pfNeutral')")))
process.eleTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfChargedNoPU"),quantity=cms.untracked.string("userFloat('pfChargedNoPU')")))
process.muTuple.variables.append(cms.PSet(tag=cms.untracked.string("tkDefault"),quantity=cms.untracked.string("userFloat('tkDefault')")))
process.muTuple.variables.append(cms.PSet(tag=cms.untracked.string("ecalDefault"),quantity=cms.untracked.string("userFloat('ecalDefault')")))
process.muTuple.variables.append(cms.PSet(tag=cms.untracked.string("hcalDefault"),quantity=cms.untracked.string("userFloat('hcalDefault')")))
process.muTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfPhoton"),quantity=cms.untracked.string("userFloat('pfPhoton')")))
process.muTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfCharged"),quantity=cms.untracked.string("userFloat('pfCharged')")))
process.muTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfNeutral"),quantity=cms.untracked.string("userFloat('pfNeutral')")))
process.muTuple.variables.append(cms.PSet(tag=cms.untracked.string("pfChargedNoPU"),quantity=cms.untracked.string("userFloat('pfChargedNoPU')")))

#  ______ _        _____                 
# |  ____| |      |  __ \                
# | |__  | | ___  | |  | | ___ _ __  ___ 
# |  __| | |/ _ \ | |  | |/ _ \ '_ \/ __|
# | |____| |  __/ | |__| |  __/ |_) \__ \
# |______|_|\___| |_____/ \___| .__/|___/
#                             | |        
#                             |_|        

for deposit in process.isoAddedElectrons.deposits[:]:
    tempDep = deposit.clone()
    tempDep.label = tempDep.label.value()+'LH90'
    for veto in tempDep.vetos[:]:
        tempDep.vetos.append('eleLH90:'+veto)
    process.isoAddedElectrons.deposits.append(tempDep)
    process.eleTuple.variables.append(
        cms.PSet( 
            tag = cms.untracked.string(tempDep.label.value()),        
            quantity = cms.untracked.string("userFloat('"+tempDep.label.value()+"')")
        )
    )


process.eleIsoSeq = cms.Sequence(
    process.eleIsoDepositTk 
#     process.eleIsoDzSeq +
#     process.eleIsoD0Seq 
)

#  __  __         _______             _               _____                 
# |  \/  |       |__   __|           | |             |  __ \                
# | \  / |_   _     | |_ __ __ _  ___| | _____ _ __  | |  | | ___ _ __  ___ 
# | |\/| | | | |    | | '__/ _` |/ __| |/ / _ \ '__| | |  | |/ _ \ '_ \/ __|
# | |  | | |_| |    | | | | (_| | (__|   <  __/ |    | |__| |  __/ |_) \__ \
# |_|  |_|\__,_|    |_|_|  \__,_|\___|_|\_\___|_|    |_____/ \___| .__/|___/
#                                                                 | |        
#                                                                 |_|        

for deposit in process.isoAddedMuons.deposits[:]:
    if deposit.label.value().find('Default') != -1: continue
    tempDep = deposit.clone()
    tempDep.label = tempDep.label.value()+'ID'
    for veto in tempDep.vetos[:]:
        tempDep.vetos.append('muID:'+veto)
    process.isoAddedMuons.deposits.append(tempDep)
    process.muTuple.variables.append(
        cms.PSet( 
            tag = cms.untracked.string(tempDep.label.value()),        
            quantity = cms.untracked.string("userFloat('"+tempDep.label.value()+"')")
        )
    )

process.muIsoSeq = cms.Sequence(
    process.muIsoDepositTk 
)



#  _____ _                              _____  ______   _____ _______ 
# / ____| |                            |  __ \|  ____| |  __ \__   __|
#| |    | |__   __ _ _ __   __ _  ___  | |__) | |__    | |__) | | |   
#| |    | '_ \ / _` | '_ \ / _` |/ _ \ |  ___/|  __|   |  ___/  | |   
#| |____| | | | (_| | | | | (_| |  __/ | |    | |      | |      | |   
# \_____|_| |_|\__,_|_| |_|\__, |\___| |_|    |_|      |_|      |_|   
#                           __/ |                                     
#                          |___/                                      

tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithCharged"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('0.01'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)

for pt in ptValues:
    # Electron 
    newLabel = "chPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithCharged"
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
    # Muon
    newLabel = "chPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedMuons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithCharged"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for pt in ptValues:
    # Electron NoPU
    newLabel = "chNoPUPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithChargedNoPU"
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
    # Muon NoPU
    newLabel = "chNoPUPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedMuons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithChargedNoPU"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))


tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithNeutral"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring(),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)

for pt in ptValues:
    # Electron 
    newLabel = "neuPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithNeutral"
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
    # Muon
    newLabel = "neuPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedMuons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithNeutral"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))


tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithPhotons"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring(),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)

for pt in ptValues:
    # Electron 
    newLabel = "gamPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithPhotons"
    process.isoAddedElectrons.deposits[-1].vetos.append('RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)')
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
    # Muon
    newLabel = "gamPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedMuons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithPhotons"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))


#  _____ _                              _____  ______   _____  _____  
# / ____| |                            |  __ \|  ____| |  __ \|  __ \ 
#| |    | |__   __ _ _ __   __ _  ___  | |__) | |__    | |  | | |__) |
#| |    | '_ \ / _` | '_ \ / _` |/ _ \ |  ___/|  __|   | |  | |  _  / 
#| |____| | | | (_| | | | | (_| |  __/ | |    | |      | |__| | | \ \ 
# \_____|_| |_|\__,_|_| |_|\__, |\___| |_|    |_|      |_____/|_|  \_\
#                           __/ |                                     
#                          |___/                                      

tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithCharged"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('Threshold(0.5)','0.01'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)


for dR in outerDRs:
    # Electron Charged
    newLabel = "chDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithCharged"
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in outerDRs:
    # Muon Charged
    newLabel = "chDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithCharged"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in outerDRs:
    # Electron NoPU Charged
    newLabel = "chNoPUDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithChargedNoPU"
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in outerDRs:
    # Muon NoPU Charged
    newLabel = "chNoPUDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithChargedNoPU"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))

tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithNeutral"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)

for dR in outerDRs:
    # Electron Neutral
    newLabel = "neuDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithNeutral"
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in outerDRs:
    # Muon Neutral
    newLabel = "neuDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithNeutral"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))

tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithNeutral"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)

for dR in outerDRs:
    # Electron Photons
    newLabel = "gamDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithPhotons"
    process.isoAddedElectrons.deposits[-1].vetos.append('RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)')
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in outerDRs:
    # Muon Photons
    newLabel = "gamDR{0:0<4}".format(dR).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dR
    process.isoAddedMuons.deposits[-1].label = newLabel
    process.isoAddedMuons.deposits[-1].src = "isoDepMuWithPhotons"
    process.muTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))



# ____                               _               _____        __      
#|  _ \                             (_)             |_   _|      / _|     
#| |_) |_ __ ___ _ __ ___  _ __ ___  _ _ __   __ _    | |  _ __ | |_ ___  
#|  _ <| '__/ _ \ '_ ` _ \| '_ ` _ \| | '_ \ / _` |   | | | '_ \|  _/ _ \ 
#| |_) | | |  __/ | | | | | | | | | | | | | | (_| |  _| |_| | | | || (_) |
#|____/|_|  \___|_| |_| |_|_| |_| |_|_|_| |_|\__, | |_____|_| |_|_| \___/ 
#                                             __/ |                       
#                                            |___/                        

tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithNeutral"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('Threshold(0.5)','RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)
# tkDeltaRs   = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
# tkStrips    = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]

for dR in tkDeltaRs[5:]:
    # Electron Photons
    newLabel = "gamInDR{0:0<5}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithPhotons"
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dR))
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in tkDeltaRs[5:]:
    # Electron Charged
    newLabel = "chInDR{0:0<5}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithCharged"
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dR))
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in tkDeltaRs[5:]:
    # Electron NoPU Charged
    newLabel = "chNoPUInDR{0:0<5}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithCharged"
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dR))
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))

tempPset = cms.PSet(
    src = cms.InputTag("isoDepElWithNeutral"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('Threshold(0.5)','0.01'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)
# tkDeltaRs   = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
# tkStrips    = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]

for dR in tkDeltaRs:
    # Electron Photons
    newLabel = "gamInDR{0:0<5}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithPhotons"
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dR))
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in tkDeltaRs:
    # Electron Charged
    newLabel = "chInDR{0:0<5}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithCharged"
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dR))
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))
for dR in tkDeltaRs:
    # Electron NoPU Charged
    newLabel = "chNoPUInDR{0:0<5}".format(dR).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].label = newLabel
    process.isoAddedElectrons.deposits[-1].src = "isoDepElWithCharged"
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dR))
    process.eleTuple.variables.append( cms.PSet( tag = cms.untracked.string(newLabel), quantity = cms.untracked.string("userFloat('"+newLabel+"')")))





#__      __                _   _            _____ ______
#\ \    / /               | | | |          |  __ \___  /
# \ \  / /_ _ _ __ _   _  | |_| |__   ___  | |  | | / / 
#  \ \/ / _` | '__| | | | | __| '_ \ / _ \ | |  | |/ /  
#   \  / (_| | |  | |_| | | |_| | | |  __/ | |__| / /__ 
#    \/ \__,_|_|   \__, |  \__|_| |_|\___| |_____/_____|
#                   __/ |                               
#                  |___/                

# Run Dzs
process.pfDZs = cms.Sequence( process.dummy )
for dz in dZValues:
    #name of this iso deposit
    isoDepName = "elChDZ{0:0<4}".format(dz).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.isoDepElWithCharged.clone())
    #change the dz
    getattr(process,isoDepName).ExtractorPSet.Diff_z = dz
    #add the VM to the boostedElectron
    process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
    process.isoAddedElectrons.deposits[-1].src = isoDepName
    process.isoAddedElectrons.deposits[-1].label = cms.string(isoDepName)
    process.isoAddedElectrons.deposits[-1].deltaR = 0.3
    process.isoAddedElectrons.deposits[-1].vetos = ['Threshold(0.5)','0.01']
    #append to the isoSeq
    process.pfDZs.replace(process.dummy,getattr(process,isoDepName)+process.dummy)
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(isoDepName),              
            quantity = cms.untracked.string("userFloat('"+isoDepName+"')"),
        )
    )
    #name of this iso deposit
    isoDepName = "muChDZ{0:0<4}".format(dz).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.isoDepMuWithCharged.clone())
    #change the dz
    getattr(process,isoDepName).ExtractorPSet.Diff_z = dz
    #add the VM to the boostedElectron
    process.isoAddedMuons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
    process.isoAddedMuons.deposits[-1].src = isoDepName
    process.isoAddedMuons.deposits[-1].label = cms.string(isoDepName)
    process.isoAddedMuons.deposits[-1].deltaR = 0.3
    process.isoAddedMuons.deposits[-1].vetos = ['Threshold(0.5)','0.01']
    #append to the isoSeq
    process.pfDZs.replace(process.dummy,getattr(process,isoDepName)+process.dummy)
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(isoDepName),              
            quantity = cms.untracked.string("userFloat('"+isoDepName+"')"),
        )
    )

process.pfDZs.remove(process.dummy)




#  ______           _   _____      _   _        _____ _          __  __ 
# |  ____|         | | |  __ \    | | | |      / ____| |        / _|/ _|
# | |__   _ __   __| | | |__) |_ _| |_| |__   | (___ | |_ _   _| |_| |_ 
# |  __| | '_ \ / _` | |  ___/ _` | __| '_ \   \___ \| __| | | |  _|  _|
# | |____| | | | (_| | | |  | (_| | |_| | | |  ____) | |_| |_| | | | |  
# |______|_| |_|\__,_| |_|   \__,_|\__|_| |_| |_____/ \__|\__,_|_| |_|  
#                                                                       


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

process.TFileService = cms.Service("TFileService",fileName = cms.string("RMMEFN"))

isBkg = RMMEBKG
# isBkg = False
if isBkg:
    process.elePath = cms.Path(~process.elFromW * ~process.elFromTau * process.otherStuff * process.eleIsoSeq * process.pfSeq * process.pfDZs * process.eleSeq)
    process.muPath  = cms.Path(~process.muFromW * ~process.muFromTau * process.otherStuff * process.muIsoSeq  * process.pfSeq * process.pfDZs * process.muSeq)
else:
    process.elePath = cms.Path(process.otherStuff * process.eleIsoSeq * process.pfSeq * process.pfDZs * process.eleSeq)
    process.muPath  = cms.Path(process.otherStuff * process.muIsoSeq  * process.pfSeq * process.pfDZs * process.muSeq)

