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
process.source.fileNames = ['file:latinosYieldSkim.root']
# from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/hww/WW_39X_ISO_V01/id101160.Flat/*.root') ]
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi")
process.load("RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi")
process.muIsoDepositTk.IOPSet.inputMuonCollection = 'boostedMuons'

process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequencePAT_cff")
process.eleIsoDepositTk.src = 'boostedElectrons'

process.dummy = cms.EDProducer("asdfadf")

dZValues    = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00, 99.9 ]
d0Values    = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00, 99.9 ]
ptValues    = [0.00, 0.2, 0.4, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5 ]
tkStrips    = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
tkDeltaRs   = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
outerDRs    = [0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10]
eneThresh   = [0.02, 0.04, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.14, 0.15, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40 ] 
etThresh    = [0.02, 0.04, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.14, 0.15, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40 ] 
crystStrips = [0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00]
crystCones  = [0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00]


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
    preselection = cms.string('pt>25 && abs(eta)<5.0'),
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

# process.genJets25 = cms.EDFilter("CandViewSelector",
#     src = cms.InputTag("highPtGenJets"),
#     filter = cms.bool(False),
#     cut = cms.string("pt > 25. && abs(eta) < 5.")
# )

process.ele80 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("boostedElectrons"),
    filter = cms.bool(False),
    cut = cms.string("electronID('eidMedium')==1 || electronID('eidMedium')==5 || electronID('eidMedium')==7 || electronID('eidMedium')==3")
)

process.ele95 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("boostedElectrons"),
    filter = cms.bool(False),
    cut = cms.string("electronID('eidLoose')==1 || electronID('eidLoose')==5 || electronID('eidLoose')==7 || electronID('eidLoose')==3")
)

from WWAnalysis.AnalysisStep.betterMuonTupleProducer_cff import isGood
from WWAnalysis.AnalysisStep.betterMuonTupleProducer_cff import isPrompt
process.muID = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("boostedMuons"),
    filter = cms.bool(False),
    cut = cms.string(isGood+'&&'+isPrompt)
)

process.goodVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag("goodPrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
    filter = cms.bool(False),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.load('RecoJets.JetProducers.kt4PFJets_cfi') 
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

process.otherStuff = cms.Sequence( process.kt6PFJets + process.goodVertices + (process.ele80 + process.ele95 + process.muID) * process.vetoAbleJets )
# process.otherStuff = cms.Sequence( process.kt6PFJets + process.genJets25 + process.goodVertices + (process.ele80 + process.ele95 + process.muID) * process.vetoAbleJets )


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
process.isoAddedElectrons.deposits[-1].deltaR = 0.3
process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.isoAddedElectrons.deposits[-1].label = cms.string('hcalDefault')
process.isoAddedElectrons.deposits[-1].deltaR = 0.3


process.load("WWAnalysis.AnalysisStep.betterElectronTupleProducer_cff")
process.eleTuple = process.betterElectronTupleProducer.clone()
process.eleTuple.vtxLabel = "goodVertices"
process.eleTuple.candCounts.nJets = cms.untracked.InputTag("vetoAbleJets")
process.eleTuple.candCounts.nEle80 = cms.untracked.InputTag("ele80")
process.eleTuple.candCounts.nEle95 = cms.untracked.InputTag("ele95")
process.eleTuple.candCounts.nMuID  = cms.untracked.InputTag("muID")
# process.eleTuple.candCounts.nGenJet  = cms.untracked.InputTag("genJets25")

process.eleSeq = cms.Sequence( 
    process.isoAddedElectrons * 
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

for deposit in process.isoAddedElectrons.deposits[:]:
    tempDep = deposit.clone()
    tempDep.label = tempDep.label.value()+'95'
    for veto in tempDep.vetos[:]:
        tempDep.vetos.append('ele95:'+veto)
    process.isoAddedElectrons.deposits.append(tempDep)
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
    process.isoAddedElectrons.deposits.append(tempDep)
    process.eleTuple.variables.append(
        cms.PSet( 
            tag = cms.untracked.string(tempDep.label.value()),        
            quantity = cms.untracked.string("userFloat('"+tempDep.label.value()+"')")
        )
    )

tempPset = cms.PSet(
    src = cms.InputTag("eleIsoDepositTk"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)


for pt in ptValues:
    newLabel = "tk000Strip{0:0<5}".format(pt).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )

tempPset.vetos = cms.vstring('Threshold(0.7)')
for strip in tkStrips:
    newLabel = "tk000Strip{0:0<5}".format(strip).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('RectangularEtaPhiVeto(-'+str(strip)+','+str(strip)+',-0.6,0.6)')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
tempPset.vetos = cms.vstring('Threshold(0.7)','0.015')
for strip in tkStrips:
    newLabel = "tk015Strip{0:0<5}".format(strip).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('RectangularEtaPhiVeto(-'+str(strip)+','+str(strip)+',-0.6,0.6)')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
tempPset.vetos = cms.vstring('Threshold(0.7)')
for dr in tkDeltaRs:
    newLabel = "tk000DR{0:0<5}".format(dr).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dr))
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
tempPset.vetos = cms.vstring('Threshold(0.7)','RectangularEtaPhiVeto(-0.015,0.015,-0.6,0.6)')
for dr in tkDeltaRs:
    newLabel = "tk015DR{0:0<5}".format(dr).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append(str(dr))
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
for dr in outerDRs:
    newLabel = "tkOuter{0:0<4}".format(dr).replace('.','')
    process.isoAddedElectrons.deposits.append(process.eleIsoFromDepsTk.deposits[0].clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dr
    process.isoAddedElectrons.deposits[-1].label = cms.string(newLabel)
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )

for dr in outerDRs:
    newLabel = "ecalOuter{0:0<4}".format(dr).replace('.','')
    process.isoAddedElectrons.deposits.append(process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dr
    process.isoAddedElectrons.deposits[-1].label = cms.string(newLabel)
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )

for dr in outerDRs:
    newLabel = "hcalOuter{0:0<4}".format(dr).replace('.','')
    process.isoAddedElectrons.deposits.append(process.eleIsoFromDepsHcalFromTowers.deposits[0].clone())
    process.isoAddedElectrons.deposits[-1].deltaR = dr
    process.isoAddedElectrons.deposits[-1].label = cms.string(newLabel)
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    

tempPset = process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].clone()
tempPset.label = cms.string('')
tempPset.vetos = [
    'NumCrystalVeto(3.0)',
    'NumCrystalEtaPhiVeto(1.5,9999.0)',
]
for energy in eneThresh:
    newLabel = "ecalEne{0:0<4}".format(energy).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('AbsThresholdFromTransverse('+str(energy)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
tempPset.vetos = [
    'NumCrystalVeto(3.0)',
    'NumCrystalEtaPhiVeto(1.5,9999.0)',
]
for et in etThresh:
    newLabel = "ecalEt{0:0<4}".format(et).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('AbsThreshold('+str(et)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
tempPset.vetos = [
    'NumCrystalVeto(3.0)',
    'EcalBarrel:AbsThresholdFromTransverse(0.08)',
    'EcalEndcaps:AbsThreshold(0.100)'
]
for strip in crystStrips:
    newLabel = "ecalStrip{0:0<4}".format(strip).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('NumCrystalEtaPhiVeto('+str(strip)+',9999.0)')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
tempPset.vetos = [
    'NumCrystalEtaPhiVeto(1.5,9999.0)',
    'EcalBarrel:AbsThresholdFromTransverse(0.08)',
    'EcalEndcaps:AbsThreshold(0.100)'
]
for dr in crystCones:
    newLabel = "ecalDR{0:0<3}".format(dr).replace('.','')
    process.isoAddedElectrons.deposits.append(tempPset.clone())
    process.isoAddedElectrons.deposits[-1].vetos.append('NumCrystalVeto('+str(dr)+')')
    process.isoAddedElectrons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    

#Run Dzs
process.eleIsoDzSeq = cms.Sequence( process.dummy )
for dz in dZValues:
    #name of this iso deposit
    isoDepName = "elTkDZ{0:0<4}".format(dz).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.eleIsoDepositTk.clone())
    #change the dz
    getattr(process,isoDepName).ExtractorPSet.Diff_z = dz
    #add the VM to the boostedElectron
    process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
    process.isoAddedElectrons.deposits[-1].src = isoDepName
    process.isoAddedElectrons.deposits[-1].label = cms.string(isoDepName)
    process.isoAddedElectrons.deposits[-1].deltaR = 0.3
    #append to the isoSeq
    process.eleIsoDzSeq.replace(process.dummy,getattr(process,isoDepName)+process.dummy)
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(isoDepName),              
            quantity = cms.untracked.string("userFloat('"+isoDepName+"')"),
        )
    )
process.eleIsoDzSeq.remove(process.dummy)

#Run D0s
process.eleIsoD0Seq = cms.Sequence( process.dummy )
for d0 in d0Values:
    #name of this iso deposit
    isoDepName = "elTkD0{0:0<4}".format(d0).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.eleIsoDepositTk.clone())
    #change the d0
    getattr(process,isoDepName).ExtractorPSet.Diff_r = d0
    #add the VM to the boostedElectron
    process.isoAddedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
    process.isoAddedElectrons.deposits[-1].src = isoDepName
    process.isoAddedElectrons.deposits[-1].label = cms.string(isoDepName)
    process.isoAddedElectrons.deposits[-1].deltaR = 0.3
    #append to the isoSeq
    process.eleIsoD0Seq.replace(process.dummy,getattr(process,isoDepName)+process.dummy)
    #add the variable to the tree
    process.eleTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(isoDepName),              
            quantity = cms.untracked.string("userFloat('"+isoDepName+"')"),
        )
    )
process.eleIsoD0Seq.remove(process.dummy)


process.eleIsoSeq = cms.Sequence(
    process.eleIsoDepositTk +
    process.eleIsoDzSeq +
    process.eleIsoD0Seq 
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
process.muTuple.vtxLabel = "goodVertices"
process.muTuple.candCounts.nJet = cms.untracked.InputTag("vetoAbleJets")
process.muTuple.candCounts.nEle80 = cms.untracked.InputTag("ele80")
process.muTuple.candCounts.nEle95 = cms.untracked.InputTag("ele95")
process.muTuple.candCounts.nMuID  = cms.untracked.InputTag("muID")
# process.muTuple.candCounts.nGenJet  = cms.untracked.InputTag("genJets25")


process.muSeq = cms.Sequence(
    process.isoAddedMuons * 
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

for deposit in process.isoAddedMuons.deposits[:]:
    if deposit.label.value().find('Dep') != -1:
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

tempPset = cms.PSet(
    src = cms.InputTag("muIsoDepositTk"),
    deltaR = cms.double(0.3),
    weight = cms.string('1'),
    label = cms.string('tkDefault95'),
    vetos = cms.vstring('0.01'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
)


for pt in ptValues:
    newLabel = "tkPt{0:0<4}".format(pt).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedMuons.deposits[-1].vetos.append('Threshold('+str(pt)+')')
    process.isoAddedMuons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )

tempPset.vetos = cms.vstring('Threshold(0.0)')
for dr in tkDeltaRs:
    newLabel = "tkDR{0:0<5}".format(dr).replace('.','')
    process.isoAddedMuons.deposits.append(tempPset.clone())
    process.isoAddedMuons.deposits[-1].vetos.append(str(dr))
    process.isoAddedMuons.deposits[-1].label = newLabel
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    
for dr in outerDRs:
    newLabel = "tkOuter{0:0<4}".format(dr).replace('.','')
    process.isoAddedMuons.deposits.append(muIsoFromDepsTk.clone())
    process.isoAddedMuons.deposits[-1].deltaR = dr
    process.isoAddedMuons.deposits[-1].label = cms.string(newLabel)
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )

for dr in outerDRs:
    newLabel = "ecalOuter{0:0<4}".format(dr).replace('.','')
    process.isoAddedMuons.deposits.append(muIsoFromDepsEcal.clone())
    process.isoAddedMuons.deposits[-1].deltaR = dr
    process.isoAddedMuons.deposits[-1].label = cms.string(newLabel)
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )

for dr in outerDRs:
    newLabel = "hcalOuter{0:0<4}".format(dr).replace('.','')
    process.isoAddedMuons.deposits.append(muIsoFromDepsHcal.clone())
    process.isoAddedMuons.deposits[-1].deltaR = dr
    process.isoAddedMuons.deposits[-1].label = cms.string(newLabel)
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(newLabel),              
            quantity = cms.untracked.string("userFloat('"+newLabel+"')"),
        )
    )
    

    
#Run Dzs
process.muIsoDzSeq = cms.Sequence( process.dummy )
for dz in dZValues:
    #name of this iso deposit
    isoDepName = "muTkDZ{0:0<4}".format(dz).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.muIsoDepositTk.clone())
    #change the dz
    getattr(process,isoDepName).ExtractorPSet.Diff_z = dz
    #add the VM to the boostedElectron
    process.isoAddedMuons.deposits.append( muIsoFromDepsTk.clone() )
    process.isoAddedMuons.deposits[-1].src = isoDepName
    process.isoAddedMuons.deposits[-1].label = cms.string(isoDepName)
    process.isoAddedMuons.deposits[-1].deltaR = 0.3
    #append to the isoSeq
    process.muIsoDzSeq.replace(process.dummy,getattr(process,isoDepName)+process.dummy)
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(isoDepName),              
            quantity = cms.untracked.string("userFloat('"+isoDepName+"')"),
        )
    )
process.muIsoDzSeq.remove(process.dummy)

#Run D0s
process.muIsoD0Seq = cms.Sequence( process.dummy )
for d0 in d0Values:
    #name of this iso deposit
    isoDepName = "muTkD0{0:0<5}".format(d0).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.muIsoDepositTk.clone())
    #change the d0
    getattr(process,isoDepName).ExtractorPSet.Diff_r = d0
    #add the VM to the boostedElectron
    process.isoAddedMuons.deposits.append( muIsoFromDepsTk.clone() )
    process.isoAddedMuons.deposits[-1].src = isoDepName
    process.isoAddedMuons.deposits[-1].label = cms.string(isoDepName)
    process.isoAddedMuons.deposits[-1].deltaR = 0.3
    #append to the isoSeq
    process.muIsoD0Seq.replace(process.dummy,getattr(process,isoDepName)+process.dummy)
    #add the variable to the tree
    process.muTuple.variables.append( 
        cms.PSet( 
            tag = cms.untracked.string(isoDepName),              
            quantity = cms.untracked.string("userFloat('"+isoDepName+"')"),
        )
    )
process.muIsoD0Seq.remove(process.dummy)


process.muIsoSeq = cms.Sequence(
    process.muIsoDepositTk +
    process.muIsoDzSeq +
    process.muIsoD0Seq 
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
    src = cms.InputTag("boostedMuons"),
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

process.TFileService = cms.Service("TFileService",fileName = cms.string("isoIpStudy.2.root"))

isBkg = False
if isBkg:
    process.elePath = cms.Path(~process.elFromW * ~process.elFromTau * process.otherStuff * process.eleIsoSeq * process.eleSeq)
    process.muPath  = cms.Path(~process.muFromW * ~process.muFromTau * process.otherStuff * process.muIsoSeq  * process.muSeq)
else:
    process.elePath = cms.Path(process.otherStuff * process.eleIsoSeq * process.eleSeq)
    process.muPath  = cms.Path(process.otherStuff * process.muIsoSeq  * process.muSeq)

