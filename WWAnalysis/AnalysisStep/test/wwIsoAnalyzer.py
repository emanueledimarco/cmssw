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
# process.source.fileNames = ['file:../../isoTest/WJetsMad.skim.root']
process.source.fileNames = ['file:wjetsSkim.root']
# from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('*.root') ]

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

weight15 = [3.05902e-07, 4.58853e-06, 3.4414e-05, 0.00017207, 0.000645263, 
            0.00193579, 0.00483947, 0.0103703, 0.0194443, 0.0324072, 0.0486108, 
            0.0662874, 0.0828592, 0.0956068, 0.102436, 0.102436, 0.0960336, 0.0847356, 
            0.070613, 0.0557471, 0.0418103, 0.0298645, 0.0203622, 0.0132797, 
            0.00829979, 0.00497988, 0.00287301, 0.00159611, 0.000855061, 0.000442273]


weight10 = [4.53999e-05, 0.000453999, 0.00227, 0.00756665, 0.0189166, 0.0378333, 
            0.0630555, 0.0900792, 0.112599, 0.12511, 0.12511, 0.113736, 0.0947803, 
            0.0729079, 0.0520771, 0.0347181, 0.0216988, 0.012764, 0.00709111, 
            0.00373216, 0.00186608, 0.00088861, 0.000403914, 0.000175615, 7.31728e-05, 
            2.92691e-05, 1.12573e-05, 4.16939e-06, 1.48907e-06, 5.13472e-07]

weight05 = [0.00673795, 0.0336897, 0.0842243, 0.140374, 0.175467, 0.175467, 0.146223, 
            0.104445, 0.065278, 0.0362656, 0.0181328, 0.00824218, 0.00343424, 0.00132086, 
            0.000471736, 0.000157245, 4.91392e-05, 1.44527e-05, 4.01464e-06, 1.05648e-06, 
            2.64121e-07, 6.2886e-08, 1.42923e-08, 3.10701e-09, 6.47295e-10, 1.29459e-10, 
            2.48959e-11, 4.61036e-12, 8.23279e-13, 1.41945e-13 ]




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
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string('hcalDefault')
process.boostedElectrons.deposits[-1].deltaR = 0.3

process.ele80 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("boostedElectrons"),
    filter = cms.bool(False),
    cut = cms.string("electronID('eidVBTFRel80')==1 || electronID('eidVBTFRel80')==5 || electronID('eidVBTFRel80')==7 || electronID('eidVBTFRel80')==3")
)

process.ele95 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("boostedElectrons"),
    filter = cms.bool(False),
    cut = cms.string("electronID('eidVBTFRel95')==1 || electronID('eidVBTFRel95')==5 || electronID('eidVBTFRel95')==7 || electronID('eidVBTFRel95')==3")
)

process.eleTuple = cms.EDAnalyzer("BetterElectronTupleProducer", 
    src = cms.untracked.InputTag("boostedElectrons"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("eid80"),              quantity = cms.untracked.string("electronID('eidVBTFRel80')")),
        cms.PSet( tag = cms.untracked.string("eid95"),              quantity = cms.untracked.string("electronID('eidVBTFRel95')")),
        cms.PSet( tag = cms.untracked.string("liklihoodID"),        quantity = cms.untracked.string("electronID('liklihoodID')")),
        cms.PSet( tag = cms.untracked.string("pt"),                 quantity = cms.untracked.string("pt")),
        cms.PSet( tag = cms.untracked.string("eta"),                quantity = cms.untracked.string("eta")),
        cms.PSet( tag = cms.untracked.string("phi"),                quantity = cms.untracked.string("phi")),
        cms.PSet( tag = cms.untracked.string("eb"),                 quantity = cms.untracked.string("isEB")),
        cms.PSet( tag = cms.untracked.string("ee"),                 quantity = cms.untracked.string("isEE")),
        cms.PSet( tag = cms.untracked.string("gap"),                quantity = cms.untracked.string("isEBEEGap")),
        cms.PSet( tag = cms.untracked.string("dxyPV"),              quantity = cms.untracked.string("userFloat('dxyPV')")),
        cms.PSet( tag = cms.untracked.string("dzPV"),               quantity = cms.untracked.string("userFloat('dzPV')")),
        cms.PSet( tag = cms.untracked.string("dcot"),               quantity = cms.untracked.string("userFloat('convValueMapProd:dcot')")),
        cms.PSet( tag = cms.untracked.string("dist"),               quantity = cms.untracked.string("userFloat('convValueMapProd:dist')")),
        cms.PSet( tag = cms.untracked.string("mom"),                quantity = cms.untracked.string("?genParticleRef().isNonnull()?genParticleRef().get().mother(0).mother(0).pdgId():9999")),
        cms.PSet( tag = cms.untracked.string("tkDefault"),          quantity = cms.untracked.string("userFloat('tkDefault')")),
        cms.PSet( tag = cms.untracked.string("ecalDefault"),        quantity = cms.untracked.string("userFloat('ecalDefault')")),
        cms.PSet( tag = cms.untracked.string("hcalDefault"),        quantity = cms.untracked.string("userFloat('hcalDefault')")),
    ), 
    weight = cms.untracked.double(1),
    vtxWeights = cms.untracked.PSet(
        vtx05 = cms.untracked.vdouble( weight05[:] ),
        vtx10 = cms.untracked.vdouble( weight10[:] ),
        vtx15 = cms.untracked.vdouble( weight15[:] ),
    ),
    vtxLabel = cms.InputTag('goodPrimaryVertices'),

)  

process.eleSeq = cms.Sequence(
    process.boostedElectrons * (
        process.ele80 + 
        process.ele95
    ) *
    process.eleTuple
)

#  ______ _        _______             _               _____                 
# |  ____| |      |__   __|           | |             |  __ \                
# | |__  | | ___     | |_ __ __ _  ___| | _____ _ __  | |  | | ___ _ __  ___ 
# |  __| | |/ _ \    | | '__/ _` |/ __| |/ / _ \ '__| | |  | |/ _ \ '_ \/ __|
# | |____| |  __/    | | | | (_| | (__|   <  __/ |    | |__| |  __/ |_) \__ \
# |______|_|\___|    |_|_|  \__,_|\___|_|\_\___|_|    |_____/ \___| .__/|___/
#                                                                 | |        
#                                                                 |_|        

#Run Dzs
process.eleIsoDzSeq = cms.Sequence( process.dummy )
for dz in dZValues:
    #name of this iso deposit
    isoDepName = "elTkDZ{0:0<5}".format(dz).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.eleIsoDepositTk.clone())
    #change the dz
    getattr(process,isoDepName).ExtractorPSet.Diff_z = dz
    #add the VM to the boostedElectron
    process.boostedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
    process.boostedElectrons.deposits[-1].src = isoDepName
    process.boostedElectrons.deposits[-1].label = cms.string(isoDepName)
    process.boostedElectrons.deposits[-1].deltaR = 0.3
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
    isoDepName = "elTkD0{0:0<5}".format(d0).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.eleIsoDepositTk.clone())
    #change the d0
    getattr(process,isoDepName).ExtractorPSet.Diff_r = d0
    #add the VM to the boostedElectron
    process.boostedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
    process.boostedElectrons.deposits[-1].src = isoDepName
    process.boostedElectrons.deposits[-1].label = cms.string(isoDepName)
    process.boostedElectrons.deposits[-1].deltaR = 0.3
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


isGood =   ("?(isGlobalMuon && isTrackerMuon &&" +
            " innerTrack.found >10 &&" +
            " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " +
            " globalTrack.normalizedChi2 <10 &&" +
            " globalTrack.hitPattern.numberOfValidMuonHits > 0 && " +
            " numberOfMatches > 1 && " +
            " abs(track.ptError / pt) < 0.10 )?1:0" ) 
isPrompt = ("?( abs(userFloat('dxyPV')) < 0.02 && " +
            "  abs(userFloat('dzPV'))  < 1.0  )?1:0  " )


process.muTuple = cms.EDAnalyzer("BetterMuonTupleProducer", 
    src = cms.untracked.InputTag("boostedMuons"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("pt"),                 quantity = cms.untracked.string("pt")),
        cms.PSet( tag = cms.untracked.string("eta"),                quantity = cms.untracked.string("eta")),
        cms.PSet( tag = cms.untracked.string("phi"),                quantity = cms.untracked.string("phi")),
        cms.PSet( tag = cms.untracked.string("isGood"),             quantity = cms.untracked.string(isGood)),
        cms.PSet( tag = cms.untracked.string("isPrompt"),           quantity = cms.untracked.string(isPrompt)),
        cms.PSet( tag = cms.untracked.string("dxyPV"),              quantity = cms.untracked.string("userFloat('dxyPV')")),
        cms.PSet( tag = cms.untracked.string("dzPV"),               quantity = cms.untracked.string("userFloat('dzPV')")),
        cms.PSet( tag = cms.untracked.string("mom"),                quantity = cms.untracked.string("?genParticleRef().isNonnull()?genParticleRef().get().mother(0).mother(0).pdgId():9999")),
        cms.PSet( tag = cms.untracked.string("tkDep"),              quantity = cms.untracked.string("userFloat('tkDep')")),
        cms.PSet( tag = cms.untracked.string("ecalDep"),            quantity = cms.untracked.string("userFloat('ecalDep')")),
        cms.PSet( tag = cms.untracked.string("hcalDep"),            quantity = cms.untracked.string("userFloat('hcalDep')")),
        cms.PSet( tag = cms.untracked.string("tkDefault"),          quantity = cms.untracked.string("userFloat('tkDefault')")),
        cms.PSet( tag = cms.untracked.string("ecalDefault"),        quantity = cms.untracked.string("userFloat('ecalDefault')")),
        cms.PSet( tag = cms.untracked.string("hcalDefault"),        quantity = cms.untracked.string("userFloat('hcalDefault')")),
    ), 
    weight = cms.untracked.double(1),
    vtxWeights = cms.untracked.PSet(
        vtx05 = cms.untracked.vdouble( weight05[:] ),
        vtx10 = cms.untracked.vdouble( weight10[:] ),
        vtx15 = cms.untracked.vdouble( weight15[:] ),
    ),
    vtxLabel = cms.InputTag('goodPrimaryVertices'),
)  

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

#Run Dzs
process.muIsoDzSeq = cms.Sequence( process.dummy )
for dz in dZValues:
    #name of this iso deposit
    isoDepName = "muTkDZ{0:0<5}".format(dz).replace('.','')
    #make the IsoDeposit Producer and add it to process
    setattr(process,isoDepName,process.muIsoDepositTk.clone())
    #change the dz
    getattr(process,isoDepName).ExtractorPSet.Diff_z = dz
    #add the VM to the boostedElectron
    process.boostedMuons.deposits.append( tempPSet.clone() )
    process.boostedMuons.deposits[-1].src = isoDepName
    process.boostedMuons.deposits[-1].label = cms.string(isoDepName)
    process.boostedMuons.deposits[-1].deltaR = 0.3
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
    process.boostedMuons.deposits.append( tempPSet.clone() )
    process.boostedMuons.deposits[-1].src = isoDepName
    process.boostedMuons.deposits[-1].label = cms.string(isoDepName)
    process.boostedMuons.deposits[-1].deltaR = 0.3
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
process.elePath = cms.Path(~process.elFromW * ~process.elFromTau * process.eleIsoSeq * process.eleSeq)
process.muPath  = cms.Path(~process.muFromW * ~process.muFromTau * process.muIsoSeq  * process.muSeq)
# process.elePath = cms.Path(process.eleIsoSeq * process.eleSeq)
# process.muPath  = cms.Path(process.muIsoSeq  * process.muSeq)

