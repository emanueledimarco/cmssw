## Jets 
from RecoJets.JetProducers.ak5PFJets_cfi import *
ak5PFJetsChargedDZ = ak5PFJets.clone( src = "reducedPFCands", jetType = "BasicJet" )

## IsoDeposits 
isoDepJetMuonsJetCone = cms.EDProducer("CandJetIsoDepMapProd",
    leptLabel  = cms.InputTag("boostedMuons"),
    vtxLabel   = cms.InputTag("goodPrimaryVertices"),
    jetLabel   = cms.InputTag("ak5PFJetsChargedDZ"),
    coneMode   = cms.string("jet"),
    drVeto     = cms.double(0.01),
)
isoDepJetMuonsFixCone = isoDepJetMuonsJetCone.clone(
    candLabel = cms.InputTag("reducedPFCands"),
    coneMode   = cms.string("fixed"),
    coneSize   = cms.double(0.7),
    dzCut      = cms.double(999), # already done in reducedPFCands
)
isoDepJetMuonsMaxCone03 = isoDepJetMuonsFixCone.clone(
    coneMode   = cms.string("max"),
    coneSize   = cms.double(0.3),
)
isoDepJetMuonsMaxCone04 = isoDepJetMuonsMaxCone03.clone(coneSize = 0.4)
isoDepJetMuonsMaxCone05 = isoDepJetMuonsMaxCone03.clone(coneSize = 0.5)
isoDepJetMuonsSumCone03 = isoDepJetMuonsFixCone.clone(
    coneMode   = cms.string("sum"),
    coneSize   = cms.double(0.3),
)
isoDepJetMuonsSumCone02 = isoDepJetMuonsSumCone03.clone(coneSize = 0.2)
isoDepJetMuonsSumCone04 = isoDepJetMuonsSumCone03.clone(coneSize = 0.4)

isoDepJetElectronsJetCone   = isoDepJetMuonsJetCone.clone(leptLabel = cms.InputTag("boostedElectrons"))
isoDepJetElectronsFixCone   = isoDepJetMuonsFixCone.clone(leptLabel = cms.InputTag("boostedElectrons"))
isoDepJetElectronsMaxCone03 = isoDepJetMuonsMaxCone03.clone(leptLabel = cms.InputTag("boostedElectrons"))
isoDepJetElectronsMaxCone04 = isoDepJetMuonsMaxCone04.clone(leptLabel = cms.InputTag("boostedElectrons"))
isoDepJetElectronsMaxCone05 = isoDepJetMuonsMaxCone05.clone(leptLabel = cms.InputTag("boostedElectrons"))
isoDepJetElectronsSumCone02 = isoDepJetMuonsSumCone02.clone(leptLabel = cms.InputTag("boostedElectrons"))
isoDepJetElectronsSumCone03 = isoDepJetMuonsSumCone03.clone(leptLabel = cms.InputTag("boostedElectrons"))
isoDepJetElectronsSumCone04 = isoDepJetMuonsSumCone04.clone(leptLabel = cms.InputTag("boostedElectrons"))

isoDepJetSeq = cms.Sequence(
    isoDepJetMuonsJetCone   + isoDepJetElectronsJetCone +
    isoDepJetMuonsFixCone   + isoDepJetElectronsFixCone +
    isoDepJetMuonsMaxCone03 + isoDepJetElectronsMaxCone03 +
    isoDepJetMuonsMaxCone04 + isoDepJetElectronsMaxCone04 +
    isoDepJetMuonsMaxCone05 + isoDepJetElectronsMaxCone05 +
    isoDepJetMuonsSumCone02 + isoDepJetElectronsSumCone02 +
    isoDepJetMuonsSumCone03 + isoDepJetElectronsSumCone03 +
    isoDepJetMuonsSumCone04 + isoDepJetElectronsSumCone04 
)

## Good leptons, for isolation vetos
from WWAnalysis.AnalysisStep.wwMuons_cfi     import MUON_ID_CUT, MUON_MERGE_ISO, MUON_MERGE_IP
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ID, ELE_MERGE_ISO, ELE_MERGE_IP, ELE_MERGE_CONV
goodMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string( 'pt > 10 && ' + MUON_ID_CUT + " && " + MUON_MERGE_IP)
)
goodElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string( 'pt > 10 && ' + ELE_MERGE_ID + " && " + ELE_MERGE_IP + " && " + ELE_MERGE_CONV )
)

## ReBoost
from WWAnalysis.AnalysisStep.isoAdding_cff import *
boostedElectrons = isoAddedElectrons.clone( electronTag = "boostedElectrons" )
boostedMuons     = isoAddedMuons.clone( muonTag = "boostedMuons" )
for algo,dr in [ ("JetCone", 9.9), 
                 ("FixCone", 0.3),   ("FixCone", 0.4), 
                 ("MaxCone03", 9.9), ("MaxCone04", 9.9), 
                 ("SumCone02", 9.9), ("SumCone03", 9.9), ("SumCone04", 9.9) ]:
    for (name, content, vetos) in [ 
                        ("Charged",         "charged",    [ 'goodMuons:0.01','goodElectrons:0.01' ]), 
                        ("ChargedNoOvRem",  "charged"   , [ ]), 
                        #("NeutralHadAll",   "neutralHad", ['goodElectrons:0.07']),
                        #("NeutralHadPt05",  "neutralHad", ['Threshold(0.5)','goodElectrons:0.07']),
                        #("NeutralHadPt1",   "neutralHad", ['Threshold(1.0)','goodElectrons:0.07']),
                        #("Photons",         "photon",     ['goodElectrons:RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)']),
                        #("PhotonsMuStrip",  "photon",     ['goodMuons:RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)', 'goodElectrons:RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)']) 
                        ]:
        algoName = "%s0%d" % (algo, 10*dr) if algo == "FixCone" else algo;
        boostedMuons.deposits.append( cms.PSet(
            label  = cms.string("jetIso%s%s" % (algoName, name)),
            src    = cms.InputTag("isoDepJetMuons"+algo,content),
            deltaR = cms.double(dr),
            vetos  = cms.vstring(vetos[:]),
            weight = cms.string('1'),
            skipDefaultVeto = cms.bool(True),
            mode            = cms.string('sum'),
        ) )
        boostedElectrons.deposits.append( cms.PSet(
            label  = cms.string("jetIso%s%s" % (algoName, name)),
            src    = cms.InputTag("isoDepJetElectrons"+algo,content),
            deltaR = cms.double(dr),
            vetos  = cms.vstring(vetos[:]),
            weight = cms.string('1'),
            skipDefaultVeto = cms.bool(True),
            mode            = cms.string('sum'),
        ) )

isoStudySequence  = cms.Sequence(
    ak5PFJetsChargedDZ +
    goodMuons + goodElectrons +
    isoDepJetSeq +
    boostedElectrons + boostedMuons
)

