import FWCore.ParameterSet.Config as cms

######################################
## UserCode/Torino/ZZAnalysis/SkimStep/test/4XY/zzPATSkim.py?revision=1.28&view=markup
## (although in part originally from Adish...)
######################################

#####Code for the eISO
ELE_ISOVETO_CUT="pt>7"
ELE_ISOVETO5_CUT="pt>5"
MUON_ISOVETO_CUT="isGlobalMuon &&pt>5"

#To be used for iso vetos
vetoMuonsZZ4L =  cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string(MUON_ISOVETO_CUT)
)
#To be used for iso vetos, with present default of pt>7
vetoElectronsZZ4L =  cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string(ELE_ISOVETO_CUT)
)
vetoElectronsZZ4L5 =  cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string(ELE_ISOVETO5_CUT)
)

from RecoMuon.MuonIsolationProducers.muIsoDeposits_cff import *
## since the iso deposit produces are not real producers but just copiers, we need to get the real ones too
from RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi import muIsoDepositTk as muIsoDepositTkForReal
from RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi import muIsoDepositCalByAssociatorTowers as muIsoDepositCalByAssociatorTowersForReal
muIsoDepositTkZZ4L = muIsoDepositTkForReal.clone()
muIsoDepositTkZZ4L.IOPSet.inputMuonCollection = "preBoostedMuons"
muIsoDepositCalByAssociatorTowersZZ4L = muIsoDepositCalByAssociatorTowersForReal.clone()
muIsoDepositCalByAssociatorTowersZZ4L.IOPSet.inputMuonCollection = "preBoostedMuons"

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
eleIsoDepositTkZZ4L = eleIsoDepositTk.clone(src = "preBoostedElectrons")
eleIsoDepositEcalFromHitsZZ4L = eleIsoDepositEcalFromHits.clone(src = "preBoostedElectrons")
eleIsoDepositHcalFromTowersZZ4L = eleIsoDepositHcalFromTowers.clone(src = "preBoostedElectrons")

muIsoFromDepsTkZZ4L = cms.PSet(
       label = cms.string("tkZZ4L"),
       mode = cms.string('sum'),
       src = cms.InputTag("muIsoDepositTkZZ4L"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuonsZZ4L:0.015',
                           'vetoElectronsZZ4L:0.015',
                           'Threshold(1.0)'),
       skipDefaultVeto = cms.bool(True)
)
eleIsoFromDepsTkZZ4L = cms.PSet(
       label = cms.string("tkZZ4L"),
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositTkZZ4L"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuonsZZ4L:0.015', 
                           'vetoElectronsZZ4L:RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)', 
                           'RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)', 
                           'Threshold(0.7)'),
       skipDefaultVeto = cms.bool(True)
)


muIsoFromDepsTkZZ4L5 = cms.PSet(
       label = cms.string("tkZZ4Lpt5"),
       mode = cms.string('sum'),
       src = cms.InputTag("muIsoDepositTkZZ4L"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuonsZZ4L:0.015',
                           'vetoElectronsZZ4L5:0.015',
                           'Threshold(1.0)'),
       skipDefaultVeto = cms.bool(True)
)


eleIsoFromDepsTkZZ4L5 = cms.PSet(
       label = cms.string("tkZZ4Lpt5"),
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositTkZZ4L"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuonsZZ4L:0.015', 
                           'vetoElectronsZZ4L5:RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)', 
                           'RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)', 
                           'Threshold(0.7)'),
       skipDefaultVeto = cms.bool(True)
)


muIsoFromDepsEcalZZ4L = cms.PSet(
       label = cms.string("ecalZZ4L"),
       mode = cms.string('sum'),
       src = cms.InputTag("muIsoDepositCalByAssociatorTowersZZ4L","ecal"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuonsZZ4L:0.015',
                           'vetoElectronsZZ4L:0.015'),
       skipDefaultVeto = cms.bool(True)
)



eleIsoFromDepsEcalFromHitsByCrystalZZ4L = cms.PSet(
       label = cms.string("ecalZZ4L"),
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositEcalFromHitsZZ4L"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('NumCrystalVeto(3.0)',
                           'EcalBarrel:NumCrystalEtaPhiVeto(1.0,9999.0)',
                           'EcalEndcaps:NumCrystalEtaPhiVeto(1.5,9999.0)',
                           'EcalBarrel:AbsThresholdFromTransverse(0.08)',
                           'EcalEndcaps:AbsThreshold(0.20)',
                           'vetoMuonsZZ4L:0.05',
                           'vetoElectronsZZ4L:NumCrystalVeto(3.0)',
                           'vetoElectronsZZ4L:NumCrystalEtaPhiVeto(1.5,15.0)'),
       skipDefaultVeto = cms.bool(True)
)

muIsoFromDepsHcalZZ4L = cms.PSet(
       label = cms.string("hcalZZ4L"),
       mode = cms.string('sum'),
       src = cms.InputTag("muIsoDepositCalByAssociatorTowersZZ4L","hcal"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuonsZZ4L:0.015',
                           'vetoElectronsZZ4L:0.015'),
       skipDefaultVeto = cms.bool(True)
)


eleIsoFromDepsHcalFromTowersZZ4L = cms.PSet(
       label = cms.string("hcalZZ4L"),
       src = cms.InputTag("eleIsoDepositHcalFromTowersZZ4L"),
       deltaR = cms.double(0.4),
       weight = cms.string('1'),
       vetos = cms.vstring('0.15', 'vetoMuonsZZ4L:0.05'),
       skipDefaultVeto = cms.bool(True),
       mode = cms.string('sum')
) 

hzz4lDetectorIsoDeposits = cms.Sequence(
    muIsoDepositTkZZ4L +
    muIsoDepositCalByAssociatorTowersZZ4L +
    eleIsoDepositTkZZ4L +
    eleIsoDepositEcalFromHitsZZ4L +
    eleIsoDepositHcalFromTowersZZ4L
)

muIsoFromDepsZZ4L = cms.VPSet(
    muIsoFromDepsTkZZ4L,
    muIsoFromDepsTkZZ4L5,
    muIsoFromDepsEcalZZ4L,
    muIsoFromDepsHcalZZ4L
)
eleIsoFromDepsZZ4L = cms.VPSet(
    eleIsoFromDepsTkZZ4L,
    eleIsoFromDepsTkZZ4L5,
    eleIsoFromDepsEcalFromHitsByCrystalZZ4L,
    eleIsoFromDepsHcalFromTowersZZ4L
)

hzz4lDetectorIsolationSequence = cms.Sequence(
    vetoMuonsZZ4L + vetoElectronsZZ4L + vetoElectronsZZ4L5 +
    hzz4lDetectorIsoDeposits
)


