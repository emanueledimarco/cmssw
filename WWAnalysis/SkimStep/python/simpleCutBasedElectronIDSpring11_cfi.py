import FWCore.ParameterSet.Config as cms

simpleCutBasedElectronID = cms.EDProducer("EleIdCutBasedExtProducer",

    src = cms.InputTag("gsfElectrons"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVertices"),
    dcsTag = cms.InputTag("scalersRawToDigi"),                                          
    algorithm = cms.string('eIDCB'),

    electronIDType  = cms.string('robust'),
    electronQuality = cms.string('80'),
    electronVersion = cms.string(''),
    etBinning = cms.bool(True),
    additionalCategories = cms.bool(True),

    robust60EleIDCuts = cms.PSet(
        #hOverE ieie dPhiIn dEtaIn e2555 e1555 tkIso ecal hcal hcal1 hcal2 tkRel ecalRel hcalRel all allRel ped pedRel minieie eOpMin eOpMax ip nhit nfirst dist dcot
        barrel = cms.vdouble(9999.,0.010,0.020,0.004,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
        endcap = cms.vdouble(9999.,0.031,0.021,0.004,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
    ),

    robust70EleIDCuts = cms.PSet(
        #hOverE ieie dPhiIn dEtaIn e2555 e1555 tkIso ecal hcal hcal1 hcal2 tkRel ecalRel hcalRel all allRel ped pedRel minieie eOpMin eOpMax ip nhit nfirst dist dcot
        barrel = cms.vdouble(9999.,0.010,0.020,0.004,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
        endcap = cms.vdouble(9999.,0.031,0.021,0.005,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
    ),

    robust80EleIDCuts = cms.PSet(
        #hOverE ieie dPhiIn dEtaIn e2555 e1555 tkIso ecal hcal hcal1 hcal2 tkRel ecalRel hcalRel all allRel ped pedRel minieie eOpMin eOpMax ip nhit nfirst dist dcot
        barrel = cms.vdouble(9999.,0.010,0.027,0.005,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
        endcap = cms.vdouble(9999.,0.031,0.021,0.006,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
    ),

    robust85EleIDCuts = cms.PSet(
        #hOverE ieie dPhiIn dEtaIn e2555 e1555 tkIso ecal hcal hcal1 hcal2 tkRel ecalRel hcalRel all allRel ped pedRel minieie eOpMin eOpMax ip nhit nfirst dist dcot
        barrel = cms.vdouble(9999.,0.010,0.039,0.005,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
        endcap = cms.vdouble(9999.,0.031,0.028,0.007,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
    ),

    robust90EleIDCuts = cms.PSet(
        #hOverE ieie dPhiIn dEtaIn e2555 e1555 tkIso ecal hcal hcal1 hcal2 tkRel ecalRel hcalRel all allRel ped pedRel minieie eOpMin eOpMax ip nhit nfirst dist dcot
        barrel = cms.vdouble(9999.,0.010,0.071,0.007,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
        endcap = cms.vdouble(9999.,0.031,0.047,0.011,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
    ),

    robust95EleIDCuts = cms.PSet(
        #hOverE ieie dPhiIn dEtaIn e2555 e1555 tkIso ecal hcal hcal1 hcal2 tkRel ecalRel hcalRel all allRel ped pedRel minieie eOpMin eOpMax ip nhit nfirst dist dcot
        barrel = cms.vdouble(9999.,0.012,0.800,0.007,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
        endcap = cms.vdouble(9999.,0.031,0.700,0.011,-1,-1,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,0.000,-9999.,9999.,9999.,0,-1,0.02,0.02),
    ),

)


