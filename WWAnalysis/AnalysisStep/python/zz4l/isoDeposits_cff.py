import FWCore.ParameterSet.Config as cms

# make the actual IsoDeposits
pfCharged = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("reducedPFCands"),
    cut = cms.string("charge() != 0"),
)
pfNeutral = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("reducedPFCands"),
    cut = cms.string("charge() == 0 && abs(pdgId) != 22"),
)
pfPhotons = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("reducedPFCands"),
    cut = cms.string("abs(pdgId) != 22"),
)
isoDepMuFromPFCharged = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("recMu"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        inputCandView = cms.InputTag("pfCharged"),
        Diff_z = cms.double(999.9),
        DR_Max = cms.double(0.5),
        Diff_r = cms.double(99999.99),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)
isoDepMuFromPFNeutral = isoDepMuFromPFCharged.clone()
isoDepMuFromPFPhotons = isoDepMuFromPFCharged.clone()
isoDepMuFromPFNeutral.ExtractorPSet.inputCandView = "pfNeutral"
isoDepMuFromPFPhotons.ExtractorPSet.inputCandView = "pfPhotons"
#isoDepElFromPFCharged = isoDepMuFromPFCharged.clone(src = "cleanedEl")
#isoDepElFromPFNeutral = isoDepMuFromPFNeutral.clone(src = "cleanedEl")
#isoDepElFromPFPhotons = isoDepMuFromPFPhotons.clone(src = "cleanedEl")
isoDepElFromPFCharged = isoDepMuFromPFCharged.clone(src = "recEl")
isoDepElFromPFNeutral = isoDepMuFromPFNeutral.clone(src = "recEl")
isoDepElFromPFPhotons = isoDepMuFromPFPhotons.clone(src = "recEl")

isoDepositSeq = cms.Sequence(
    pfCharged + pfNeutral + pfPhotons +
    isoDepMuFromPFCharged + isoDepElFromPFCharged +
    isoDepMuFromPFPhotons + isoDepElFromPFPhotons +
    isoDepMuFromPFNeutral + isoDepElFromPFNeutral
)


