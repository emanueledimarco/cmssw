import FWCore.ParameterSet.Config as cms

fsrPhotonsNoEleSC = cms.EDProducer("PATPFParticleCleaner",
    src = cms.InputTag("boostedFsrPhotons"),
    preselection = cms.string("pt > 2 && abs(eta) < 2.4"),
    checkOverlaps = cms.PSet(
        eleVeto = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string("abs(pdgId) == 11"),
            deltaR = cms.double(2.05), # wider than sqrt(2^2 + 0.05*2)
            pairCut = cms.string("abs(deltaPhi(ele.phi,pf.phi)) < 2.0 && abs(pf.eta - ele.eta) < 0.05 || deltaR(ele.eta,ele.phi,pf.eta,pf.phi) < 0.15"),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False), ## Will Veto these photons later
        ),
    ),
    finalCut = cms.string(""),
)

## Definition of PF relative isolation with deltaBeta correction
_netural_iso_uncorr = "pf.userFloat('fsrPhotonPFIsoPhoton03') + pf.userFloat('fsrPhotonPFIsoNHad03')"
_neutral_iso_dbeta  = "max(0, %s - 0.5*pf.userFloat('fsrPhotonPFIsoChHadPU03pt02'))" % _netural_iso_uncorr
_comb_rel_iso_dbeta = "(pf.userFloat('fsrPhotonPFIsoChHad03pt02') + %s)/pf.pt" % _neutral_iso_dbeta
## Definition of PF relative isolation including pileup
_comb_iso_all = "+".join(["pf.userFloat('%s')"%x for x in ('fsrPhotonPFIsoChHad03pt02','fsrPhotonPFIsoPhoton03','fsrPhotonPFIsoNHad03','fsrPhotonPFIsoChHadPU03pt02')])
_comb_rel_iso_all = "(%s)/pf.pt" % _comb_iso_all

fsrPhotonsStep2 = cms.EDProducer("PATPFParticleCleaner",
    src = cms.InputTag("fsrPhotonsNoEleSC"),
    preselection = cms.string(""),
    checkOverlaps = cms.PSet(
        ## Attach by deltaR, pT and isolation
        goodLepIso = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            pairCut = cms.string("deltaR < 0.07 || (pf.pt > 4 && %s < 1.0)" % _comb_rel_iso_all),
            deltaR = cms.double(0.5),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False),
        ),
        ## Attach by deltaR, pT and isolation with deltabeta (old)
        goodLepIsoDBeta = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            pairCut = cms.string("deltaR < 0.07 || (pf.pt > 4 && %s < 1.0)" % _comb_rel_iso_dbeta),
            deltaR = cms.double(0.5),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False),
        ),
        ## Just attach by deltaR and pT, for debugging and isolation studies
        goodLepNoIso = cms.PSet(
            src = cms.InputTag("goodLep"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            pairCut = cms.string("deltaR < 0.07 || pf.pt > 4"),
            deltaR = cms.double(0.5),
            checkRecoComponents = cms.bool(False),
            requireNoOverlaps   = cms.bool(False),
        ),
    ),
    finalCut = cms.string(""),
)

fsrPhotonSeq = cms.Sequence(fsrPhotonsNoEleSC * fsrPhotonsStep2)

fsrPhotonsCR = fsrPhotonsNoEleSC.clone()
fsrPhotonsCR.checkOverlaps.eleVeto.src = 'looseLepCR'
fsrPhotonsCR.checkOverlaps.goodLepIso = fsrPhotonsStep2.checkOverlaps.goodLepIso.clone(src = 'looseLepCR')
fsrPhotonsCR.checkOverlaps.goodLepNoIso = fsrPhotonsStep2.checkOverlaps.goodLepNoIso.clone(src = 'looseLepCR')


