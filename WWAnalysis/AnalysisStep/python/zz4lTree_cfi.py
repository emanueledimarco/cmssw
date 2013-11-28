import FWCore.ParameterSet.Config as cms

zz4lTree = cms.EDFilter("ProbeTreeProducer",
    src = cms.InputTag("skimEvent4L"),
    cut = cms.string("mass > 0"),
    addRunLumiInfo = cms.bool(True),
    variables = cms.PSet(
        hypo = cms.string("hypo()"),
        genhiggsmass = cms.string("getGenHiggsMass()"),
        genhiggsmassweight = cms.string("getHiggsMassWeight()"),
        channel = cms.string("channel()"),
        mass = cms.string("mass"),
        m4l  = cms.string("m4l"),
        pt   = cms.string("pt"),
        rap  = cms.string("rapidity"),
        ## ----------- Z and di-lepton variables ---------------
        z1mass = cms.string("z(0).mass"),
        z1mll  = cms.string("mll(0,0, 0,1)"),
        z1pt   = cms.string("z(0).pt"),
        z1eta  = cms.string("z(0).eta"),
        z1rap  = cms.string("z(0).rapidity"),
        z2mass = cms.string("z(1).mass"),
        z2mll  = cms.string("mll(1,0, 1,1)"),
        z2pt   = cms.string("z(1).pt"),
        z2eta  = cms.string("z(1).eta"),
        z2rap  = cms.string("z(1).rapidity"),
        mll13  = cms.string("mll(0,0, 1,0)"), qll13 = cms.string("lq(0,0) + lq(1,0)"),
        mll14  = cms.string("mll(0,0, 1,1)"), qll14 = cms.string("lq(0,0) + lq(1,1)"),
        mll23  = cms.string("mll(0,1, 1,0)"), qll23 = cms.string("lq(0,1) + lq(1,0)"),
        mll24  = cms.string("mll(0,1, 1,1)"), qll24 = cms.string("lq(0,1) + lq(1,1)"),
        ## ------------- Other variables ------------------
        elPtMin = cms.string("elePtMin"),
        muPtMin  = cms.string("muPtMin"),
        pfmet    = cms.string("pfMet.pt()"),
        rho    = cms.string("getRho()"),
        worsePairCombRelIsoBaseline = cms.string("worsePairCombRelIsoBaseline()"),
        intimeSimVertices = cms.string("intimeSimVertices()"),
        numTrueInteractions = cms.string("numTrueInteractions"),
        recoVertices = cms.string("getNumRecoVertices()"),
        ## ------------- Angular variables ------------------
        melaCosThetaStar = cms.string("getCosThetaStar"),
        melaCosTheta1    = cms.string("getCosTheta1"),
        melaCosTheta2    = cms.string("getCosTheta2"),
        melaPhi          = cms.string("getPhi"),
        melaPhiStar1         = cms.string("getPhiStar1"),
        #
        # other MEKDs 
        ME_SMH_ZZ          = cms.string("userFloat('ME_SMH_ZZ')"),
        #ME_SMH_ZZ_status   = cms.string("userInt('ME_SMH_ZZ_status')"),
        ME_SMH_0minus      = cms.string("userFloat('ME_SMH_0-')"),
        #ME_SMH_0minus_status = cms.string("userInt('ME_SMH_0-_status')"),
        ME_SMH_0hplus        = cms.string("userFloat('ME_SMH_0+h')"),
        #ME_SMH_0hplus_status = cms.string("userInt('ME_SMH_0+h_status')"),
        ME_SMH_1minus      = cms.string("userFloat('ME_SMH_1-')"),
        #ME_SMH_1minus_status = cms.string("userInt('ME_SMH_1-_status')"),
        ME_SMH_1plus        = cms.string("userFloat('ME_SMH_1+')"),
        #ME_SMH_1plus_status = cms.string("userInt('ME_SMH_1+_status')"),
        ME_SMH_2mplus_gg          = cms.string("userFloat('ME_SMH_2+m_gg')"),
        #ME_SMH_2mplus_gg_status   = cms.string("userInt('ME_SMH_2+m_gg_status')"),
        ME_SMH_2mplus_qq          = cms.string("userFloat('ME_SMH_2+m_qq')"),
        #ME_SMH_2mplus_qq_status   = cms.string("userInt('ME_SMH_2+m_qq_status')"),
        # BDT 
        bdtScalarVsBkg_125 = cms.string("userFloat('BDT_ScalarVsBkg_125')"),
        ##  ------------ Mass resolution ------------------
        massErr  = cms.string("userFloat('massErr')"),
        l1massErr = cms.string("userFloat('massErr[0]')"),
        l2massErr = cms.string("userFloat('massErr[1]')"),
        l3massErr = cms.string("? z(0).numberOfDaughters == 2 ? userFloat('massErr[2]') :  userFloat('massErr[3]')"),
        l4massErr = cms.string("? z(0).numberOfDaughters == 2 ? userFloat('massErr[3]') :  userFloat('massErr[4]')"),
        pho1massErr = cms.string("? z(0).numberOfDaughters == 2 ? 0 : userFloat('massErr[2]')"),
        pho2massErr = cms.string("? z(1).numberOfDaughters == 2 ? 0 : "+ ## 4 if Z1 has no FSR, 5 if Z1 has FSR
                                 " (3-z(0).numberOfDaughters)*userFloat('massErr[4]') + (z(0).numberOfDaughters-2)*userFloat('massErr[5]')"),
        z1massErr = cms.string("z(0).masterClone().userFloat('massErr')"),
        z1l1massErr = cms.string("z(0).masterClone().userFloat('massErr[0]')"),
        z1l2massErr = cms.string("z(0).masterClone().userFloat('massErr[1]')"),
        ##  ------------ Z1 Mass refit info ------------------
        z1Refitted = cms.string("userInt('z1Refitted')"),
        z1Refit_chi2In = cms.string("userFloat('z1Refit_chi2In')"),
        z1Refit_chi2Out = cms.string("userFloat('z1Refit_chi2Out')"),
        z1Refit_l1 = cms.string("userFloat('z1Refit_l1')"),
        z1Refit_l2 = cms.string("userFloat('z1Refit_l2')"),
        z1Refit_z1mllBefore = cms.string("userFloat('z1Refit_z1mllBefore')"),
        z1Refit_z1massBefore = cms.string("userFloat('z1Refit_z1massBefore')"),
        z1Refit_z1mllAfter = cms.string("userFloat('z1Refit_z1mllAfter')"),
        z1Refit_z1massAfter = cms.string("userFloat('z1Refit_z1massAfter')"),
        z1Refit_massBefore = cms.string("userFloat('z1Refit_massBefore')"),
        z1Refit_massErrBefore = cms.string("userFloat('z1Refit_massErrBefore')"),
        z1Refit_massErrL12 = cms.string("userFloat('z1Refit_massErrL12')"),
        ## ------------- Photon variables ------------------
        pho1pt   = cms.string("? z(0).numberOfDaughters == 2 ? 0 : l(0,2).pt"),
        pho1eta  = cms.string("? z(0).numberOfDaughters == 2 ? 0 : l(0,2).eta"),
        pho1phi  = cms.string("? z(0).numberOfDaughters == 2 ? 0 : l(0,2).phi"),
        pho1dr   = cms.string("? z(0).numberOfDaughters == 2 ? 0 : lval(0,2,'deltaR(overlaps(\"goodLepNoIso\")[0].eta, overlaps(\"goodLepNoIso\")[0].phi, eta, phi)')"),
        pho1iso  = cms.string("? z(0).numberOfDaughters == 2 ? 0 : lval(0,2,'(userFloat(\"fsrPhotonPFIsoChHad03pt02\")+userFloat(\"fsrPhotonPFIsoPhoton03\")+userFloat(\"fsrPhotonPFIsoNHad03\")+userFloat(\"fsrPhotonPFIsoChHadPU03pt02\"))/pt')"),
        pho2pt   = cms.string("? z(1).numberOfDaughters == 2 ? 0 : l(1,2).pt"),
        pho2eta  = cms.string("? z(1).numberOfDaughters == 2 ? 0 : l(1,2).eta"),
        pho2phi  = cms.string("? z(1).numberOfDaughters == 2 ? 0 : l(1,2).phi"),
        pho2dr   = cms.string("? z(1).numberOfDaughters == 2 ? 0 : lval(1,2,'deltaR(overlaps(\"goodLepNoIso\")[0].eta, overlaps(\"goodLepNoIso\")[0].phi, eta, phi)')"),
        pho2iso  = cms.string("? z(1).numberOfDaughters == 2 ? 0 : lval(1,2,'(userFloat(\"fsrPhotonPFIsoChHad03pt02\")+userFloat(\"fsrPhotonPFIsoPhoton03\")+userFloat(\"fsrPhotonPFIsoNHad03\")+userFloat(\"fsrPhotonPFIsoChHadPU03pt02\"))/pt')"),
        ## --------- Jet variables ------
        njets30 = cms.string("njets(30,4.7)"),
        jet1pt  = cms.string("? njets(0,4.7) > 0 ? jet(0,0.0,4.7).pt  : 0"),
        jet1eta = cms.string("? njets(0,4.7) > 0 ? jet(0,0.0,4.7).eta : 0"),
        jet1phi = cms.string("? njets(0,4.7) > 0 ? jet(0,0.0,4.7).phi : 0"),
        jet2pt  = cms.string("? njets(0,4.7) > 1 ? jet(1,0.0,4.7).pt  : 0"),
        jet2eta = cms.string("? njets(0,4.7) > 1 ? jet(1,0.0,4.7).eta : 0"),
        jet2phi = cms.string("? njets(0,4.7) > 1 ? jet(1,0.0,4.7).phi : 0"),
        mjj     = cms.string("? njets(0,4.7) > 1 ? mjj(0,1,0.0,4.7) : 0"),
    ),
    flags = cms.PSet(
        # ===========================
        threeMassCut12SF = cms.string("hypo == 4 || hypo == 5 || nGoodPairs(\"mass > 12\", 0) >= 3"),
        fourMassCut4SF   = cms.string("hypo == 4 || hypo == 5 || nGoodPairs(\"mass > 4\", 0) >= 4"),
        fourMassCut4Any  = cms.string("nGoodPairs(\"mass > 4\", 0) >= 4"),
        fourMassCut4AnyS = cms.string("nGoodPairs(\"mass > 4\", 1) >= 6"),
        # ===========================
        mc_l1 = cms.string("genl(0,0).isNonnull"),
        mc_l2 = cms.string("genl(0,1).isNonnull"),
        mc_l3 = cms.string("genl(1,0).isNonnull"),
        mc_l4 = cms.string("genl(1,1).isNonnull"),
        mc_z1 = cms.string("genz(0).isNonnull"),
        mc_z2 = cms.string("genz(0).isNonnull"),
        mc_4lany = cms.string("nMatchedLeptons == 4"),
        mc_4l    = cms.string("hasCorrectAssignment"),
    ),
)
def addLepVar(zz4l, name, expr) :
    setattr(zz4l.variables, "l1"+name, cms.string(expr % (0, 0)));
    setattr(zz4l.variables, "l2"+name, cms.string(expr % (0, 1)));
    setattr(zz4l.variables, "l3"+name, cms.string(expr % (1, 0)));
    setattr(zz4l.variables, "l4"+name, cms.string(expr % (1, 1)));
def addLepCut(zz4l, name, expr) :
    setattr(zz4l.flags, "l1"+name, cms.string(expr % (0, 0)));
    setattr(zz4l.flags, "l2"+name, cms.string(expr % (0, 1)));
    setattr(zz4l.flags, "l3"+name, cms.string(expr % (1, 0)));
    setattr(zz4l.flags, "l4"+name, cms.string(expr % (1, 1)));
addLepVar(zz4lTree, "pt", "lpt(%d,%d)")
addLepVar(zz4lTree, "eta", "leta(%d,%d)")
addLepVar(zz4lTree, "phi", "lphi(%d,%d)")
addLepVar(zz4lTree, "q", "lq(%d,%d)")
addLepVar(zz4lTree, "dz", "ldz(%d,%d)")
addLepVar(zz4lTree, "ip2d", "lip2d(%d,%d)")
addLepVar(zz4lTree, "sip3d", "lsip3d(%d,%d)")
addLepVar(zz4lTree, "pdgId", "lpdgId(%d,%d)")
addLepVar(zz4lTree, "trig",  "lfiresTrigger(%d,%d)")
#addLepVar(zz4lTree, "chIso", "lisoTrkCustom(%d,%d)")
#addLepVar(zz4lTree, "neuIso", "lisoNeuCustom(%d,%d)")
#addLepVar(zz4lTree, "phoIso", "lisoPhoCustom(%d,%d)")
#addLepVar(zz4lTree, "trkIsoBaseline", "lisoTrkBaseline(%d,%d)")
#addLepVar(zz4lTree, "ecalIsoBaseline", "lisoEcalBaseline(%d,%d)")
#addLepVar(zz4lTree, "hcalIsoBaseline", "lisoHcalBaseline(%d,%d)")
#addLepVar(zz4lTree, "ecalIsoBaselineRaw", "luserFloat(%d,%d, 'ecalZZ4L')")
#addLepVar(zz4lTree, "hcalIsoBaselineRaw", "luserFloat(%d,%d, 'hcalZZ4L')")
#addLepVar(zz4lTree, "combRelIsoBaseline", "lisoCombRelBaseline(%d,%d)")
#addLepVar(zz4lTree, "combRelIsoBaseline", "lisoCombRelBaseline(%d,%d)")
addLepVar(zz4lTree, "idNew", "luserInt(%d,%d,'newID')")
addLepVar(zz4lTree, "idPRL", "luserInt(%d,%d,'prlID')")
addLepVar(zz4lTree, "pfIsoChHad04",      "luserFloat(%d,%d, 'pfChHadIso04')")
addLepVar(zz4lTree, "pfPhotonIso04",      "luserFloat(%d,%d, 'pfPhotonIso04')")
addLepVar(zz4lTree, "pfNHadIso04",      "luserFloat(%d,%d, 'pfNHadIso04')")
addLepVar(zz4lTree, "pfChHadPU04",      "luserFloat(%d,%d, 'pfChHadPU04')")
addLepVar(zz4lTree, "pfIsoComb04EACorr", "luserFloat(%d,%d, 'pfCombIso04EACorr')")
addLepVar(zz4lTree, "pfCombRelIso04dBCorr",      "luserFloat(%d,%d, 'pfCombRelIso04dBCorr')")
addLepVar(zz4lTree, "bdtIso", "luserFloat(%d,%d, 'bdtIso')")
addLepVar(zz4lTree, "nmisshits", "nmisshits(%d,%d)")
