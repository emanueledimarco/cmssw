import FWCore.ParameterSet.Config as cms

MU_ID_PF="userInt('pfMuId')>0"
MU_ID_PRL="isGlobalMuon && track.numberOfValidHits > 10"

#EL_BDT="bdtnontrig"
EL_BDT="bdtnontrigAfterEScale"
EL_ID_BDT=("((pt<7 && ((isEB && sigmaIetaIeta<0.01) || (isEE && sigmaIetaIeta<0.03)) ) || "+
            "(pt <= 10 && pt > 7) && (       abs(superCluster.eta) <  0.8    && userFloat('bdtID') > 0.47  ||"+
                        " 0.8 <= abs(superCluster.eta) <  1.479  && userFloat('bdtID') > 0.004 || "+
                        "        abs(superCluster.eta) >= 1.479  && userFloat('bdtID') > 0.295) || "+
            "pt > 10 && (        abs(superCluster.eta) <  0.8    && userFloat('bdtID') > -0.34  ||"+
                        " 0.8 <= abs(superCluster.eta) <  1.479  && userFloat('bdtID') > -0.65 || "+
                        "        abs(superCluster.eta) >=  1.479 && userFloat('bdtID') > 0.6))")
EL_CONV="gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1"
EL_ID_NEW=(EL_ID_BDT+" && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")
EL_ID_PRL=("test_bit(electronID('cicTight'),0) == 1 && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1")

## ==== Alternative set of cuts optimized for MVA iso (and electron ID) from Si and company
MU_MVA_ISO=("(isGlobalMuon && isTrackerMuon && (" +
                 "pt <  10 && abs(eta) <  1.5 && userFloat('bdtIso') > -0.593 ||"+
                 "pt >= 10 && abs(eta) <  1.5 && userFloat('bdtIso') >  0.337 ||"+
                 "pt <  10 && abs(eta) >= 1.5 && userFloat('bdtIso') > -0.767 ||"+
                 "pt >= 10 && abs(eta) >= 1.5 && userFloat('bdtIso') >  0.410) ||"+
              " !isGlobalMuon && isTrackerMuon && userFloat('bdtIso') > -0.989 ||"+
              "isGlobalMuon &&  !isTrackerMuon && userFloat('bdtIso') > -0.995)")


EL_MVA_ID=("(pt <  10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.369 ||"+
                        " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') > -0.025 || "+
                        "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.531 ) || "+
            "pt >= 10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.735 ||"+
                        " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') >  0.467 || "+
                        "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.795 ))")
EL_MVA_ISO=("(pt < 10 && (       abs(eta) <  0.8    && userFloat('bdtIso') >  0.385 ||"+
                         "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') > -0.083 || "+
                         "       abs(eta) >= 1.479  && userFloat('bdtIso') > -0.573) || "+
             "pt >= 10 && (      abs(eta) <  0.8    && userFloat('bdtIso') >  0.413||"+
                         "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') >  0.271|| "+
                         "       abs(eta) >= 1.479  && userFloat('bdtIso') >  0.135))")

EL_MVA_ID_TIGHT=("(pt <  10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.093 ||"+
                              " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') >  0.451 || "+
                              "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.595 ) || "+
                  "pt >= 10 && (       abs(eta) <  0.8    && userFloat('bdtID') >  0.881 ||"+
                              " 0.8 <= abs(eta) <  1.479  && userFloat('bdtID') >  0.731 || "+
                              "        abs(eta) >= 1.479  && userFloat('bdtID') >  0.819 ))")
EL_MVA_ISO_TIGHT=("(pt < 10 && (       abs(eta) <  0.8    && userFloat('bdtIso') >  0.553 ||"+
                               "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') > -0.237 || "+
                               "       abs(eta) >= 1.479  && userFloat('bdtIso') > -0.573) || "+
                   "pt >= 10 && (      abs(eta) <  0.8    && userFloat('bdtIso') >  0.521||"+
                               "0.8 <= abs(eta) <  1.479  && userFloat('bdtIso') >  0.531|| "+
                               "       abs(eta) >= 1.479  && userFloat('bdtIso') >  0.493))")



 
##================================================================================================0

## Compute EA-corrected isolations
boostedElectronsID = cms.EDProducer("PatElectronBoosterBDTID",
    src = cms.InputTag("boostedElectrons"),
    postfix = cms.string("AfterEScale"),
    doTrigMVA = cms.bool(False)
)

boostedElectronsEAPFIso = cms.EDProducer("PatElectronEffAreaIso",
    src = cms.InputTag("boostedElectronsID"),
    rho = cms.string("rhoElActiveArea"),
    deltaR = cms.string("04"),
    label = cms.string("pfCombIso04EACorr"),
    chargedOption = cms.string(""),  # postfix to the userFloat values
    neutralsOption = cms.string(""), # postfix to the userFloat value
    effectiveAreaTarget = cms.string("Data2011"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
)
boostedMuonsEAPFIso = cms.EDProducer("PatMuonEffAreaIso",
    src = cms.InputTag("boostedMuons"),
    rho = cms.string("rhoMuActiveArea"),
    deltaR = cms.string("04"),
    chargedOption = cms.string(""),  # postfix to the userFloat values
    neutralsOption = cms.string("pt05"), # postfix to the userFloat value
    label = cms.string("pfCombIso04EACorr"),
    effectiveAreaTarget = cms.string("Data2011"),
    separatePhotonAndHadronEAs = cms.bool(False), # use total EA
    truncateAtZero = cms.string("yes"), # (yes|no) for total EA, (both|sum|no) for separate EA
)

##================================================================================================0

## Embed final values and discriminators. Beyond here, fixups should be invisible
boostedElectronsStep2 = cms.EDProducer("PatElectronUserFloatAdder",
    src = cms.InputTag("boostedElectronsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfCombRelIso04dBCorr = cms.string("(userFloat('electronPFIsoChHad04') + max(userFloat('electronPFIsoNHad04') + userFloat('electronPFIsoPhoton04') - 0.5*userFloat('electronPFIsoChHadPU04'), 0.))/pt"),
        pfChHadRelIso04 = cms.string("userFloat('electronPFIsoChHad04')/pt"),
        pfChHadPU04 = cms.string("userFloat('electronPFIsoChHadPU04')"),
        pfChHadIso04 = cms.string("userFloat('electronPFIsoChHad04')"),
        pfNHadIso04 = cms.string("userFloat('electronPFIsoNHad04')"),
        pfPhotonIso04 = cms.string("userFloat('electronPFIsoPhoton04')"),
        bdtID     = cms.string("userFloat('%s')" % EL_BDT),
        bdtIDTrig = cms.string("userFloat('bdttrig')"),
        bdtIso = cms.string("userFloat('bdtisonontrig')"),
        sip   = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    flags = cms.PSet(
        convRej = cms.string(EL_CONV),
        bdtID = cms.string(EL_ID_BDT), 
        newID = cms.string(EL_ID_NEW), 
        prlID = cms.string(EL_ID_PRL),
        mvaID  = cms.string(EL_MVA_ID), 
        mvaIso = cms.string(EL_MVA_ISO), 
        mvaIDTight  = cms.string(EL_MVA_ID_TIGHT), 
        mvaIsoTight = cms.string(EL_MVA_ISO_TIGHT), 
    )
)
boostedMuonsStep2 = cms.EDProducer("PatMuonUserFloatAdder",
    src = cms.InputTag("boostedMuonsEAPFIso"),
    variables = cms.PSet(
        pfCombRelIso04EACorr = cms.string("userFloat('pfCombIso04EACorr')/pt"),
        pfCombRelIso04dBCorr = cms.string("(userFloat('muonPFIsoChHad04') + max(userFloat('muonPFIsoNHad04pt05') + userFloat('muonPFIsoPhoton04pt05') - 0.5*pfIsolationR04().sumPUPt, 0.))/pt"),
        pfChHadRelIso04 = cms.string("userFloat('muonPFIsoChHad04')/pt"),
        pfChHadIso04 = cms.string("userFloat('muonPFIsoChHad04')"),
        pfNHadIso04 = cms.string("userFloat('muonPFIsoNHad04pt05')"),
        pfPhotonIso04 = cms.string("userFloat('muonPFIsoPhoton04pt05')"),
        pfChHadPU04 = cms.string("userFloat('muonPFIsoChHadPU04')"),
        bdtIso = cms.string("userFloat('bdtisonontrigPFNOPU')"),
        sip    = cms.string("userFloat('ip')/userFloat('ipErr')"),
    ),
    flags = cms.PSet(
        prlID = cms.string(MU_ID_PRL),
        newID = cms.string(MU_ID_PF),
        mvaIso = cms.string(MU_MVA_ISO), 
    )
)

from WWAnalysis.AnalysisStep.zz4l.fsr_cff import *

#reboosting = cms.Sequence(
#    boostedElectronsID * boostedElectronsEAPFIso   *  boostedElectrons +
#    boostedMuonsEAPFIso * boostedMuons 
#)

reboosting = cms.Sequence(
    boostedElectronsID +
    boostedMuonsEAPFIso * boostedMuonsStep2 
)


#### CUT FLOW BUILDING BLOCKS

SKIM_SEQUENCE = 'reskimNoOS'   # their default (20/10  and a 40 SF pair with mll > 40 but no OS request)
TRIGGER_FILTER = None 

EL_PT_MIN=7;  EL_PT_MIN_LOOSE=7
MU_PT_MIN=5;  MU_PT_MIN_LOOSE=5

MU_PRESELECTION = ("abs(eta) < 2.4 && (isGlobalMuon || numberOfMatches > 0)") ## the PT cut
EL_PRESELECTION = ("abs(eta) < 2.5")                                          ## is below
SINGLE_SIP_CUT        = "userFloat('sip') < 4"
SINGLE_DXYZ_CUT_LOOSE = "abs(userFloat('dzPV')) < 1.0 && userFloat('tip') < 0.5"
SINGLE_SIP_CUT_LOOSE = "userFloat('sip') < 100"
PAIR_SIP_CUT       = "luserFloat(0,'sip') < 4   && luserFloat(1,'sip') < 4  "
PAIR_SIP_CUT_LOOSE = "luserFloat(0,'sip') < 100 && luserFloat(1,'sip') < 100"

SINGLE_ID_NEW = "userInt('newID')"
SINGLE_ID_MVA = "userInt('mvaID')" # to go jointly with MVA iso 
SINGLE_ID_MVA_TIGHT = "userInt('mvaIDTight')" # to go jointly with MVA iso 
PAIR_ID_NEW   = "luserInt(0, 'newID') && luserInt(1, 'newID')"

SINGLE_PFISO_1D_LOOSE="userFloat('pfChHadRelIso04') < 0.7"
#SINGLE_PFISO_1D="userFloat('pfCombRelIso04EACorr') < 0.4"
SINGLE_PFISO_1D="(abs(pdgId)==11 && userFloat('pfCombRelIso04EACorr') < 0.4) || (abs(pdgId)==13 && userFloat('pfCombRelIso04dBCorr') < 0.4)"
SINGLE_MVA_ISO="userInt('mvaIso')"
SINGLE_MVA_ISO_TIGHT="userInt('mvaIsoTight')"
PAIR_PFISO_1D="luserFloat(0,'pfCombRelIso04EACorr') < 0.4 && luserFloat(1,'pfCombRelIso04EACorr') < 0.4"
PAIR_PFISO_2D="luserFloat(0,'pfCombRelIso04EACorr') + luserFloat(1,'pfCombRelIso04EACorr') < 0.35"
PAIR_DETISO_2D="combinedPairRelativeIso() < 0.35"

JET_ETA_MAX = 4.7;
JET_ID   = ( "neutralMultiplicity+chargedMultiplicity > 0 && " +
               "neutralEmEnergyFraction < 0.99 && neutralHadronEnergyFraction < 0.99 && "+
               "(abs(eta) >= 2.4 || "+
                    "chargedEmEnergyFraction < 0.99 && "+
                    "chargedHadronEnergyFraction > 0 && "+
                    "chargedMultiplicity > 0)" )
JET_PUID = "test_bit(userInt('jetId'),2)"

#### CUT FLOW
MUID_LOOSE_NO_PT_CUT = " && ".join([MU_PRESELECTION, SINGLE_DXYZ_CUT_LOOSE])
ELID_LOOSE_NO_PT_CUT = " && ".join([EL_PRESELECTION, SINGLE_DXYZ_CUT_LOOSE, EL_CONV])
MUID_LOOSE = "pt > %f && %s" % (MU_PT_MIN_LOOSE, MUID_LOOSE_NO_PT_CUT)
ELID_LOOSE = "pt > %f && %s" % (EL_PT_MIN_LOOSE, ELID_LOOSE_NO_PT_CUT)

# for later
#MUID_MID = "(" + (" || ".join([SINGLE_PFISO_1D_LOOSE])) + ")"
#ELID_MID = "(" + (" || ".join([SINGLE_PFISO_1D_LOOSE + " && "+ SINGLE_ID_MVA, SINGLE_MVA_ISO])) + ")"

MUID_GOOD = " && ".join(["pt > %f" % MU_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])
ELID_GOOD = " && ".join(["pt > %f" % EL_PT_MIN, SINGLE_SIP_CUT, SINGLE_ID_NEW, SINGLE_PFISO_1D])
JETID_GOOD = " && ".join(["abs(eta) < %f" % JET_ETA_MAX, JET_ID, JET_PUID])

SEL_ANY_Z = "daughter(0).pdgId = - daughter(1).pdgId"

ARBITRATE_EARLY = True # True = PRL-logic; False = keep all candidates until the end
                       # Baseline logic is True; Current synch twiki is False

FOUR_LEPTON_FILTER_PRE_Z  = False # plug a 4-lepton count filter before making Z's
FOUR_LEPTON_FILTER_POST_Z = True # plug a 4-lepton count filter after  making Z's

DO_FSR_RECOVERY = False
FSR_MATCH = "goodLepIso"

SEL_BEST_Z1 = "40 < mass < 120" # Not used if ARBITRATE_EARLY = False

SEL_ZZ4L_STEP_1 = "4 < mz(1) < 120"
SEL_ZZ4L_STEP_2 = "lByPt(0).pt > 20 && lByPt(1).pt > 10"
#SEL_ZZ4L_STEP_3 = "nGoodPairs(\"mass > 4\", 1) >= 6"
SEL_ZZ4L_STEP_3 = "nGoodPairs(\"mass > 4\", 0) >= 4"
SEL_ZZ4L_STEP_4 = ""
SEL_ZZ4L_STEP_5 = "mz(1) > 12"
SEL_ZZ4L_STEP_6 = ""

SEL_ZZ4L_ARBITRATION_1 = "-abs(mz(0)-91.188)"
SEL_ZZ4L_ARBITRATION_2 = "daughter(1).daughter(0).pt + daughter(1).daughter(1).pt"

SEL_ZZ4L_VBF_STEP_1 = "njets(30,4.7) >= 1"
SEL_ZZ4L_VBF_STEP_2 = "njets(30,4.7) >= 2"
SEL_ZZ4L_VBF_STEP_3 = "njets(30,4.7) == 2"
SEL_ZZ4L_VBF_STEP_4 = "abs(jet(0,30,4.7).eta-jet(1,30,4.7).eta) > 3 && mjj(0,1,30,4.7) > 300"
#### CUTS RELATIVE TO CONTROL REGION ONLY
MUID_LOOSE_CR = " && ".join(["pt > %f" % MU_PT_MIN, MU_PRESELECTION, SINGLE_SIP_CUT])
ELID_LOOSE_CR = " && ".join(["pt > %f" % EL_PT_MIN, EL_PRESELECTION, SINGLE_SIP_CUT])

SEL_ANY_CR_Z = "mass > 4"
