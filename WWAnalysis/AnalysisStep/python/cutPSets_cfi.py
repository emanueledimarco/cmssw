import FWCore.ParameterSet.Config as cms
import math

defaultWW = cms.VPSet(
    cms.PSet(
        label = cms.string("Skim"),
        cut = cms.string("(1)"),
        simple = cms.string("skim"),
    ),
    cms.PSet(
        label = cms.string("Trigger Path"),
        cut = cms.string("triggerBitsCut('MC')"),
#         cut = cms.string("(1)"),
        simple = cms.string("triggerPath"),
    ),
    cms.PSet(
        label = cms.string("Lepton Selection"),
        cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMin() > 10 && ptMax() > 20"),
        simple = cms.string("fiducial"),
    ),
    cms.PSet(
        label = cms.string("Lepton Veto"),
        cut = cms.string("nExtraLep(10) == 0"),
        simple = cms.string("extraLep"),
    ),
#     cms.PSet(
#         label = cms.string("Trigger Match"),
#         cut = cms.string("triggerMatchingCut('MC')"),
#         simple = cms.string("triggerMatch"),
#     ),
    cms.PSet(
        label = cms.string("#slash{E}_{T}"),
        cut = cms.string("pfMet() > 20"),
        simple = cms.string("pfMet"),
    ),
    cms.PSet(
        label = cms.string("m_{ll} > 12"),
        cut = cms.string("mll() > 12"),
        simple = cms.string("mll"),
    ),
    cms.PSet(
        label = cms.string("Z-veto"),
        cut = cms.string("abs(mll()-91.1876)>15.0") ,
        simple = cms.string("zVeto"),
    ),
    cms.PSet(
        label = cms.string("minProj#slash{E}_{T}"),
        cut = cms.string("min(projPfMet(),projChargedMet()) > 35"),
        simple = cms.string("minPMet"),
    ),
    cms.PSet(
        label = cms.string("Jet Veto"),
        cut = cms.string("nCentralJets( 30.0, 5.0, 1 ) == 0"),
        simple = cms.string("jetVeto"),
    ),
    cms.PSet(
        label = cms.string("Soft #mu Veto"),
        cut = cms.string("nSoftMu( 3.0 ) == 0"),
        simple = cms.string("softMu"),
    ),
    cms.PSet(
        label = cms.string("WW Selection"),
        cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) == 0"),
        simple = cms.string("bTag"),
    ),
)

oneJetWW = cms.VPSet(
    cms.PSet(
        label = cms.string("Skim"),
        cut = cms.string("(1)"),
        simple = cms.string("skim"),
    ),
    cms.PSet(
        label = cms.string("Trigger Path"),
        cut = cms.string("triggerBitsCut('MC')"),
#         cut = cms.string("(1)"),
        simple = cms.string("triggerPath"),
    ),
    cms.PSet(
        label = cms.string("Lepton Selection"),
        cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMin() > 10 && ptMax() > 20"),
        simple = cms.string("fiducial"),
    ),
    cms.PSet(
        label = cms.string("Lepton Veto"),
        cut = cms.string("nExtraLep(10) == 0"),
        simple = cms.string("extraLep"),
    ),
#     cms.PSet(
#         label = cms.string("Trigger Match"),
#         cut = cms.string("triggerMatchingCut('MC')"),
#         simple = cms.string("triggerMatch"),
#     ),
    cms.PSet(
        label = cms.string("#slash{E}_{T}"),
        cut = cms.string("pfMet() > 20"),
        simple = cms.string("pfMet"),
    ),
    cms.PSet(
        label = cms.string("m_{ll} > 12"),
        cut = cms.string("mll() > 12"),
        simple = cms.string("mll"),
    ),
    cms.PSet(
        label = cms.string("Z-veto"),
        cut = cms.string("abs(mll()-91.1876)>15.0") ,
        simple = cms.string("zVeto"),
    ),
    cms.PSet(
        label = cms.string("minProj#slash{E}_{T}"),
        cut = cms.string("min(projPfMet(),projChargedMet()) > 35"),
        simple = cms.string("minPMet"),
    ),
    cms.PSet(
        label = cms.string("One Jet"),
        cut = cms.string("nCentralJets( 30.0, 5.0, 1 ) == 1"),
        simple = cms.string("oneJet"),
    ),
    cms.PSet(
        label = cms.string("Jet Phi"),
        cut = cms.string("dPhiJetllInDegrees(0) < 165."),
        simple = cms.string("jetPhi"),
    ),
    cms.PSet(
        label = cms.string("High b-Jet Veto"),
        cut = cms.string("bTaggedJetsOver( 30.0, 2.1 ) == 0"),
        simple = cms.string("bTagOver"),
    ),
    cms.PSet(
        label = cms.string("Soft #mu Veto"),
        cut = cms.string("nSoftMu( 3.0 ) == 0"),
        simple = cms.string("softMu"),
    ),
    cms.PSet(
        label = cms.string("b-tag Veto"),
        cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) == 0"),
        simple = cms.string("bTagUnder"),
    ),
)

ttBar = cms.VPSet(
    cms.PSet(
        label = cms.string("Skim"),
        cut = cms.string("(1)"),
        simple = cms.string("skim"),
    ),
    cms.PSet(
        label = cms.string("Trigger Path"),
        cut = cms.string("triggerBitsCut('MC')"),
#         cut = cms.string("(1)"),
        simple = cms.string("triggerPath"),
    ),
    cms.PSet(
        label = cms.string("Lepton Selection"),
        cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMin() > 10 && ptMax() > 20"),
        simple = cms.string("fiducial"),
    ),
    cms.PSet(
        label = cms.string("Lepton Veto"),
        cut = cms.string("nExtraLep(10) == 0"),
        simple = cms.string("extraLep"),
    ),
    cms.PSet(
        label = cms.string("Trigger"),
        cut = cms.string("triggerMatchingCut('MC')"),
        simple = cms.string("trigger"),
    ),
    cms.PSet(
        label = cms.string("#slash{E}_{T}"),
        cut = cms.string("pfMet() > 20"),
        simple = cms.string("pfMet"),
    ),
    cms.PSet(
        label = cms.string("m_{ll} > 12"),
        cut = cms.string("mll() > 12"),
        simple = cms.string("mll"),
    ),
    cms.PSet(
        label = cms.string("Z-veto"),
        cut = cms.string("abs(mll()-91.1876)>15.0") ,
        simple = cms.string("zVeto"),
    ),
    cms.PSet(
        label = cms.string("pMin#slash{E}_{T}"),
        cut = cms.string("min(projPfMet(),projChargedMet()) > 35"),
        simple = cms.string("minPMet"),
    ),
)

highPtB = cms.PSet(
    label = cms.string("b-Jet"),
    cut = cms.string("bTaggedJetsOver( 30.0, 2.1 ) == 1 && nCentralJets( 30.0, 5.0, 1) == 1"),
    simple = cms.string("oneBOver"),
)

softMu = cms.PSet(
    label = cms.string("Soft #mu Veto"),
    cut = cms.string("nSoftMu( 3.0 ) > 0"),
    simple = cms.string("hasSoftMu"),
)

bTag = cms.PSet(
    label = cms.string("b-tag Veto"),
    cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) > 0"),
    simple = cms.string("bUnder"),
)

NoBtag = cms.PSet(
    label = cms.string("b-tag Veto"),
    cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) == 0"),
    simple = cms.string("noBTag"),
)

topTag = cms.PSet(
    label = cms.string("Top Veto"),
    cut = cms.string("bTaggedJetsUnder( 30.0,  2.1 ) > 0 || nSoftMu( 3.0 ) > 0"),
    simple = cms.string("hasTopTag"),
)


def cloneVPSet(oldvpset):
    newvpset = cms.VPSet()
    for pset in oldvpset:
        newvpset.append(pset.clone())
    return newvpset

def switchToOppoFlavor(vp):
    vp[7].cut = "(1)"
    vp[8].cut = "min(projPfMet(),projChargedMet()) > 20.0"

def switchToOppoFlavor1Jet(vp):
    vp[7].cut = "(1)"
    vp[8].cut = "min(projPfMet(),projChargedMet()) > 20.0"
    vp[10].cut = "(1)"

def switchToTrailingElectron(vp):
    vp[2].cut = cms.string("q(0)*q(1) < 0 && !isSTA(0) && !isSTA(1) && leptEtaCut(2.4,2.5) && ptMin() > 15 && ptMax() > 20")

#swaps for PDs
MC =             "MC"
SingleMuon =     "SingleMuon"
DoubleMuon =     "DoubleMuon"
MuEG =           "MuEG"
DoubleElectron = "DoubleElectron"

def switchTrigger(vp,value):
    vp[1].cut = "triggerBitsCut('{0}')".format(value)
    vp[4].cut = "triggerMatchingCut('{0}')".format(value)

def addMassDependentCuts(vpa,vpb):
    for pset in vpb:
        vpa.append(pset)

def injectTightIsolationCut(vp):
    vp.insert( 2, cms.PSet(
        label = cms.string("Isolation"),
        cut   = cms.string("(ptMin > 20 || hypo == 3 || hypo == 5 || allIsoByPt(1) < 0.10)" ),
        simple = cms.string("isolation"),
    ))

def injectIPCut(vp):
    vp.insert( 2, cms.PSet(
        label = cms.string("IP"),
        cut = cms.string("( passesIP() )"),
        simple = cms.string("ip"),
    ))

