import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *

MUON_ID_CUT=("(( (isGlobalMuon() && "
             "    globalTrack.normalizedChi2 <10 &&" +
             "    globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
             "    numberOfMatches > 1 ) || " + 
             "   (isTrackerMuon() && muonID('TMLastStationTight')) ) && " + 
             " innerTrack.found >10 &&" +
             " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
             " abs(track.ptError / pt) < 0.10 )")

MUON_KINKID_CUT=("(( (isGlobalMuon() && "
                 "    globalTrack.normalizedChi2 <10 &&" +
                 "    globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
                 "    numberOfMatches > 1 ) || " + 
                 "   (isTrackerMuon() && muonID('TMLastStationTight')) ) && " + 
                 " combinedQuality.trkKink < 20 &&" +
                 " innerTrack.found >10 &&" +
                 " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
                 " abs(track.ptError / pt) < 0.10 )")

MUON_ID_CUT_OLD=("(isGlobalMuon && isTrackerMuon &&" +
                 " innerTrack.found >10 &&" +
                 " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
                 " globalTrack.normalizedChi2 <10 &&" +
                 " globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
                 " numberOfMatches > 1 && " + 
                 " abs(track.ptError / pt) < 0.10 )")

#SMURF_ISO = ("( userFloat('smurfCharged') + userFloat('smurfPhoton') + userFloat('smurfNeutral') )/ pt")
SMURF_ISO = ("( userFloat('muSmurfPF') )/ pt")
MUON_MERGE_ISO  =   ("( (abs(eta) < 1.479 && pt >  20 && " + SMURF_ISO + " < 0.13) || ( abs(eta) >= 1.479 && pt >  20 && " + SMURF_ISO + " < 0.09 ) || " + 
                     "  (abs(eta) < 1.479 && pt <= 20 && " + SMURF_ISO + " < 0.06) || ( abs(eta) >= 1.479 && pt <= 20 && " + SMURF_ISO + " < 0.05 ) )  ")

MUON_MERGE_IP  = ("( ( (pt >= 20 && abs(userFloat('tip')) < 0.02) || (pt < 20 && abs(userFloat('tip')) < 0.01) ) && " +
                  "  abs(userFloat('dzPV'))  < 0.1 )" )
              
              

#MUON_ID_LOOSE = ( "( " + MUON_ID_CUT + " &&" +
MUON_ID_LOOSE = ( "( " + MUON_KINKID_CUT + " &&" +
                  "  ( userFloat('muSmurfPF') < 0.4 * pt ) && abs(userFloat('tip')) < 0.2 && abs(userFloat('dzPV')) < 0.1 )" )


MUON_ISO_CUT = ("(isolationR03().emEt +" +
                " isolationR03().hadEt +" +
                " isolationR03().sumPt - userFloat('rhoMu')*3.14159265*0.3*0.3)/pt < 0.15 ");

MUON_ISO_CUT_TIGHT=("( ( pt > 20 && (isolationR03().emEt + isolationR03().hadEt + " +
                        " isolationR03().sumPt - userFloat('rhoMu')*3.14159265*0.3*0.3)/pt < 0.15 ) || " + 
                    "  ( pt <= 20 && (isolationR03().emEt + isolationR03().hadEt +" +
                        " isolationR03().sumPt - userFloat('rhoMu')*3.14159265*0.3*0.3)/pt < 0.12 ) )");

MUON_ISOPF_CUT = ("( (userFloat('pfCharged')+userFloat('pfPhoton')+userFloat('pfNeutral')-userFloat('rhoMuNoPU')*3.14159265*0.4*0.4) / pt < 0.20)")


MUON_IP_CUT=( "( abs(userFloat('tip2')) < 0.01 && " +
              "  abs(userFloat('dzPV'))  < 0.05    )" )


MUON_ID_CUT_4VETO=("(isTrackerMuon &&" +
                   " muonID('TMLastStationAngTight') &&" +
                   " innerTrack.found >10 && abs(userFloat('tip')) < 0.2 && abs(userFloat('dzPV')) < 0.1 &&" +
                   " ( (pt <= 20) || " +
                   "   (pt >  20  && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt" +
                   "                 )/pt > 0.10) ) )")

MUON_ID_CUT_4VETO_NEW =("(isTrackerMuon &&" +
                        " muonID('TMLastStationAngTight') &&" +
                        " innerTrack.found >10 && abs(userFloat('tip')) < 0.2 && abs(userFloat('dzPV')) < 0.2 &&" +
                        " ( (pt <= 20) || (pt >  20  && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt)/pt > 0.10) ) )")

PRESEL_MU=("pt > 10 && abs(eta)<2.4")

selectedRefPatMuons = cms.EDFilter("PATMuonRefSelector",
   src = cms.InputTag("boostedMuons"),
   filter = cms.bool(False),
   cut = cms.string("")
)

wwMuMatch = selectedRefPatMuons.clone()
wwMuMatch.cut = (PRESEL_MU)
#wwMuMatch.cut = ( "pt > 10 && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==13 && abs(genParticleRef.get().mother().mother().pdgId()) ==24")


wwMuonsLoose = selectedRefPatMuons.clone()
wwMuonsLoose.cut = (PRESEL_MU +"&&"+ MUON_ID_LOOSE )

wwMuonsID = selectedRefPatMuons.clone()
wwMuonsID.cut = (PRESEL_MU +"&&"+ MUON_ID_CUT )

wwMuonsMergeID = selectedRefPatMuons.clone()
wwMuonsMergeID.cut = (PRESEL_MU +"&&"+ MUON_ID_CUT )

wwMuonsISO = selectedRefPatMuons.clone()
wwMuonsISO.cut = (PRESEL_MU +"&&"+ MUON_ID_CUT +"&&"+ MUON_ISO_CUT )

wwMuonsISOT = selectedRefPatMuons.clone()
wwMuonsISOT.cut = (PRESEL_MU +"&&"+ MUON_ID_CUT +"&&"+ MUON_ISO_CUT_TIGHT )

wwMuonsMergeISO = selectedRefPatMuons.clone()
wwMuonsMergeISO.cut = (PRESEL_MU +"&&"+ MUON_ID_CUT +"&&"+ MUON_MERGE_ISO )

wwMuonsISOPF = selectedRefPatMuons.clone()
wwMuonsISOPF.cut = (PRESEL_MU +"&&"+  MUON_ID_CUT +"&&"+ MUON_ISOPF_CUT )

wwMuonsIP = selectedRefPatMuons.clone()
wwMuonsIP.cut = ( PRESEL_MU +"&&"+  MUON_ID_CUT +"&&"+ MUON_IP_CUT )

#wwMuonsIPT = selectedRefPatMuons.clone()
#wwMuonsIPT.src = "wwMuonsISOT"
#wwMuonsIPT.filter = cms.bool(False)
#wwMuonsIPT.cut = ( MUON_IP_CUT )

wwMuonsMergeCONV = selectedRefPatMuons.clone()
wwMuonsMergeCONV.cut = (PRESEL_MU +"&&"+ MUON_ID_CUT +"&&"+ MUON_MERGE_ISO )


wwMuonsMergeIP = selectedRefPatMuons.clone()
wwMuonsMergeIP.cut = (PRESEL_MU +"&&"+  MUON_ID_CUT +"&&"+ MUON_MERGE_ISO+"&&"+MUON_MERGE_IP )

#wwMuonsIPPF = selectedRefPatMuons.clone()
#wwMuonsIPPF.src = "wwMuonsISOPF"
#wwMuonsIPPF.filter = cms.bool(False)
#wwMuonsIPPF.cut = ( MUON_IP_CUT )

wwMuons4Veto = selectedRefPatMuons.clone()
wwMuons4Veto.cut = ( "pt > 3 && " +
                     MUON_ID_CUT_4VETO )

wwMuons4VetoNew = selectedRefPatMuons.clone()
wwMuons4VetoNew.cut = ( "pt > 3 && " +
                         MUON_ID_CUT_4VETO_NEW )

#       _                 _                         
#      | |               | |                        
#      | | __ _ _ __ ___ | |__   ___  _ __ ___  ___ 
#  _   | |/ _` | '_ ` _ \| '_ \ / _ \| '__/ _ \/ _ \
# | |__| | (_| | | | | | | |_) | (_) | | |  __/  __/
#  \____/ \__,_|_| |_| |_|_.__/ \___/|_|  \___|\___|
#                                                   
#                                                   
#   _____                           _           
#  / ____|                         (_)          
# | (___   ___ ___ _ __   __ _ _ __ _  ___  ___ 
#  \___ \ / __/ _ \ '_ \ / _` | '__| |/ _ \/ __|
#  ____) | (_|  __/ | | | (_| | |  | | (_) \__ \
# |_____/ \___\___|_| |_|\__,_|_|  |_|\___/|___/
#                                               

# Scenario 1 is LP selection:
Scenario1_LP_MUONS    =  PRESEL_MU +"&&"+  MUON_ID_CUT +"&&"+ MUON_MERGE_ISO+"&&"+MUON_MERGE_IP 
wwMuScenario1 = selectedRefPatMuons.clone( cut = Scenario1_LP_MUONS )
wwMu4VetoScenario1 = wwMuons4Veto.clone()

# Scenario 2 is LP selection + muon kink finder, no change in the electrons
Scenario2_KINK_MUONS  =  PRESEL_MU +"&&"+  MUON_KINKID_CUT +"&&"+ MUON_MERGE_ISO+"&&"+MUON_MERGE_IP 
wwMuScenario2 = selectedRefPatMuons.clone( cut = Scenario2_KINK_MUONS )
wwMu4VetoScenario2 = wwMuons4VetoNew.clone()

# Scenario 3 is LP selection + muon kink finder switching to new LH ele ID
# the loose ID is applied before the LH ID 
Scenario3_KINK_MUONS  =  PRESEL_MU +"&&"+  MUON_KINKID_CUT +"&&"+ MUON_MERGE_ISO+"&&"+MUON_MERGE_IP 
wwMuScenario3 = selectedRefPatMuons.clone( cut = Scenario3_KINK_MUONS )
wwMu4VetoScenario3 = wwMuons4VetoNew.clone()

# Scenario 4 is LP selection + muon kink finder switching to new BDT ele ID 
# also have to change the source because boris is going to add a new filter to make the new BDT passing electrons
# alhtough they still need the loose cuts applied
Scenario4_KINK_MUONS  =  PRESEL_MU +"&&"+  MUON_KINKID_CUT +"&&"+ MUON_MERGE_ISO+"&&"+MUON_MERGE_IP 
wwMuScenario4 = selectedRefPatMuons.clone( cut = Scenario4_KINK_MUONS )
wwMu4VetoScenario4 = wwMuons4VetoNew.clone()

# Scenario 5 is a loose loose tree to do the Wjets studies
Scenario5_LOOSE_MUONS =  MUON_ID_LOOSE
wwMuScenario5 = selectedRefPatMuons.clone( cut = Scenario5_LOOSE_MUONS )
wwMu4VetoScenario5 = wwMuons4VetoNew.clone()

# 
#  ______           _ 
# |  ____|         | |
# | |__   _ __   __| |
# |  __| | '_ \ / _` |
# | |____| | | | (_| |
# |______|_| |_|\__,_|
#                     
#                     


wwMuonSequence = cms.Sequence( 
    wwMuMatch * 
    wwMuonsID * 
    wwMuonsISO * 
    wwMuonsIP * 
    wwMuonsMergeID * 
    wwMuonsMergeISO * 
    wwMuonsMergeCONV * 
    wwMuonsMergeIP * 
#     wwMuonsISOT * 
#     wwMuonsIPT * 
#     wwMuonsISOPF * 
#     wwMuonsIPPF * 
    wwMuons4Veto
)

