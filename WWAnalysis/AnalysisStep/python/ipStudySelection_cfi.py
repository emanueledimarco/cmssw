import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *

eleMatch = selectedPatElectrons.clone()
eleMatch.src = "patElectronsWithTrigger"
eleMatch.filter = cms.bool(False)
#eleMatch.cut = ( "genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==11 && abs(genParticleRef.get().mother().mother().pdgId()) ==24  ")
eleMatch.cut = ( "pt>10.  ")

muMatch = selectedPatMuons.clone()
muMatch.src = "patMuonsWithTrigger"
muMatch.filter = cms.bool(False)
#muMatch.cut = ( "pt>10. && genParticleRef.isNonnull() && abs(genParticleRef.get().pdgId())==13 &&  abs(genParticleRef.get().mother().mother().pdgId()) ==24 ")
muMatch.cut = ( "pt>10.")

# ----
boostedElectrons = cms.EDProducer("PatElectronBooster",
  electronTag = cms.untracked.InputTag("eleMatch"),
  vertexTag = cms.InputTag("goodPrimaryVertices"),
  electronIDSources = cms.PSet(),
  deposits = cms.VPSet(),
)

boostedMuons = cms.EDProducer("PatMuonBooster",
  muonTag = cms.untracked.InputTag("muMatch"),
  vertexTag = cms.InputTag("goodPrimaryVertices"),
  deposits = cms.VPSet(),
)


countBoostedMuons  = cms.EDFilter("CandViewCountFilter", 
   src = cms.InputTag("boostedMuons"), 
   minNumber = cms.uint32(1)
)

countBoostedElectrons  = cms.EDFilter("CandViewCountFilter", 
   src = cms.InputTag("boostedElectrons"), 
   minNumber = cms.uint32(1)
)

# ----
ELE_ID_CUT=("(" +
            "((isEB && sigmaIetaIeta < 0.01 &&" +           
            " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
            " hadronicOverEm < 0.04) ||" +
            " ( (!isEB) && sigmaIetaIeta < 0.03  &&  " +
            " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
            " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 &&" +
            " hadronicOverEm < 0.025 )) )" ) 
#            " abs(userFloat('dxyPV')) < 0.02 && " +
#            " abs(userFloat('dzPV'))  < 1.0  )  " )
#            " (abs(userFloat('convValueMapProd:dist')) > 0.02   ||" +
#            "  abs(userFloat('convValueMapProd:dcot')) > 0.02 ) )")
#            " userInt('expectedHitsEle') == 0 && " + 

MU_ID_CUT=("(isGlobalMuon && isTrackerMuon &&" +
           " innerTrack.found >10 &&" +
           " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
           " globalTrack.normalizedChi2 <10 &&" +
           " globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
           " numberOfMatches > 1 && " + 
           " abs(track.ptError / pt) < 0.10)")


# ----
ELE_ISO_CUT=("(( isEB && (dr03TkSumPt +" +
             " max(0,dr03EcalRecHitSumEt - 1.0) +" +
             " dr03HcalTowerSumEt)/pt < 0.10) ||"+
             "( isEE && (dr03TkSumPt +" +
             " dr03EcalRecHitSumEt +" +
             " dr03HcalTowerSumEt)/pt < 0.10) ) ")

MU_ISO_CUT=("(isolationR03().emEt +" +
            " isolationR03().hadEt +" +
            " isolationR03().sumPt)/pt < 0.15 ");


ELE_CONV_CUT=("( abs(userFloat('convValueMapProd:dist')) > 0.02   ||"+
              "  abs(userFloat('convValueMapProd:dcot')) > 0.02 ) &&" +
              "userInt('expectedHitsEle:in') == 0")


# ----
eleID = selectedPatElectrons.clone()
eleID.src = "boostedElectrons"
eleID.filter = cms.bool(False)
eleID.cut = ( "pt > 10 && " + ELE_ID_CUT)

muID = selectedPatMuons.clone()
muID.src = "boostedMuons"
muID.filter = cms.bool(False)
muID.cut = ( "pt > 10 && " + MU_ID_CUT)



# ----
eleISO = selectedPatElectrons.clone()
eleISO.src = "eleID"
eleISO.filter = cms.bool(False)
eleISO.cut = ( ELE_ISO_CUT)

muISO = selectedPatMuons.clone()
muISO.src = "muID"
muISO.filter = cms.bool(False)
muISO.cut = ( MU_ISO_CUT)



# ----
eleNoConv = selectedPatElectrons.clone()
eleNoConv.src = "eleISO"
eleNoConv.filter = cms.bool(False)
eleNoConv.cut = ( ELE_CONV_CUT )
                    

                    



ipStudyEleSelectionSequence = cms.Sequence(
    eleMatch*boostedElectrons*countBoostedElectrons*eleID*eleISO*eleNoConv
)

ipStudyMuSelectionSequence = cms.Sequence(
    muMatch*boostedMuons*countBoostedMuons*muID*muISO
)



