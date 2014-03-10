// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class HWWTreeDumper
//      Analyzer module that takes the Candidate Collections from
//      the analysis producers and dumps an ntuple
//      
//-----------------------------------------------------------------------



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMuonFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFTauFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFPreIdFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPhotonFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConversionFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCalibElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFlowElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsSuperClusterFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsBasicClusterFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGenInfoFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTrackFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGsfTrackFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsVertexFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJPTJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCaloTowerFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEcalRecHitFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsV0CandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthExtraTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsRunInfoFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsHcalNoiseFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPdfWeightFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"


HWWTreeDumper::HWWTreeDumper(const edm::ParameterSet& iConfig)
{
  
  nameFile_      = iConfig.getUntrackedParameter<std::string>("nameFile", "RootOutput.root");
  nameTree_      = iConfig.getUntrackedParameter<std::string>("nameTree", "BaseTree");
  dumpTree_      = iConfig.getUntrackedParameter<bool>("dumpTree", false);
  dumpMCTruth_   = iConfig.getUntrackedParameter<bool>("dumpMCTruth", false);
  dumpMCTruthExtra_ = iConfig.getUntrackedParameter<bool>("dumpMCTruthExtra", false);
  is8TeV_        = iConfig.getUntrackedParameter<bool>("is8TeV", true);
  
  // control level of Reco Adapters in the tree
  saveTrk_        = iConfig.getUntrackedParameter<bool>("saveTrk", false);
  saveEcal_       = iConfig.getUntrackedParameter<bool>("saveEcal", false);
  saveHcal_       = iConfig.getUntrackedParameter<bool>("saveHcal", false);
  saveDT_         = iConfig.getUntrackedParameter<bool>("saveDT", false);
  saveCSC_        = iConfig.getUntrackedParameter<bool>("saveCSC", false);
  saveRPC_        = iConfig.getUntrackedParameter<bool>("saveRPC", false);
  saveFatTrk_     = iConfig.getUntrackedParameter<bool>("saveFatTrk", false);
  saveTrackDeDx_  = iConfig.getUntrackedParameter<bool>("saveTrackDeDx", false);
  saveTrackImpactParameters_ = iConfig.getUntrackedParameter<bool>("saveTrackImpactParameters",true);
  saveFatEcal_    = iConfig.getUntrackedParameter<bool>("saveFatEcal", false);
  saveFatHcal_    = iConfig.getUntrackedParameter<bool>("saveFatHcal", false);
  saveFatDT_      = iConfig.getUntrackedParameter<bool>("saveFatDT", false);
  saveFatCSC_     = iConfig.getUntrackedParameter<bool>("saveFatCSC", false);
  saveFatRPC_     = iConfig.getUntrackedParameter<bool>("saveFatRPC", false);
  saveJetBTag_    = iConfig.getUntrackedParameter<bool>("saveJetBTag", false);

  //electron pflow
  savePFEleBasic_  = iConfig.getUntrackedParameter<bool>("savePFEleBasic",  true);
  savePFEleIsoDep_ = iConfig.getUntrackedParameter<bool>("savePFEleIsoDep", true);

  // particle identification
  saveEleID_    = iConfig.getUntrackedParameter<bool>("saveEleID", false);

  // basic kinematic informations
  saveCand_  = iConfig.getUntrackedParameter<bool>("saveCand", true);

  // Candidate Collections
  dumpPreselInfo_     = iConfig.getUntrackedParameter<bool>("dumpPreselInfo", false);
  dumpGenInfo_        = iConfig.getUntrackedParameter<bool>("dumpGenInfo", false);
  dumpLHE_            = iConfig.getUntrackedParameter<bool>("dumpLHE", false);
  dumpElectrons_      = iConfig.getUntrackedParameter<bool>("dumpElectrons", false);
  dumpCalibratedElectrons_ = iConfig.getUntrackedParameter<bool>("dumpCalibratedElectrons", false);
  dumpPhotons_        = iConfig.getUntrackedParameter<bool>("dumpPhotons", false);
  dumpConversions_    = iConfig.getUntrackedParameter<bool>("dumpConversions", false);
  dumpPFlowElectrons_ = iConfig.getUntrackedParameter<bool>("dumpPFlowElectrons", false);
  dumpSCs_            = iConfig.getUntrackedParameter<bool>("dumpSCs", false);
  dumpBCs_            = iConfig.getUntrackedParameter<bool>("dumpBCs", false);
  dumpPhoPFClusters_  = iConfig.getUntrackedParameter<bool>("dumpPhoPFClusters", false);
  dumpLinkedTracks_   = iConfig.getUntrackedParameter<bool>("dumpLinkedTracks", false);
  dumpTracks_         = iConfig.getUntrackedParameter<bool>("dumpTracks", false);
  dumpGsfTracks_      = iConfig.getUntrackedParameter<bool>("dumpGsfTracks", false);
  dumpMuonTracks_     = iConfig.getUntrackedParameter<bool>("dumpMuonTracks", false);
  dumpMuons_          = iConfig.getUntrackedParameter<bool>("dumpMuons", false);
  dumpJets_           = iConfig.getUntrackedParameter<bool>("dumpJets", false);
  dumpGenJets_        = iConfig.getUntrackedParameter<bool>("dumpGenJets", false);
  dumpPUcorrPFJet_    = iConfig.getUntrackedParameter<bool>("dumpPUcorrPFJet", false);
  dumpCHSPFJet_       = iConfig.getUntrackedParameter<bool>("dumpCHSPFJet", false);
  dumpMet_            = iConfig.getUntrackedParameter<bool>("dumpMet", false);
  dumpGenMet_         = iConfig.getUntrackedParameter<bool>("dumpGenMet", false);
  dumpVertices_       = iConfig.getUntrackedParameter<bool>("dumpVertices", false);
  dumpK0s_            = iConfig.getUntrackedParameter<bool>("dumpK0s", false);
  dumpCaloTowers_     = iConfig.getUntrackedParameter<bool>("dumpCaloTowers", false);
  dumpEcalRecHits_    = iConfig.getUntrackedParameter<bool>("dumpEcalRecHits", false);
  dumpLogErrorFlags_  = iConfig.getUntrackedParameter<bool>("dumpLogErrorFlags", false); 
  dumpHcalNoiseFlags_ = iConfig.getUntrackedParameter<bool>("dumpHcalNoiseFlags", false);
  aodHcalNoiseFlags_  = iConfig.getUntrackedParameter<bool>("AODHcalNoiseFlags", true);

  // Particle Flow objects
  dumpParticleFlowObjects_ = iConfig.getUntrackedParameter<bool>("dumpParticleFlowObjects",false);
  dumpPFpreId_             = iConfig.getUntrackedParameter<bool>("dumpPFpreId", false);  

  // data run informations
  dumpRunInfo_ = iConfig.getUntrackedParameter<bool>("dumpRunInfo",false);

  // eigenvalues for PDF systematics
  dumpPdfWeight_ = iConfig.getUntrackedParameter<bool>("dumpPdfWeight",false);

  electronCollection_      = iConfig.getParameter<edm::InputTag>("electronCollection");
  calibElectronCollection_ = iConfig.getParameter<edm::InputTag>("calibElectronCollection");
  pflowElectronCollection_ = iConfig.getParameter<edm::InputTag>("pflowElectronCollection");
  photonCollection_        = iConfig.getParameter<edm::InputTag>("photonCollection");
  muonCollection_          = iConfig.getParameter<edm::InputTag>("muonCollection");
  calibMuonCollection_     = iConfig.getParameter<edm::InputTag>("calibMuonCollection");
  PFCandidateCollection_   = iConfig.getParameter<edm::InputTag>("PFCandidateCollection");
  PFPUCandidateCollection_  = iConfig.getParameter<edm::InputTag>("PFPUCandidateCollection");
  PFNoPUCandidateCollection_  = iConfig.getParameter<edm::InputTag>("PFNoPUCandidateCollection");
  ecalSCCollection_        = iConfig.getParameter<edm::InputTag>("ecalSCCollection");
  ecalElePFClusterCollection_ = iConfig.getParameter<edm::InputTag>("ecalElePFClusterCollection");
  ecalPhoPFClusterCollection_ = iConfig.getParameter<edm::InputTag>("ecalPhoPFClusterCollection");
  ecalBCCollection_        = iConfig.getParameter<edm::InputTag>("ecalBCCollection");
  ecalBarrelRecHits_       = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHits");
  ecalEndcapRecHits_       = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHits");
  esRecHits_               = iConfig.getParameter<edm::InputTag>("esRecHits");
  calotowersForIsolationProducer_ = iConfig.getParameter<edm::InputTag>("calotowersForIsolationProducer");
  conversions_             = iConfig.getParameter<edm::InputTag>("conversionCollection");
  trackCollection_         = iConfig.getParameter<edm::InputTag>("trackCollection");
  generalTrackCollection_  = iConfig.getParameter<edm::InputTag>("generalTrackCollection");
  refittedForDeDxTrackCollection_ = iConfig.getParameter<edm::InputTag>("refittedForDeDxTrackCollection");
  gsfTrackCollection_      = iConfig.getParameter<edm::InputTag>("gsfTrackCollection");
  globalMuonTrackCollection_ = iConfig.getParameter<edm::InputTag>("globalMuonTrackCollection");
  standAloneMuonTrackCollection_ = iConfig.getParameter<edm::InputTag>("standAloneMuonTrackCollection");
  vertexCollection_        = iConfig.getParameter<edm::InputTag>("vertexCollection");
  K0sCollection_           = iConfig.getParameter<edm::InputTag>("K0sCollection");
  genJetCollection_        = iConfig.getParameter<edm::InputTag>("genJetCollection");
  jetCollection1_          = iConfig.getParameter<edm::InputTag>("jetCollection1");
  jetCollection2_          = iConfig.getParameter<edm::InputTag>("jetCollection2");
  PFjetCollection1_        = iConfig.getParameter<edm::InputTag>("PFjetCollection1");
  PFpuCorrJetCollection1_  = iConfig.getParameter<edm::InputTag>("PFpuCorrJetCollection1");
  PFJetCorrectionService_  = iConfig.getParameter<std::string>("PFJetCorrectionService");
  JetCorrectionService_    = iConfig.getParameter<std::string>("JetCorrectionService");
  JPTjetCollection1_       = iConfig.getParameter<edm::InputTag>("JPTjetCollection1");
  JPTjetCollection2_       = iConfig.getParameter<edm::InputTag>("JPTjetCollection2");

  // btag collections
  PFJetsBTags_              = iConfig.getUntrackedParameter<edm::ParameterSet>("PFJetsBTags");
  PFPUcorrJetsBTags_        = iConfig.getUntrackedParameter<edm::ParameterSet>("PFPUcorrJetsBTags");

  // MVA based jet id collection
  puJetIDAlgos_            = iConfig.getParameter<std::vector<edm::ParameterSet> >("puJetIDAlgos");
  metCollection_           = iConfig.getParameter<edm::InputTag>("metCollection");
  // corrmetCollection_       = iConfig.getParameter<edm::InputTag>("corrmetCollection");
  TCmetCollection_         = iConfig.getParameter<edm::InputTag>("TCmetCollection");
  PFmetCollection_         = iConfig.getParameter<edm::InputTag>("PFmetCollection");
  PFChMetCollection_       = iConfig.getParameter<edm::InputTag>("PFChMetCollection");                       
  leptonLinkedPFCandidates_ = iConfig.getParameter<edm::InputTag>("leptonLinkedPFCandidates");
  genMetCollection_        = iConfig.getParameter<edm::InputTag>("genMetCollection");
  chargedMetCollection_    = iConfig.getParameter<edm::InputTag>("chargedMetCollection");
  mcTruthCollection_       = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
  hepMcCollection_         = iConfig.getParameter<edm::InputTag>("hepMcCollection");
  genInfoCollection_       = iConfig.getParameter<edm::InputTag>("genInfoCollection");
  genWeightCollection_     = iConfig.getUntrackedParameter<std::string>("genWeightCollection");
  PFpreIdCollection_       = iConfig.getParameter<edm::InputTag>("PFpreIdCollection");

  // calotowers collections
  calotowerCollection_ = iConfig.getParameter<edm::InputTag>("calotowerCollection");
  hbheLabel_  = iConfig.getParameter<edm::InputTag>("hbheInput");
  hoLabel_    = iConfig.getParameter<edm::InputTag>("hoInput");
  hfLabel_    = iConfig.getParameter<edm::InputTag>("hfInput");
  ecalLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("ecalInputs");

  // ECAL rechits collections
  ebRHLabel_ = iConfig.getParameter<edm::InputTag>("EBRecHits");
  eeRHLabel_ = iConfig.getParameter<edm::InputTag>("EERecHits");

  // trigger Collections
  dumpTriggerResults_  = iConfig.getUntrackedParameter<bool>("dumpTriggerResults");
  dumpHLTObject_       = iConfig.getUntrackedParameter<bool>("dumpHLTObjects");
  hltParms_            = iConfig.getUntrackedParameter<edm::ParameterSet>("HLTObjectsInfo");

  // dump PFCandidates
  dumpPFCandidates_  = iConfig.getUntrackedParameter<bool>("dumpPFCandidates");

  // Hcal collections
  hcalNoiseSummaryLabel_ = iConfig.getParameter<edm::InputTag>("hcalNoiseSummary");

  // PDF sets
  pdfSet1_ = iConfig.getParameter<edm::InputTag>("pdfSet1");
  pdfSet2_ = iConfig.getParameter<edm::InputTag>("pdfSet2");
  pdfSet3_ = iConfig.getParameter<edm::InputTag>("pdfSet3");

  namePdf1_ = iConfig.getUntrackedParameter<std::string>("namepdf1", "pdfSet1");
  namePdf2_ = iConfig.getUntrackedParameter<std::string>("namepdf2", "pdfSet2");
  namePdf3_ = iConfig.getUntrackedParameter<std::string>("namepdf3", "pdfSet3");

}



HWWTreeDumper::~HWWTreeDumper() { 
  if(hltObjectFiller_) delete(hltObjectFiller_);
}



// ------------ method called to for each event  ------------
void HWWTreeDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  /// fill the run info (run number, event, ...)
  if(dumpRunInfo_) {
    CmsRunInfoFiller runFiller( tree_, dumpMCTruth_, is8TeV_ );
    runFiller.dumpL1Trigger(dumpTriggerResults_);
    runFiller.dumpLogErrorFlags(dumpLogErrorFlags_);
    runFiller.writeRunInfoToTree(iEvent,iSetup,false);
  }

  // get MC truth
  CmsMcTruthTreeFiller treeFill(tree_);
  treeFill.saveLHEComments(dumpLHE_);

  if(dumpMCTruth_) {

    bool firstEvent = (jevt_>1) ? false : true;
    treeFill.writeCollectionToTree( mcTruthCollection_, LHEComments_, iEvent, 100, firstEvent );
  }

  // get MC truth for tH analysis 
  CmsMcTruthExtraTreeFiller treeFillExtra(tree_);
  treeFillExtra.saveLHEComments(false); // we do not care for this now                                                                                                
  if(dumpMCTruthExtra_) {
    bool firstEvent = (jevt_>1) ? false : true;
    treeFillExtra.writeCollectionToTree( mcTruthCollection_, LHEComments_, iEvent, 100, firstEvent );
  }


  // fill the PDF weights
  if(dumpPdfWeight_) {
    CmsPdfWeightFiller pdfWeightFiller( tree_);
    pdfWeightFiller.writePdfWeightToTree(pdfSet1_, iEvent, iSetup, "", namePdf1_.c_str(), false);
    pdfWeightFiller.writePdfWeightToTree(pdfSet2_, iEvent, iSetup, "", namePdf2_.c_str(), false);
    pdfWeightFiller.writePdfWeightToTree(pdfSet3_, iEvent, iSetup, "", namePdf3_.c_str(), false);
  }

  jevtInRun_++;
  // fill the trigger paths info
  if(dumpTriggerResults_) {

    CmsTriggerTreeFiller triggerTreeFill (tree_, hltParms_) ;
    std::string prefix ("") ;
    std::string suffix ("Trg") ;
    triggerTreeFill.writeTriggerToTree (iEvent, prefix, suffix) ;

    /// fill the trigger mask in the Event tree
    bool firstEvent = true;
    if(jevt_>1) firstEvent=false; 
    CmsConditionsFiller conditionsFiller( tree_, hltParms_, trgNames_ );
    conditionsFiller.writeConditionsToTree(iEvent,  firstEvent);
    jevt_++;
  }
  if(dumpHLTObject_) {
    //forward beginRun to HLTObjectFiller
    if(jevtInRun_==1) hltObjectFiller_->beginRun( iEvent, iEvent.getRun(), iSetup);
    hltObjectFiller_->writeHLTObjectToTree(iEvent);
  }

  // fill preselection output
  if (dumpPreselInfo_) {

    Handle<bool> selected;
    iEvent.getByLabel("preselectionMarker", selected );
    bool isSelected = *selected;
    tree_->column ("evtPresel", isSelected, false, "Pres");

  }


  // fill Electrons block
  if(dumpElectrons_) {

    CmsElectronFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Ele");
    treeFill.saveCand(saveCand_);
    treeFill.saveTrk(saveTrk_);
    treeFill.saveEcal(saveEcal_);
    treeFill.saveFatTrk(saveFatTrk_);
    treeFill.saveFatEcal(saveFatEcal_);
    treeFill.setGeneralTracks(trackCollection_);
    // if set to false, it will not save any PF isolation-related variable (ALCARECO has not pflowcands)
    treeFill.savePFlowIsolations(dumpParticleFlowObjects_);
    treeFill.setEcalSuperClusters(ecalSCCollection_);
    treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
    treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
    // for custom isolation
    treeFill.setTracksProducer(generalTrackCollection_);
    treeFill.setCalotowersProducer(calotowersForIsolationProducer_);
    // for custom MVA id
    treeFill.setVertexCollection(vertexCollection_);
    treeFill.setEleIdMVAs(fEleIdMVATrig,fEleIdMVATrigIdIsoCombined,fEleIdMVANonTrig);
    treeFill.setPFCandidateCollection(PFNoPUCandidateCollection_);
    treeFill.saveEleID(true);
    // for full vertex fit conversion veto
    treeFill.setConversionsProdcer(conversions_);
    // for the fully calibrated electrons
    treeFill.setCalibElectronCollection(calibElectronCollection_);

    treeFill.writeCollectionToTree(electronCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  // fill minimal electrons block after energy corrections
  if(dumpCalibratedElectrons_) {
    CmsCalibElectronFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("CalibEle");
    treeFill.writeCollectionToTree(calibElectronCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  if(dumpPFlowElectrons_) {
    CmsPFlowElectronFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("PFEle");
    treeFill.saveCand(saveCand_);
    treeFill.setGeneralTracks(trackCollection_);
    treeFill.savePFEleBasic(savePFEleBasic_);
    treeFill.savePFEleIsoDep(savePFEleIsoDep_);
    treeFill.writeCollectionToTree(pflowElectronCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  if(dumpPFpreId_) {
    CmsPFPreIdFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("PFpreId");
    treeFill.writeCollectionToTree(PFpreIdCollection_, trackCollection_, iEvent, iSetup, prefix, suffix, false);  
  }

  // fill Photons block
  if(dumpPhotons_) {

    CmsPhotonFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Pho");
    treeFill.saveCand(saveCand_);
    treeFill.setEcalSuperClusters(ecalSCCollection_);
    treeFill.setPFCandidates(PFCandidateCollection_);
    treeFill.setPrimaryVertices(vertexCollection_);
    // for full vertex fit conversion veto
    treeFill.setConversionsProdcer(conversions_);
    treeFill.writeCollectionToTree(photonCollection_, iEvent, iSetup, prefix, suffix, false);
  }
  if(dumpConversions_) {

    CmsConversionFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Conv");
    treeFill.writeCollectionToTree(conversions_, iEvent, iSetup, prefix, suffix, false);
  }

  // fill SC block
  if (dumpSCs_) {

      CmsSuperClusterFiller treeFill(tree_, 1000);
      std::string prefix("");
      std::string suffix("SC");
      treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
      treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
      treeFill.setESRecHits(esRecHits_);
      treeFill.setCalotowers(calotowersForIsolationProducer_);
      treeFill.writeCollectionToTree(ecalSCCollection_, iEvent, iSetup, prefix, suffix, false);
      
      if(dumpParticleFlowObjects_) {
        CmsSuperClusterFiller treeFillPF(tree_, 1000);
        suffix = "PFSC";
        treeFillPF.setEcalBarrelRecHits(ecalBarrelRecHits_);
        treeFillPF.setEcalEndcapRecHits(ecalEndcapRecHits_);
        treeFillPF.setESRecHits(esRecHits_);
        treeFillPF.setCalotowers(calotowersForIsolationProducer_);
        treeFillPF.writeCollectionToTree(ecalElePFClusterCollection_, iEvent, iSetup, prefix, suffix, false);
        
        if(dumpPhoPFClusters_) {
          CmsSuperClusterFiller treeFillPhoPF(tree_, 1000);
          suffix = "PhoPFSC";
          treeFillPhoPF.setEcalBarrelRecHits(ecalBarrelRecHits_);
          treeFillPhoPF.setEcalEndcapRecHits(ecalEndcapRecHits_);
          treeFillPhoPF.setESRecHits(esRecHits_);
          treeFillPhoPF.setCalotowers(calotowersForIsolationProducer_);
          treeFillPhoPF.writeCollectionToTree(ecalPhoPFClusterCollection_, iEvent, iSetup, prefix, suffix, false);
        }
      }
  }

  // fill BC block
  if (dumpBCs_) {

      CmsBasicClusterFiller treeFill(tree_, 350);
      std::string prefix("");
      std::string barrelSuffix("BC");
      treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
      treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
      treeFill.setEcalSuperClusters(ecalSCCollection_);
      treeFill.writeCollectionToTree(ecalBCCollection_, iEvent, iSetup, prefix, barrelSuffix, false);

      if(dumpParticleFlowObjects_) {
        CmsBasicClusterFiller treeFillPF(tree_, 350);
        prefix="";
        barrelSuffix="PFBC";
        treeFillPF.setEcalBarrelRecHits(ecalBarrelRecHits_);
        treeFillPF.setEcalEndcapRecHits(ecalEndcapRecHits_);
        treeFillPF.setEcalSuperClusters(ecalElePFClusterCollection_);
        treeFillPF.writeCollectionToTree(ecalElePFClusterCollection_, iEvent, iSetup, prefix, barrelSuffix, false);
        
        if(dumpPhoPFClusters_) {
          CmsBasicClusterFiller treeFillPhoPF(tree_, 350);
          prefix="";
          barrelSuffix="PhoPFBC";
          treeFillPF.setEcalBarrelRecHits(ecalBarrelRecHits_);
          treeFillPF.setEcalEndcapRecHits(ecalEndcapRecHits_);
          treeFillPF.setEcalSuperClusters(ecalPhoPFClusterCollection_);
          treeFillPF.writeCollectionToTree(ecalPhoPFClusterCollection_, iEvent, iSetup, prefix, barrelSuffix, false);
        }
      }
  }

  // fill general track block
  if(dumpTracks_) {

    CmsTrackFiller treeFiller(tree_, vertexCollection_, true);
    treeFiller.saveFatTrk(saveFatTrk_);

    treeFiller.isGsf(false);
    treeFiller.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    treeFiller.saveDeDx(saveTrackDeDx_);

    treeFiller.saveVtxTrk(saveTrackImpactParameters_);

    std::string prefix("");
    std::string suffix("GeneralTrack");

    treeFiller.writeCollectionToTree(generalTrackCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  // fill the leptonLinkedTracks
  if(dumpLinkedTracks_) {
    CmsTrackFiller leptonLLTreeFiller(tree_, vertexCollection_, true);
    leptonLLTreeFiller.saveFatTrk(saveFatTrk_);
    leptonLLTreeFiller.isGsf(false);
    leptonLLTreeFiller.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    leptonLLTreeFiller.saveDeDx(saveTrackDeDx_);
    leptonLLTreeFiller.saveVtxTrk(saveTrackImpactParameters_);
    std::string prefix("");
    std::string suffix("Track");
    leptonLLTreeFiller.writeCollectionToTree(trackCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  if(dumpGsfTracks_) {

    CmsGsfTrackFiller treeFiller(tree_, vertexCollection_, true);
    treeFiller.saveFatTrk(saveFatTrk_);
    treeFiller.setRefittedTracksForDeDxProducer(gsfTrackCollection_);
    treeFiller.isGsf(true);
    treeFiller.saveDeDx(false);
    treeFiller.saveVtxTrk(saveTrackImpactParameters_);


    std::string prefix("");
    std::string suffix("GsfTrack");

    treeFiller.writeCollectionToTree(gsfTrackCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  if(dumpMuonTracks_) {

    // dump tracks reconstructed in both tracked and muon detector
    CmsTrackFiller treeFillerGlobalMuonTrack(tree_, vertexCollection_, true);
    treeFillerGlobalMuonTrack.saveFatTrk(saveFatTrk_);
    treeFillerGlobalMuonTrack.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    treeFillerGlobalMuonTrack.saveDeDx(false);
    treeFillerGlobalMuonTrack.isGsf(false);
    treeFillerGlobalMuonTrack.saveVtxTrk(saveTrackImpactParameters_);

    std::string prefix("");
    std::string suffix1("GlobalMuonTrack");
    treeFillerGlobalMuonTrack.writeCollectionToTree(globalMuonTrackCollection_, iEvent, iSetup, prefix, suffix1, false);

    // dump tracks reconstructed in the muon detector only
    CmsTrackFiller treeFillerSTAMuonTrack(tree_, vertexCollection_, true);
    treeFillerSTAMuonTrack.saveFatTrk(saveFatTrk_);
    treeFillerSTAMuonTrack.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    treeFillerSTAMuonTrack.saveDeDx(false);
    treeFillerSTAMuonTrack.isGsf(false);
    treeFillerSTAMuonTrack.saveVtxTrk(false);

    std::string suffix2("STAMuonTrack");
    treeFillerSTAMuonTrack.writeCollectionToTree(standAloneMuonTrackCollection_, iEvent, iSetup, prefix, suffix2, false);

  }

  //fill Primary Vertex and associated tracks
  if(dumpVertices_){
    CmsVertexFiller treeFillerVertices(tree_, true);
    treeFillerVertices.setChargedMet(chargedMetCollection_);
    treeFillerVertices.setGeneralTracksCollection(generalTrackCollection_);
    std::string prefix("");
    std::string suffix("PV");
    treeFillerVertices.writeCollectionToTree(vertexCollection_, iEvent, iSetup, prefix, suffix);
  }

  //fill V0 candidates and associated daughter tracks indices
  if(dumpK0s_){
    CmsV0CandidateFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("K0s");
    treeFill.saveCand(saveCand_);
    treeFill.writeCollectionToTree(K0sCollection_, iEvent, iSetup, prefix, suffix);
  }

  // fill muons block
  if(dumpMuons_) {
    CmsMuonFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Muon");
    treeFill.setGeneralTracks(trackCollection_);
    treeFill.setVertexCollection(vertexCollection_);
    treeFill.setMuonIsoMVA(fMuonIsoMVA);
    treeFill.saveCand(saveCand_);
    treeFill.saveFatTrk(saveFatTrk_);
    // for the fully calibrated muons
    treeFill.setCalibMuonCollection(calibMuonCollection_);
    treeFill.writeCollectionToTree(muonCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  // PF candidates
  if(dumpPFCandidates_) {
    CmsPFCandidateFiller treeFill(tree_);
    std::string prefix("");
    std::string suffix("PFCand");
    treeFill.saveCand(saveCand_);
    treeFill.setGeneralTracks(trackCollection_);
    // for isolation
    treeFill.setMinPtPFCand(10.);
    treeFill.dumpChargeOnly(true);
    treeFill.writeCollectionToTree(PFNoPUCandidateCollection_, PFPUCandidateCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  // PF candidates. Only the ones in a cone 0.5 from leptons to correct ChMET. In Candidate format to save space
  if(dumpParticleFlowObjects_) {
    CmsCandidateFiller treeFill(tree_);
    std::string prefix("");
    std::string suffix("ReducedPFCand");
    treeFill.writeCollectionToTree(leptonLinkedPFCandidates_, iEvent, iSetup, prefix, suffix, false);
  }

  // fill CaloTower block
  if(dumpCaloTowers_){
    CmsCaloTowerFiller treeFill(tree_, hbheLabel_, hoLabel_, hfLabel_, ecalLabels_, true);
    std::string prefix("");
    std::string suffix("CaloTowers");
    treeFill.saveCand(dumpCaloTowers_);
    treeFill.saveCaloTowerExtras(dumpCaloTowers_);
    treeFill.writeCollectionToTree(calotowerCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  // fill ECAL rechits collections
  if(dumpEcalRecHits_) {
    CmsEcalRecHitFiller ebFill(tree_);
    std::string prefix("");
    std::string suffix("EBRecHits");
    ebFill.writeCollectionToTree(ebRHLabel_, iEvent, iSetup, prefix, suffix, false);

    CmsEcalRecHitFiller eeFill(tree_);
    suffix = std::string("EERecHits");
    ebFill.writeCollectionToTree(eeRHLabel_, iEvent, iSetup, prefix, suffix, false);
  }

  // fill MET block
  if(dumpMet_) {

    // Calo MET
    CmsMetFiller treeRecoFill1(tree_, true);
    std::string prefix("");
    std::string suffix("Met");
    treeRecoFill1.saveCand(saveCand_);
    treeRecoFill1.isData(!dumpMCTruth_);
    treeRecoFill1.writeCollectionToTree(metCollection_, iEvent, iSetup, prefix, suffix, false);

    // Corrected CALO MET
    // CmsCandidateFiller treeRecoFill1bis(tree_, true);
    // suffix = "CorrMet";
    // treeRecoFill1bis.saveCand(saveCand_);
    // treeRecoFill1bis.writeCollectionToTree(corrmetCollection_, iEvent, iSetup, prefix, suffix, false);

    // Track-Corrected MET
    // CmsMetFiller treeRecoFill2(tree_, true);
    // suffix = "TCMet";
    // treeRecoFill2.isData(false); // the met flags are per event, dumped for caloMET
    // treeRecoFill2.saveCand(saveCand_);
    // treeRecoFill2.writeCollectionToTree(TCmetCollection_, iEvent, iSetup, prefix, suffix, false);

    // particle flow met
    // [0] = uncorrected PF met
    // [0] = Type-0 corrected PFMET
    // [1] = Type-0 and Type-I corrected PFMET
    // [2] = Type-0, Type-I, and Type-II corrected PFMET
    if ( dumpParticleFlowObjects_ ) {
      CmsMetFiller pfMetFiller(tree_, 4, 4, true);
      suffix = "PFMet";
      pfMetFiller.saveCand(saveCand_);
      pfMetFiller.isData(false); // the met flags are per event, dumped for caloMET
      pfMetFiller.writeCollectionToTree(PFmetCollection_, iEvent, iSetup, prefix, suffix, false);

      // charged PF MET HWW version
      CmsMetFiller treeRecoFill1(tree_, true);
      std::string prefix("");
      std::string suffix("PFChMet");
      treeRecoFill1.isData(false); // the met flags are per event, dumped for caloMET
      treeRecoFill1.saveCand(saveCand_);
      treeRecoFill1.writeCollectionToTree(PFChMetCollection_, iEvent, iSetup, prefix, suffix, false);
    }

    // dump generated MET
    if(dumpGenMet_) {

      CmsCandidateFiller treeGenFill(tree_, true);
      std::string suffix("GenMet");
      treeGenFill.writeCollectionToTree(genMetCollection_, iEvent, iSetup, prefix, suffix, false);

    }

  }


  // fill JET block
  if(dumpJets_) {

    // Calo jets: not used for the moment
    CmsJetFiller caloJetFiller(tree_, true);
    std::string prefix("");
    std::string suffix("AK5Jet");
    caloJetFiller.saveCand(saveCand_);
    caloJetFiller.saveJetExtras(true);
    caloJetFiller.saveJetBTag(saveJetBTag_);
    caloJetFiller.setJetCorrectionService(JetCorrectionService_);
    caloJetFiller.writeCollectionToTree(jetCollection1_, iEvent, iSetup, prefix, suffix, false, jetCollection2_);

  }

  if ( dumpParticleFlowObjects_ && dumpCHSPFJet_ ) { 

    // particle flow jets. These take as input the PFnoPU candidates. 
    CmsPFJetFiller pfJetFiller(tree_, true);
    std::string prefix("");
    std::string suffix = "AK5PFNoPUJet";
    pfJetFiller.saveCand(saveCand_);
    pfJetFiller.saveJetBTag(saveJetBTag_);
    pfJetFiller.setBTags(PFJetsBTags_);
    pfJetFiller.setjetMVAAlgos(_jetId_algos);
    pfJetFiller.setVertexCollection(vertexCollection_);
    pfJetFiller.setJetCorrectionService(PFJetCorrectionService_);
    pfJetFiller.writeCollectionToTree(PFjetCollection1_, iEvent, iSetup, prefix, suffix, false);

  }

  // particle flow jets with correction for pileup
  if ( dumpParticleFlowObjects_ && dumpPUcorrPFJet_ ) {

    CmsPFJetFiller pfPUcorrJetFiller(tree_, true);
    std::string prefix("");
    std::string suffix = "AK5PFPUcorrJet";
    pfPUcorrJetFiller.saveCand(saveCand_);
    pfPUcorrJetFiller.saveJetBTag(saveJetBTag_);
    pfPUcorrJetFiller.setBTags(PFPUcorrJetsBTags_);
    pfPUcorrJetFiller.setjetMVAAlgos(_jetId_algos);
    pfPUcorrJetFiller.setVertexCollection(vertexCollection_);
    pfPUcorrJetFiller.setJetCorrectionService(PFJetCorrectionService_);
    pfPUcorrJetFiller.writeCollectionToTree(PFpuCorrJetCollection1_, iEvent, iSetup, prefix, suffix, false);

  }

  // Jet Plus Tracks jets: not used for the moment
  //     CmsJPTJetFiller jptJetFiller(tree_, true);
  //     suffix = "AK5JPTJet";
  //     jptJetFiller.saveCand(saveCand_);
  //     jptJetFiller.saveJetBTag(saveJetBTag_);
  //     jptJetFiller.writeCollectionToTree(JPTjetCollection1_, iEvent, iSetup, prefix, suffix, false, JPTjetCollection2_);
  
  // dump generated JETs
  if(dumpGenJets_) {
    
    CmsJetFiller genJetFiller(tree_, true);
    std::string prefix("");
    std::string suffix = "AK5GenJet";
    genJetFiller.saveJetExtras(false);
    genJetFiller.saveJetBTag(false);
    genJetFiller.setJetCorrectionService(JetCorrectionService_);
    genJetFiller.isGenJets(true);
    genJetFiller.writeCollectionToTree(genJetCollection_, iEvent, iSetup, prefix, suffix, false);
    
  }

  // dump infos on MC production 
  if (dumpGenInfo_) {

    Handle<GenEventInfoProduct> gei;
    iEvent.getByLabel( "generator", gei );

    CmsGenInfoFiller treeFill(tree_);
    treeFill.writeGenInfoToTree( gei );

  }
  
  
  // dump Hcal noise flags
  if(dumpHcalNoiseFlags_) {

    CmsHcalNoiseFiller treeFill(tree_, true);

    treeFill.writeHcalNoiseSummaryToTree(hbheLabel_, hfLabel_, hcalNoiseSummaryLabel_, iEvent, iSetup,
      aodHcalNoiseFlags_);
  }

 
  if(dumpTree_) tree_->dumpData();

}



// ------------ method called once each job just before starting event loop  ------------
void HWWTreeDumper::beginJob() {
  
  fileOut_ = TFile::Open(nameFile_.c_str(), "RECREATE");

  tree_  = new  CmsTree(nameTree_.c_str(),nameTree_.c_str());

  jevt_ = 1;

  //HLTObject Filler needs to exist before beginRun is called
  if(dumpHLTObject_)
    hltObjectFiller_ = new CmsHLTObjectFiller(tree_,hltParms_);
  else hltObjectFiller_ = 0;

  // this pointer MUST survive until tree is closed
  trgNames_ = new vector<std::string>;
  LHEComments_ = new vector<std::string>;

  // initialize MVAs...
  // electron IDs
  std::vector<std::string> myManualCatWeigths;
  myManualCatWeigths.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml").fullPath());
  myManualCatWeigths.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml").fullPath());

  Bool_t manualCat = true;
  
  fEleIdMVANonTrig = new EGammaMvaEleEstimator();
  fEleIdMVANonTrig->initialize("BDT",
                               EGammaMvaEleEstimator::kNonTrig,
                               manualCat, 
                               myManualCatWeigths);
  
  // NOTE: it is better if you copy the MVA weight files locally. See the previous remark
  std::vector<std::string> myManualCatWeigthsTrig;
  myManualCatWeigthsTrig.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV0_Cat1.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV0_Cat2.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV0_Cat3.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV0_Cat4.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV0_Cat5.weights.xml").fullPath());
  myManualCatWeigthsTrig.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV0_Cat6.weights.xml").fullPath());

  fEleIdMVATrig = new EGammaMvaEleEstimator();
  fEleIdMVATrig->initialize("BDT",
                        EGammaMvaEleEstimator::kTrig,
                        manualCat,
                        myManualCatWeigthsTrig);

  std::vector<std::string> myManualCatWeigthsTrigIdIsoCombined;
  myManualCatWeigthsTrigIdIsoCombined.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV1_Cat1.weights.xml").fullPath());
  myManualCatWeigthsTrigIdIsoCombined.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV1_Cat2.weights.xml").fullPath());
  myManualCatWeigthsTrigIdIsoCombined.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV1_Cat3.weights.xml").fullPath());
  myManualCatWeigthsTrigIdIsoCombined.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV1_Cat4.weights.xml").fullPath());
  myManualCatWeigthsTrigIdIsoCombined.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV1_Cat5.weights.xml").fullPath());
  myManualCatWeigthsTrigIdIsoCombined.push_back(edm::FileInPath("HiggsAnalysis/HiggsToWW2e/data/Electrons_BDTG_TrigV1_Cat6.weights.xml").fullPath());

  fEleIdMVATrigIdIsoCombined = new EGammaMvaEleEstimator();
  fEleIdMVATrigIdIsoCombined->initialize("BDT",
                                         EGammaMvaEleEstimator::kTrigIDIsoCombined,
                                         manualCat,
                                         myManualCatWeigthsTrigIdIsoCombined);


  // muon isolation
  fMuonIsoMVA = new MuonMVAEstimator();
  vector<string> muoniso_weightfiles;
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml").fullPath());                 
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml").fullPath());
  muoniso_weightfiles.push_back(edm::FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml").fullPath());
  fMuonIsoMVA->initialize("MuonIso_BDTG_IsoRings",
                          MuonMVAEstimator::kIsoRings,
                          true,
                          muoniso_weightfiles);
  fMuonIsoMVA->SetPrintMVADebug(false);

  // jet ID
  _jetId_algos.resize(puJetIDAlgos_.size());
  for(unsigned int imva=0; imva<puJetIDAlgos_.size(); imva++){
    _jetId_algos[imva] = new PileupJetIdAlgo((puJetIDAlgos_.at(imva)));
  }

}

void HWWTreeDumper::beginRun( const Run & iRun, const EventSetup & iSetup )
{
  jevtInRun_ = 0;
}



// ------------ method called once each job just after ending the event loop  ------------
void  HWWTreeDumper::endJob() {

  fileOut_->cd();

  TTree* treeEventsOut = tree_->getTree();
  treeEventsOut->Write();

  fileOut_->Close();

  delete fEleIdMVATrig;
  delete fEleIdMVATrigIdIsoCombined;
  delete fEleIdMVANonTrig;
  delete fMuonIsoMVA;

}


