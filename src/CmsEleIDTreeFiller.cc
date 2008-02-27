#include <memory>
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleTrackerIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleCaloIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"



//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------

CmsEleIDTreeFiller::CmsEleIDTreeFiller(CmsTree *cmsTree, int maxTracks, bool noOutputIfLimitsReached ):
  CmsCandidateFiller(cmsTree,maxTracks,500,noOutputIfLimitsReached),
  privateData_(new CmsEleIDTreeFillerData)
{
  cmstree=cmsTree;
  maxTracks_=maxTracks;
  trkIndexName_ = new std::string("n");
  
  standalone_ = true;

  privateData_->initialise();
}

//--------------
// Destructor --
//--------------
CmsEleIDTreeFiller::~CmsEleIDTreeFiller() {
  // delete here the vector ptr's
  delete privateData_->eleClass;
  delete privateData_->eleHoE;
  delete privateData_->eleNotCorrEoP;
  delete privateData_->eleCorrEoP;
  delete privateData_->eleNotCorrEoPout;
  delete privateData_->eleCorrEoPout;
  delete privateData_->eleDeltaEtaAtVtx;
  delete privateData_->eleDeltaEtaAtCalo;
  delete privateData_->eleDeltaPhiAtVtx;
  delete privateData_->eleDeltaPhiAtCalo;
  delete privateData_->eleFullCorrE;
  delete privateData_->eleCaloCorrE;
  delete privateData_->eleNxtalCorrE;
  delete privateData_->eleRawE;
  delete privateData_->eleTrackerP;
  delete privateData_->eleTrackerIso_sumPt;
  delete privateData_->eleCaloIso_sumPt;
  delete privateData_->eleLik;
  delete privateData_->eleIdCutBasedDecision;
  delete privateData_->eleTip;
}


//-------------
// Methods   --
//-------------
void CmsEleIDTreeFiller::setStandalone(bool what) { standalone_=what; }


void CmsEleIDTreeFiller::writeCollectionToTree(edm::InputTag collectionTag,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {

  // used for the general tree dump
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsCandidateFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  // used only for the electron ID association map
  try { iEvent.getByLabel(collectionTag, explicitElectronCollectionHandle_); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsCandidateFiller") << "Can't get candidate collection: " << collectionTag; }

  privateData_->clearTrkVectors();
  
  if(collection) {

    // Cluster shape variables --- 
    Handle<BasicClusterShapeAssociationCollection> barrelClShpHandle;
    try { iEvent.getByLabel(EcalBarrelClusterShapes_, barrelClShpHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get ECAL barrel Cluster Shape Collection" << EcalBarrelClusterShapes_; }
    const reco::BasicClusterShapeAssociationCollection& barrelClShpMap = *barrelClShpHandle;

    Handle<BasicClusterShapeAssociationCollection> endcapClShpHandle;
    try { iEvent.getByLabel(EcalEndcapClusterShapes_, endcapClShpHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get ECAL endcap Cluster Shape Collection" << EcalEndcapClusterShapes_; }
    const reco::BasicClusterShapeAssociationCollection& endcapClShpMap = *endcapClShpHandle;

    // Read in electron ID association map
    edm::Handle<reco::ElectronIDAssociationCollection> electronIDAssocHandle;
    try { iEvent.getByLabel(electronIDAssocProducer_, electronIDAssocHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get electronId Assoc Producer" << electronIDAssocProducer_; }
    const reco::ElectronIDAssociationCollection& eleIdAssoc = *electronIDAssocHandle;

    // Read the tracker and HCAL isolation association vectors
    try { iEvent.getByLabel(tkIsolationProducer_, tkIsolationHandle_); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get tracker isolation product" << tkIsolationProducer_; }

    try { iEvent.getByLabel(towerIsolationProducer_, towerIsolationHandle_); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get HCAL tower isolation product" << towerIsolationProducer_; }

    int index=0;
    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      const PixelMatchGsfElectron *electron = dynamic_cast< const PixelMatchGsfElectron * > ( &(*cand) );
      if ( electron != 0 )
	writeEleInfo(electron,index,iEvent,iSetup,barrelClShpMap,endcapClShpMap,eleIdAssoc);
      else edm::LogInfo("CmsEleIDTreeFiller") << "Warning! The collection seems to be not made by "
					 << "electrons, electron-specific infos will be set to default.";
      index++;
    }

  }

  // if used standalone, it is necessary to initialize the size of the block event by event
  if(standalone_) {
    std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
    cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  }

  treeEleInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
}



void CmsEleIDTreeFiller::writeEleInfo(const PixelMatchGsfElectron *electron, int index,
				      const edm::Event& iEvent, const edm::EventSetup& iSetup,
				      const reco::BasicClusterShapeAssociationCollection& barrelClShpMap, 
				      const reco::BasicClusterShapeAssociationCollection& endcapClShpMap,
				      const reco::ElectronIDAssociationCollection& eleIdAssoc) { 

  // --------------------------------------
  // collections needed for isolation - FIXME: take once per event, not once per electron: slow
  //
  // get reconstructed tracks
//   edm::Handle<TrackCollection> tracks;
//   try { iEvent.getByLabel("ctfWithMaterialTracks","",tracks); }
//   catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get collection: " << "ctfWithMaterialTracks"; }

//   /// get hcal cells
//   edm::Handle<HBHERecHitCollection> hcalrhits;
//   try { iEvent.getByLabel("hbhereco", hcalrhits); }
//   catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get collection: " << "hbhereco"; }

//   // taking the calo geometry
//   edm::ESHandle<CaloGeometry> pG;
//   iSetup.get<IdealGeometryRecord>().get(pG);
//   caloGeo = pG.product();
  //--------------------------------------

  SuperClusterRef sclusRef = electron->superCluster();
  // ele corr - notcorr energy
  float myEleTrackerP   = electron->trackMomentumAtVtx().r();
  float myEleFullCorrE  = electron->energy();
  float myEleCaloCorrE  = electron->caloEnergy();
  float myEleNxtalCorrE = 0.;
  float myEleRawE       = 0.;
  float mySeedE         = 0.;
  float mySeedCorrE     = 0.;            
  if(&sclusRef) {
    myEleNxtalCorrE = sclusRef->energy();
    myEleRawE       = sclusRef->rawEnergy();
    mySeedE         = sclusRef->seed()->energy();
    if ((int(electron->classification()/10) == 3) || (int(electron->classification()/10) == 13) ){
      double Ebrem = 0.;  
      basicCluster_iterator bc;
      for(bc = sclusRef->clustersBegin(); bc!=sclusRef->clustersEnd(); bc++) {
	Ebrem = Ebrem +(*bc)->energy();
      }
      Ebrem = Ebrem - mySeedE;
      mySeedCorrE = myEleNxtalCorrE - Ebrem;
    }
    else {mySeedCorrE = myEleFullCorrE;}
  }

  // transverse impact parameter
  GsfTrackRef trRef = electron->gsfTrack();
  float myTip = sqrt((trRef->vertex().x())*(trRef->vertex().x()) + (trRef->vertex().y())*(trRef->vertex().y()));

  // eleID
  privateData_->eleFullCorrE     ->push_back(myEleFullCorrE);  
  privateData_->eleCaloCorrE     ->push_back(myEleCaloCorrE);  
  privateData_->eleNxtalCorrE    ->push_back(myEleNxtalCorrE); 
  privateData_->eleRawE          ->push_back(myEleRawE); 
  privateData_->eleTrackerP      ->push_back(myEleTrackerP);
  privateData_->eleClass         ->push_back(electron->classification());
  privateData_->eleHoE           ->push_back(electron->hadronicOverEm());
  privateData_->eleCorrEoP       ->push_back(electron->eSuperClusterOverP());
  privateData_->eleNotCorrEoP    ->push_back(myEleRawE/myEleTrackerP);
  privateData_->eleCorrEoPout    ->push_back(electron->eSeedClusterOverPout());
  privateData_->eleNotCorrEoPout ->push_back(electron->eSeedClusterOverPout()*(mySeedCorrE/mySeedE));
  privateData_->eleDeltaEtaAtVtx ->push_back(electron->deltaEtaSuperClusterTrackAtVtx());
  privateData_->eleDeltaPhiAtVtx ->push_back(electron->deltaPhiSuperClusterTrackAtVtx());
  privateData_->eleDeltaEtaAtCalo->push_back(electron->deltaEtaSeedClusterTrackAtCalo());
  privateData_->eleDeltaPhiAtCalo->push_back(electron->deltaPhiSeedClusterTrackAtCalo());
  privateData_->eleTip           ->push_back(myTip);

  // electron ID (cut-based and likelihood)
  // for the cut-based, store the decision
  // for the likelihood, store the output of the algorithm
  
  bool hasBarrel=true;
  bool hasEndcap=true;
  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
  seedShpItr = barrelClShpMap.find(sclusRef->seed());

  if(seedShpItr==barrelClShpMap.end()) {
    hasBarrel=false;
    seedShpItr=endcapClShpMap.find(sclusRef->seed());
    if(seedShpItr==endcapClShpMap.end()) hasEndcap=false;
  }

  if(hasBarrel || hasEndcap) {

    edm::Ref< reco::PixelMatchGsfElectronCollection > electronRef(explicitElectronCollectionHandle_, index);

    reco::ElectronIDAssociationCollection::const_iterator electronIDAssocItr;

    if( !(electronRef.isNull()) ) {

      // electron ID
      electronIDAssocItr = eleIdAssoc.find( electronRef );
      
      if ( electronIDAssocItr==eleIdAssoc.end() ) edm::LogWarning("CmsEleIDTreeFiller") << "cannot find the electron id associated with electron";
      
      const reco::ElectronIDRef& id = electronIDAssocItr->val;
      
      privateData_->eleIdCutBasedDecision->push_back( id->cutBasedDecision() );
      privateData_->eleLik->push_back( id->likelihood() );


    }
    
    else {

      edm::LogWarning("CmsEleIDTreeFiller") << "cannot find the electron ref in the electron collection";
      privateData_->eleIdCutBasedDecision->push_back( false );
      privateData_->eleLik->push_back( -1.);
      
    }

  }
  else {
    edm::LogWarning("CmsEleIDTreeFiller") << "cannot find cluster shapes in ECAL barrel or endcap "
					  << " setting electron ID results to the default";
    privateData_->eleIdCutBasedDecision->push_back( false );
    privateData_->eleLik->push_back( -1.);
  }


  // isolations
  double sumPt = (*tkIsolationHandle_)[index].second;
  double sumEt = (*towerIsolationHandle_)[index].second;
  privateData_->eleTrackerIso_sumPt->push_back( sumPt );
  privateData_->eleCaloIso_sumPt->push_back( sumEt );


  // tracker isolation
//   const TrackCollection tracksC = *(tracks.product());
//   hwwEleTrackerIsolation trackIsolation(electron, tracksC);
//   trackIsolation.setExtRadius(0.2);    
//   trackIsolation.setIntRadius(0.015);    
//   float minDR_tracker     = trackIsolation.minDeltaR(0.15);  
//   float minDRveto_tracker = trackIsolation.minDeltaR_withVeto(0.15);  
//   float sumPt_tracker     = trackIsolation.getPtTracks();  

//   std::cout << "sumPt = " << sumPt << "\tsumPt_tracker = " << sumPt_tracker << std::endl;

//   privateData_->eleTrackerIso_minDR->push_back(minDR_tracker);
//   privateData_->eleTrackerIso_minDR_veto->push_back(minDRveto_tracker);
//   privateData_->eleTrackerIso_sumPt->push_back(sumPt_tracker);

  // calo isolation
//   const HBHERecHitCollection hcalRecHits = *(hcalrhits.product());
//   hwwEleCaloIsolation caloIsolation(electron, hcalRecHits, caloGeo);
//   //float minDR_calo = caloIsolation.minDeltaR(0.15);  
//   float minDR_calo = -1000.;
//   caloIsolation.setExtRadius(0.2);    
//   float sumEt_calo = caloIsolation.getEtHcal();  

//   std::cout << "sumEt = " << sumEt << "\tsumEt_calo = " << sumEt_calo << std::endl;

//   privateData_->eleCaloIso_minDR->push_back(minDR_calo);
//   privateData_->eleCaloIso_sumPt->push_back(sumEt_calo);
  
  

}




void CmsEleIDTreeFiller::treeEleInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"eleFullCorrE"+colSuffix).c_str(), *privateData_->eleFullCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCaloCorrE"+colSuffix).c_str(), *privateData_->eleCaloCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNxtalCorrE"+colSuffix).c_str(), *privateData_->eleNxtalCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleRawE"+colSuffix).c_str(), *privateData_->eleRawE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackerP"+colSuffix).c_str(), *privateData_->eleTrackerP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleClass"+colSuffix).c_str(), *privateData_->eleClass, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleHoE"+colSuffix).c_str(), *privateData_->eleHoE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCorrEoP"+colSuffix).c_str(), *privateData_->eleCorrEoP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNotCorrEoP"+colSuffix).c_str(), *privateData_->eleNotCorrEoP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCorrEoPout"+colSuffix).c_str(), *privateData_->eleCorrEoPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNotCorrEoPout"+colSuffix).c_str(), *privateData_->eleNotCorrEoPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaEtaAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaPhiAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaEtaAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaPhiAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackerIso_sumPt"+colSuffix).c_str(), *privateData_->eleTrackerIso_sumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCaloIso_sumPt"+colSuffix).c_str(), *privateData_->eleCaloIso_sumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIdCutBased"+colSuffix).c_str(), *privateData_->eleIdCutBasedDecision, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleLikelihood"+colSuffix).c_str(), *privateData_->eleLik, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTip"+colSuffix).c_str(), *privateData_->eleTip, nCandString.c_str(), 0, "Reco");
}




void CmsEleIDTreeFillerData::initialise() {
  initialiseCandidate();
  eleClass                 = new vector<int>;
  eleHoE                   = new vector<float>;
  eleNotCorrEoP            = new vector<float>;
  eleCorrEoP               = new vector<float>;
  eleNotCorrEoPout         = new vector<float>;
  eleCorrEoPout            = new vector<float>;
  eleDeltaEtaAtVtx         = new vector<float>;
  eleDeltaEtaAtCalo        = new vector<float>;
  eleDeltaPhiAtVtx         = new vector<float>;
  eleDeltaPhiAtCalo        = new vector<float>;
  eleFullCorrE             = new vector<float>;
  eleCaloCorrE             = new vector<float>;
  eleNxtalCorrE            = new vector<float>;
  eleRawE                  = new vector<float>;
  eleTrackerP              = new vector<float>;
  eleTrackerIso_sumPt      = new vector<float>;
  eleCaloIso_sumPt         = new vector<float>;
  eleIdCutBasedDecision    = new vector<bool>;
  eleLik                   = new vector<float>;
  eleTip                   = new vector<float>;
}




void CmsEleIDTreeFillerData::clearTrkVectors() {
  clearTrkVectorsCandidate();
  eleClass                 ->clear();
  eleHoE                   ->clear();
  eleNotCorrEoP            ->clear();
  eleCorrEoP               ->clear();
  eleNotCorrEoPout         ->clear();
  eleCorrEoPout            ->clear();
  eleDeltaEtaAtVtx         ->clear();
  eleDeltaEtaAtCalo        ->clear();
  eleDeltaPhiAtVtx         ->clear();
  eleDeltaPhiAtCalo        ->clear();
  eleFullCorrE             ->clear();
  eleCaloCorrE             ->clear();
  eleNxtalCorrE            ->clear();
  eleRawE                  ->clear();
  eleTrackerP              ->clear();
  eleTrackerIso_sumPt      ->clear();
  eleCaloIso_sumPt         ->clear();
  eleIdCutBasedDecision    ->clear();
  eleLik                   ->clear();
  eleTip                   ->clear();
}

