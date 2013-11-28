#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWMetCorrector.h"

// Lorentz vector
typedef math::XYZTLorentzVector LorentzVector;

HWWMetCorrector::HWWMetCorrector(const edm::ParameterSet& iConfig)
{
  caloMetLabel_    = iConfig.getParameter<edm::InputTag>("CaloMetLabel");
  muonLabel_       = iConfig.getParameter<edm::InputTag>("MuonLabel");

  corrmetMin_      = iConfig.getParameter<double>("CorrMetMin");   
  muonPtMin_       = iConfig.getParameter<double>("MuonPtMin"); 
  muonEtaMax_      = iConfig.getParameter<double>("MuonEtaMax"); 
  muonTrackD0Max_  = iConfig.getParameter<double>("MuonTrackD0Max"); 
  muonTrackDzMax_  = iConfig.getParameter<double>("MuonTrackDzMax"); 
  muonNHitsMin_    = iConfig.getParameter<double>("MuonNHitsMin"); 
  muonDPtMax_      = iConfig.getParameter<double>("MuonDPtMax");
  muonChiSqMax_    = iConfig.getParameter<double>("MuonChiSqMax");

  muonCorrection_  = iConfig.getParameter<bool>("MuonCorrection");
  muonDepositCor_  = iConfig.getParameter<bool>("MuonDepositCor");

  //register your products
  produces<reco::CaloMETCollection>("");
}


HWWMetCorrector::~HWWMetCorrector()
{
}

// ------------ method called to produce the data  ------------
void
HWWMetCorrector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
 
  // get met from CMSSW
  Handle<reco::CaloMETCollection> calomet;
  iEvent.getByLabel(caloMetLabel_,calomet);

  // For muon corrections get muons 
  Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(muonLabel_,muons); 

  // define ethz collections to be stored

  auto_ptr<reco::CaloMETCollection> pOutMets (new reco::CaloMETCollection);

  double DeltaPx = 0.0;
  double DeltaPy = 0.0;
  double DeltaSumET = 0.0;

  double DeltaExDep = 0.0;
  double DeltaEyDep = 0.0;
  double DeltaSumETDep = 0.0;

  for(reco::CaloMETCollection::const_iterator met=calomet->begin(); met != calomet->end(); ++met) {

    // first we check if we have to correct for muons
    if(muonCorrection_) {
      for(reco::MuonCollection::const_iterator muon=muons->begin(); muon != muons->end(); ++muon) {
	bool goodmuon = false;

	TrackRef track;
	bool hasTrackerTrack = false;
	
	if ( muon->track().isNonnull() ) {
	  hasTrackerTrack = true;
	  track = muon->get<TrackRef>();
	}
	
	//if ( hasTrackerTrack && &track != 0 ) {
	if ( hasTrackerTrack && track.isNonnull() ) {
	  goodmuon = track->pt() > muonPtMin_ &&
	    fabs(track->eta()) < muonEtaMax_ &&
	    fabs(track->d0())  < muonTrackD0Max_ &&
	    fabs(track->dz())  < muonTrackDzMax_ && 
	    track->numberOfValidHits() > muonNHitsMin_;
	} else {
	  goodmuon = muon->pt() > muonPtMin_ &&
	    fabs(muon->eta()) < muonEtaMax_;
	}
	if( goodmuon ) {
	  bool goodtrackfit = false;
	  //if ( hasTrackerTrack && &track != 0 ) {
	  if ( hasTrackerTrack && track.isNonnull() ) {
	    float dpt_track = track->error(0)/(track->qoverp());
	    float chisq = track->normalizedChi2();
	    if (dpt_track < muonDPtMax_  && chisq < muonChiSqMax_) goodtrackfit = true;
	  } else goodtrackfit = true;
	  
	  if ( goodtrackfit ) {
	    DeltaPx +=  muon->px();
	    DeltaPy +=  muon->py();
	    DeltaSumET += muon->et();
	  }

	  //----------- Calculate muon energy deposition in the calorimeters

	  /*   if (muonDepositCor_) {
	       TrackRef mu_track = muon->combinedMuon();
	       TrackDetMatchInfo info =  trackAssociator.associate(iEvent, iSetup,trackAssociator.getFreeTrajectoryState(iSetup, *mu_track),
	       trackAssociatorParameters);
	       double ene = info.crossedEnergy(TrackDetMatchInfo::TowerTotal);
	       DeltaExDep += ene*sin((*mu_track).theta())*cos((*mu_track).phi());
	       DeltaEyDep += ene*sin((*mu_track).theta())*sin((*mu_track).phi());
	       DeltaSumETDep += ene*sin((*mu_track).theta());
	       }*/
        }
      }
    }
        						      
    //----------------- Calculate and set deltas for new MET correction
    CorrMETData delta;
    delta.mex     = - DeltaPx;  //correction to MET (from muons) is negative,    
    delta.mey     = - DeltaPy;    //since MET points in direction of muons
    delta.sumet   = DeltaSumET; 
    delta.mex    += DeltaExDep;    //correction to MET (from muon depositions) is positive,	
    delta.mey    += DeltaEyDep;    //since MET points in opposite direction of muons
    delta.sumet  -= DeltaSumETDep;    
     
    //----------------- Fill holder with corrected MET (= uncorrected + delta) values

    double CorrPx = met->px()+delta.mex;
    double CorrPy = met->py()+delta.mey;
    double CorrSumEt = met->sumEt()+delta.sumet;

    double CorrPt=sqrt(CorrPx*CorrPx+CorrPy*CorrPy);
  
    //     if(CorrPt < corrmetMin_) continue;
     
    LorentzVector p4(CorrPx,CorrPy,0,sqrt(CorrPt*CorrPt));

    //----------------- get previous corrections and push into new corrections 

    std::vector<CorrMETData> corrections = met->mEtCorr();
    corrections.push_back( delta );

    //----------------- Push onto MET Collection

    reco::CaloMET* MET = new reco::CaloMET(met->getSpecific(), CorrSumEt,corrections, p4, met->vertex());

    pOutMets->push_back(*MET);

    iEvent.put(pOutMets,"");

    delete MET;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
HWWMetCorrector::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HWWMetCorrector::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HWWMetCorrector);
