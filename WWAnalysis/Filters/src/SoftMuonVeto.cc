// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"

// from ROOT
#include "Math/VectorUtil.h"

//
// class declaration
//

class SoftMuonVeto : public edm::EDFilter {
   public:
      explicit SoftMuonVeto(const edm::ParameterSet&);
      ~SoftMuonVeto();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag inputCompCandsName_;
  edm::InputTag inputMuons_;
  int minTrkHits_;
  double minPt_,d0Cut_;
  

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SoftMuonVeto::SoftMuonVeto(const edm::ParameterSet& iConfig):
  inputCompCandsName_(iConfig.getParameter<edm::InputTag>("srcCompCands")),
  inputMuons_(iConfig.getParameter<edm::InputTag>("srcMuons")),
  minTrkHits_(iConfig.getParameter<int>("minTrkHits")),
  minPt_(iConfig.getParameter<double>("minPt")),
  d0Cut_(iConfig.getParameter<double>("d0Cut"))
{
  //now do what ever initialization is needed

}


SoftMuonVeto::~SoftMuonVeto()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SoftMuonVeto::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   
   /*
   cout << "this event arrived to the softMuonVeto cut run/ls/evt: "
	<< iEvent.eventAuxiliary().run() << " / " 
	<< iEvent.eventAuxiliary().luminosityBlock() << " / " 
	<< iEvent.eventAuxiliary().event() << endl;
   */

   Handle<VertexCollection> vertices;
   iEvent.getByLabel("offlinePrimaryVertices",vertices);

   Handle<View<Muon> > muons;
   iEvent.getByLabel(inputMuons_,muons);

   Handle<vector<CompositeCandidate> > compCandidates;
   iEvent.getByLabel(inputCompCandsName_,compCandidates);


    Vertex::Point pvPos = vertices->front().position();

    for( View<Muon>::const_iterator mu = muons->begin(); 
         mu != muons->end(); ++ mu ){
      if ( mu->pt() > minPt_         &&
	   mu->isTrackerMuon()       &&
	   muon::isGoodMuon( *mu, muon::TMLastStationAngTight ) &&
	   mu->innerTrack()->found() > 10 &&
	   fabs( mu->innerTrack()->dxy(pvPos) ) < d0Cut_
	   ) {
	bool doesNotMatch(false);
	for(vector<CompositeCandidate>::const_iterator it=compCandidates->begin(),
	      ed=compCandidates->end(); it!=ed; ++it)
	  {	    
	    CandidateBaseRef ref1(it->daughter(0)->masterClone());
	    CandidateBaseRef ref2(it->daughter(1)->masterClone());

	    /// --- poor-man matching between muon
	    if( (ref1->pt() != mu->pt() ||
		 ref1->eta() != mu->eta() ) &&
		(ref2->pt() != mu->pt() ||
		 ref2->eta() != mu->eta() ) ) doesNotMatch = true ;
	  }
	if(doesNotMatch) return false;
      }          
    }
   return  true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SoftMuonVeto::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SoftMuonVeto::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SoftMuonVeto);
