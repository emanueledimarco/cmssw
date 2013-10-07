// -*- C++ -*-
//
// Package:    PatMuIPSelector
// Class:      PatMuIPSelector

//
// Original Author:  Boris Mangano
//         Created:  Wed Jun 30 12:21:29 CEST 2010
// $Id: PatMuIPSelector.cc,v 1.1 2010/06/30 13:54:18 mangano Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>



//
// class declaration
//

class PatMuIPSelector : public edm::EDProducer {
   public:
      explicit PatMuIPSelector(const edm::ParameterSet& pset);
      ~PatMuIPSelector();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag src_;
  edm::InputTag vtxs_;
  double d0Cut_;
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
PatMuIPSelector::PatMuIPSelector(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<edm::InputTag>("src")),
  vtxs_(iConfig.getParameter<edm::InputTag>("vertexs")),
  d0Cut_( iConfig.getParameter<double>( "d0Cut" ) )

{
  produces< pat::MuonCollection>(); 
}


PatMuIPSelector::~PatMuIPSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PatMuIPSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  // --- input
  Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(src_,muons);

  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vtxs_,vertices);

  reco::Vertex::Point pvPos = vertices->front().position();

  // --- output
  std::auto_ptr< pat::MuonCollection>  selected(new pat::MuonCollection);

  for(pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
    if ( mu->type() != 8         &&
	 fabs( mu->innerTrack()->dxy(pvPos) ) < d0Cut_
	 ) selected->push_back( *mu );

  }



  // --- eventually put in the event
  iEvent.put(selected);

 
}

// ------------ method called once each job just before starting event loop  ------------
void 
PatMuIPSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PatMuIPSelector::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatMuIPSelector);
