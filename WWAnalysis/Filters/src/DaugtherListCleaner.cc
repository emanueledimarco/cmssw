// -*- C++ -*-
//
// Package:    DaugtherListCleaner
// Class:      DaugtherListCleaner

//
// Original Author:  Boris Mangano
//         Created:  Wed Jun 30 12:21:29 CEST 2010
// $Id: DaugtherListCleaner.cc,v 1.1 2010/06/30 13:54:18 mangano Exp $
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

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"


//
// class declaration
//

class DaugtherListCleaner : public edm::EDProducer {
   public:
      explicit DaugtherListCleaner(const edm::ParameterSet& pset);
      ~DaugtherListCleaner();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag inputCollName_;
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
DaugtherListCleaner::DaugtherListCleaner(const edm::ParameterSet& iConfig):
  inputCollName_(iConfig.getParameter<edm::InputTag>("src"))
{
  produces< std::vector<edm::Ptr<reco::Candidate> > >(); 
}


DaugtherListCleaner::~DaugtherListCleaner()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DaugtherListCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  // --- input
  Handle<reco::CandidateView> compCandidates;
  iEvent.getByLabel(inputCollName_,compCandidates);
  
  // --- output
  std::auto_ptr< vector< Ptr<Candidate> > > daugthVect(new vector< Ptr<Candidate> >);



  set< Ptr<Candidate> > daugthSet; //use a std::set in order to prune away duplicates

  for(reco::CandidateView::const_iterator it=compCandidates->begin(),
	ed=compCandidates->end(); it!=ed; ++it)
    {
      const CompositeCandidate* comp = dynamic_cast<const CompositeCandidate*>(&*it);
      if(!comp) {
	throw cms::Exception("DaugtherListCleaner") 
	  << "input collection doesn't contain composite candidates" ;
      }else{
	CandidateBaseRef ref1(comp->daughter(0)->masterClone());
	CandidateBaseRef ref2(comp->daughter(1)->masterClone());
	
	Ptr<Candidate> ptr1(ref1.id(), ref1.get(), ref1.key());
	Ptr<Candidate> ptr2(ref2.id(), ref2.get(), ref2.key());
	
	daugthSet.insert(ptr1); 
        daugthSet.insert(ptr2); 
      }
   }
  
  for(set< Ptr<Candidate> >::const_iterator it=daugthSet.begin(),
	ed=daugthSet.end(); it!=ed; ++it){
    daugthVect->push_back(*it);
  }

  // --- eventually put in the event
  iEvent.put(daugthVect);

 
}

// ------------ method called once each job just before starting event loop  ------------
void 
DaugtherListCleaner::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DaugtherListCleaner::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DaugtherListCleaner);
