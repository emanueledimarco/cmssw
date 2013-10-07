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

#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"


// from ROOT
#include "Math/VectorUtil.h"

//
// class declaration
//

class ExtraLeptonVeto : public edm::EDFilter {
   public:
      explicit ExtraLeptonVeto(const edm::ParameterSet&);
      ~ExtraLeptonVeto();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag inputCompCandsName_;
  edm::InputTag inputMuons_;
  edm::InputTag inputEles_;
  double dR_;
  

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
ExtraLeptonVeto::ExtraLeptonVeto(const edm::ParameterSet& iConfig):
  inputCompCandsName_(iConfig.getParameter<edm::InputTag>("srcCompCands")),
  inputMuons_(iConfig.getParameter<edm::InputTag>("srcMuons")),
  inputEles_(iConfig.getParameter<edm::InputTag>("srcElectrons")),
  dR_(iConfig.getParameter<double>("dR"))
{
  //now do what ever initialization is needed

}


ExtraLeptonVeto::~ExtraLeptonVeto()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ExtraLeptonVeto::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   
   Handle<reco::CandidateView> compCandidates;
   iEvent.getByLabel(inputCompCandsName_,compCandidates);

   Handle<View<Muon> > muons;
   iEvent.getByLabel(inputMuons_,muons);

   Handle<View<GsfElectron> > eles;
   iEvent.getByLabel(inputEles_,eles);
   
   for(reco::CandidateView::const_iterator it=compCandidates->begin(),
	 ed=compCandidates->end(); it!=ed; ++it)
     {
       const CompositeCandidate* comp = dynamic_cast<const CompositeCandidate*>(&*it);
       if(!comp) 
	 throw cms::Exception("DaugtherListCleaner") 
	   << "input collection doesn't contain composite candidates" ;
	    
       CandidateBaseRef ref1(comp->daughter(0)->masterClone());
       CandidateBaseRef ref2(comp->daughter(1)->masterClone());

       bool foundMatch(false);
       /// --- poor-man matching between muon
       for(View<Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
	 if( (ref1->pt() != mu->pt() || ref1->eta() != mu->eta()) &&
	     (ref2->pt() != mu->pt() || ref2->eta() != mu->eta()) ) foundMatch = true;
       }

       for(View<GsfElectron>::const_iterator ele=eles->begin(); ele!=eles->end(); ++ele){
	 double deltaR1 = ROOT::Math::VectorUtil::DeltaR( ref1->p4(),ele->p4());
	 double deltaR2 = ROOT::Math::VectorUtil::DeltaR( ref2->p4(),ele->p4());
	 if( ((ref1->pt() != ele->pt() || ref1->eta() != ele->eta()) && deltaR1 >=dR_) &&
	     ((ref2->pt() != ele->pt() || ref2->eta() != ele->eta()) && deltaR2 >=dR_) )	   
	   foundMatch = true;
       }
      
       //I'm not sure if the logic is correct when there is more than a single compositeCandidate
       if(foundMatch) return false;       
     }
   return  true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ExtraLeptonVeto::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExtraLeptonVeto::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExtraLeptonVeto);
