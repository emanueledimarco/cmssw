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

#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"

#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"


// from ROOT
#include "Math/VectorUtil.h"

//
// class declaration
//

class JetVetoFilter : public edm::EDFilter {
   public:
      explicit JetVetoFilter(const edm::ParameterSet&);
      ~JetVetoFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag inputJets_;
  edm::InputTag inputMuons_;
  edm::InputTag inputEles_;
  double etaMax_,etMax_,dRcut_;

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
JetVetoFilter::JetVetoFilter(const edm::ParameterSet& iConfig):
  inputJets_(iConfig.getParameter<edm::InputTag>("srcJets")),
  inputMuons_(iConfig.getParameter<edm::InputTag>("srcMuons")),
  inputEles_(iConfig.getParameter<edm::InputTag>("srcElectrons")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  etMax_(iConfig.getParameter<double>("etMax")),
  dRcut_(iConfig.getParameter<double>("deltaR"))
{
   //now do what ever initialization is needed

}


JetVetoFilter::~JetVetoFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
JetVetoFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;

   /*
   cout << "this event arrived to the JetVetoFilter cut run/ls/evt: "
	<< iEvent.eventAuxiliary().run() << " / " 
	<< iEvent.eventAuxiliary().luminosityBlock() << " / " 
	<< iEvent.eventAuxiliary().event() << endl;
   */

   Handle<View<Jet> > jets;
   iEvent.getByLabel(inputJets_,jets);

   Handle<View<Muon> > muons;
   iEvent.getByLabel(inputMuons_,muons);

   Handle<View<GsfElectron> > eles;
   iEvent.getByLabel(inputEles_,eles);


   bool foundJet(false);
   for(View<Jet>::const_iterator it=jets->begin(); it!=jets->end(); ++it){
     if( fabs(it->eta()) >= etaMax_ ) continue;
     //if( it->et() <= etMax_ ) continue;
     if( it->pt() <= etMax_ ) continue;
     // if here, found a jet passing jet-selection


     Candidate::LorentzVector jetP4 = it->p4();

     bool thisJetIsLepton(false);
     for(View<Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
       double dR = fabs(ROOT::Math::VectorUtil::DeltaR(jetP4,mu->p4()) );
       if(dR < dRcut_){
	 thisJetIsLepton = true;
	 break;
       }
     }
     
     for(View<GsfElectron>::const_iterator ele=eles->begin(); ele!=eles->end(); ++ele){
       double dR = fabs(ROOT::Math::VectorUtil::DeltaR(jetP4,ele->p4()) );
       if(dR < dRcut_){
	 thisJetIsLepton = true;
	 break;
       }
     }
     
     if(!thisJetIsLepton) {
       /*
       cout <<"found a veto jet!" << endl;
       cout << "jet et,pt,eta,phi: " 
	    << it->et() << " , " << jetP4.pt() << " , " 
	    << jetP4.eta() << " , " << jetP4.phi() << endl;
       */
       foundJet = true;
       break;
     }
   }

   return !foundJet ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetVetoFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetVetoFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetVetoFilter);
