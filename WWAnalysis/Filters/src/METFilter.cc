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

#include <DataFormats/METReco/interface/MET.h>
#include <DataFormats/METReco/interface/METFwd.h>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"

// from ROOT
#include "Math/VectorUtil.h"

//
// class declaration
//

class METFilter : public edm::EDFilter {
   public:
      explicit METFilter(const edm::ParameterSet&);
      ~METFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag inputMET_;
  edm::InputTag inputLeptons_;
  double projMetCut_,MetCut_;
  bool skipMetCut_;

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
METFilter::METFilter(const edm::ParameterSet& iConfig):
  inputMET_(iConfig.getParameter<edm::InputTag>("srcMET")),
  inputLeptons_(iConfig.getParameter<edm::InputTag>("srcLeptons")),
  projMetCut_(iConfig.getParameter<double>("projectedMETCut")),
  MetCut_(iConfig.getParameter<double>("METCut")),
  skipMetCut_(iConfig.getParameter<bool>("skipMEtCut"))
{
   //now do what ever initialization is needed

}


METFilter::~METFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
METFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   
   Handle<View<MET> > met;
   iEvent.getByLabel(inputMET_,met);

   Handle< vector< Ptr<Candidate> > > cands;
   iEvent.getByLabel(inputLeptons_,cands);

   Candidate::LorentzVector metP4 = met->front().p4();



   double smallestDphi(999.);
   double projectedMET(0);
   for (vector< Ptr<Candidate> >::const_iterator it = cands->begin(), ed = cands->end();
	it != ed; ++it){
     double dphi = abs(ROOT::Math::VectorUtil::DeltaPhi(metP4,(*it)->p4()) );
     if( dphi < smallestDphi){
       if(dphi < M_PI/2.)
	 projectedMET = metP4.pt()*sin(dphi);
       else
	 projectedMET = metP4.pt();
       
       smallestDphi = dphi;
     }
   }

   if(skipMetCut_){
     return (projectedMET >projMetCut_);
   }else{
     return (projectedMET >projMetCut_ && metP4.pt() > MetCut_);
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
METFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METFilter);
