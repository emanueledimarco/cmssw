// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

// from ROOT
#include "Math/VectorUtil.h"

//
// class declaration
//

class PileUpSummaryFilter : public edm::EDFilter {
   public:
      explicit PileUpSummaryFilter(const edm::ParameterSet&);
      ~PileUpSummaryFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  // ----------member data ---------------------------
  //edm::InputTag inputMET_;
  int minPU,maxPU;

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
PileUpSummaryFilter::PileUpSummaryFilter(const edm::ParameterSet& iConfig):
  minPU(iConfig.getParameter<int>("minPU")),
  maxPU(iConfig.getParameter<int>("maxPU"))
{
   //now do what ever initialization is needed
}


PileUpSummaryFilter::~PileUpSummaryFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
PileUpSummaryFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   
   Handle<std::vector< PileupSummaryInfo > >  PuInfo;
   iEvent.getByLabel("addPileupInfo", PuInfo);

   int numInTimePU = PuInfo->begin()->getPU_NumInteractions();

   return (numInTimePU>=minPU && numInTimePU<=maxPU);
}

// ------------ method called once each job just before starting event loop  ------------
void 
PileUpSummaryFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PileUpSummaryFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PileUpSummaryFilter);
