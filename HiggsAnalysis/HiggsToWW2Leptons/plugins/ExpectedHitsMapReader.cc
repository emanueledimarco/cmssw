// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/Common/interface/ValueMap.h>
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"


//
// class declaration
//

class ExpectedHitsMapReader : public edm::EDAnalyzer {
   public:
      explicit ExpectedHitsMapReader(const edm::ParameterSet&);
      ~ExpectedHitsMapReader();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
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
ExpectedHitsMapReader::ExpectedHitsMapReader(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


ExpectedHitsMapReader::~ExpectedHitsMapReader()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ExpectedHitsMapReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

   Handle<GsfElectronCollection> electrons;
   iEvent.getByLabel("gsfElectrons",electrons);

   Handle<MuonCollection> muons;
   iEvent.getByLabel("muons",muons);

   edm::Handle<edm::ValueMap<int> > vmEl;
   iEvent.getByLabel("expectedHitsEle",vmEl);     
   
   edm::Handle<edm::ValueMap<int> > vmMu;
   iEvent.getByLabel("expectedHitsMu",vmMu);     


   for(unsigned int i=0;i<electrons->size();++i) {
     reco::GsfElectronRef gsf(electrons, i);     
     cout << "Electron pt, # of missing hits on expected-crossed-layers: " 
	  << gsf->pt() << " , " <<  (*vmEl)[gsf] << endl;
   }

   for(unsigned int i=0;i<muons->size();++i) {
     reco::MuonRef mu(muons, i);     
     cout << "Muon pt, # of missing hits on expected-crossed-layers: " 
	  << mu->pt() << " , " <<  (*vmMu)[mu] << endl;
   }


}


// ------------ method called once each job just before starting event loop  ------------
void 
ExpectedHitsMapReader::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExpectedHitsMapReader::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExpectedHitsMapReader);
