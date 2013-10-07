// -*- C++ -*-
//
// Package:    SkimEventCombiner
// Class:      SkimEventCombiner
// 
/**\class SkimEventCombiner SkimEventCombiner.cc WWAnalysis/SkimEventCombiner/src/SkimEventCombiner.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu Jun  9 19:06:34 CEST 2011
// $Id$
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

#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

//
// class declaration
//

class SkimEventCombiner : public edm::EDProducer {
   public:
      explicit SkimEventCombiner(const edm::ParameterSet&);
      ~SkimEventCombiner();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::InputTag collATag_;
  edm::InputTag collBTag_;
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
SkimEventCombiner::SkimEventCombiner(const edm::ParameterSet& iConfig):
  collATag_(iConfig.getParameter<edm::InputTag>("inputA")),
  collBTag_(iConfig.getParameter<edm::InputTag>("inputB"))
{
  using namespace reco;
   //register your products
    produces<std::vector<reco::SkimEvent> >().setBranchAlias(iConfig.getParameter<std::string>("@module_label"));  
}


SkimEventCombiner::~SkimEventCombiner()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SkimEventCombiner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;
   std::auto_ptr<std::vector<reco::SkimEvent> > mergedColl(new std::vector<reco::SkimEvent> );


   edm::Handle<vector<SkimEvent > > collA;
   iEvent.getByLabel(collATag_,collA);

   edm::Handle<vector<SkimEvent > > collB;
   iEvent.getByLabel(collBTag_,collB);

   for(vector<SkimEvent >::const_iterator it=collA->begin(); it!=collA->end(); ++it){
     mergedColl->push_back( *it);
   }

   for(vector<SkimEvent >::const_iterator it=collB->begin(); it!=collB->end(); ++it){
     mergedColl->push_back( *it);
   }

    iEvent.put(mergedColl);
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
SkimEventCombiner::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimEventCombiner::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
SkimEventCombiner::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SkimEventCombiner::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SkimEventCombiner::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SkimEventCombiner::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SkimEventCombiner::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimEventCombiner);
