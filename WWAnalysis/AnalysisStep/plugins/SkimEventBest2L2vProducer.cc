// -*- C++ -*-
//
// Package:    SkimEventBest2L2vProducer
// Class:      SkimEventBest2L2vProducer
// 
/**\class SkimEventBest2L2vProducer SkimEventBest2L2vProducer.cc WWAnalysis/SkimEventBest2L2vProducer/src/SkimEventBest2L2vProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Fri Nov 23 18:14:42 CET 2012
// $Id: SkimEventBest2L2vProducer.cc,v 1.4 2012/12/19 10:15:05 amassiro Exp $
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

class SkimEventBest2L2vProducer : public edm::EDProducer {
   public:
      explicit SkimEventBest2L2vProducer(const edm::ParameterSet&);
      ~SkimEventBest2L2vProducer();

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
      edm::InputTag mumuTag_;
      edm::InputTag elelTag_;
      edm::InputTag elmuTag_;
      edm::InputTag muelTag_;
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
SkimEventBest2L2vProducer::SkimEventBest2L2vProducer(const edm::ParameterSet& cfg ) :
//     triggerTag_(cfg.getParameter<edm::InputTag>("triggerTag")),
    mumuTag_(cfg.getParameter<edm::InputTag>("mumuEvents")),
    elelTag_(cfg.getParameter<edm::InputTag>("elelEvents")),
    elmuTag_(cfg.getParameter<edm::InputTag>("elmuEvents")),
    muelTag_(cfg.getParameter<edm::InputTag>("muelEvents"))
//     mumuTag_( cfg.getParameter<edm::InputTag> >("mumuEvents") ), 
//     elelTag_( cfg.getParameter<edm::InputTag> >("elelEvents") ), 
//     elmuTag_( cfg.getParameter<edm::InputTag> >("elmuEvents") ), 
//     muelTag_( cfg.getParameter<edm::InputTag> >("muelEvents") )
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
    produces<std::vector<reco::SkimEvent> >().setBranchAlias(cfg.getParameter<std::string>("@module_label"));
}


SkimEventBest2L2vProducer::~SkimEventBest2L2vProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SkimEventBest2L2vProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<reco::SkimEvent> > skimEvent(new std::vector<reco::SkimEvent> );

    edm::Handle<edm::View<reco::SkimEvent> > elel;
    iEvent.getByLabel(elelTag_,elel);
    edm::Handle<edm::View<reco::SkimEvent> > elmu;
    iEvent.getByLabel(elmuTag_,elmu);
    edm::Handle<edm::View<reco::SkimEvent> > muel;
    iEvent.getByLabel(muelTag_,muel);
    edm::Handle<edm::View<reco::SkimEvent> > mumu;
    iEvent.getByLabel(mumuTag_,mumu);


//     std::cout << " mumu = " << mumu->size() << ";  muel = " << muel->size() << "; elmu = " << elmu->size() << "; elel = " << elel->size() << std::endl;

    //---- importance order : mm, mu, em, ee ----
    const reco::SkimEvent* best = 0x0;
    for (size_t i=0; i<mumu->size(); ++i) {//---- mm
     const reco::SkimEvent* candidate =  mumu->ptrAt(i).get();
     if ( !best ) { 
//       std::cout << " mumu " << std::endl;
      best = candidate;
      continue;
     } else {
      double bpt = (best->pt(0)+best->pt(1));
      double cpt = (candidate->pt(0)+candidate->pt(1));
      best = ( bpt > cpt  ? best : candidate);
     }
    }
    for (size_t i=0; i<muel->size(); ++i) {//---- me
     const reco::SkimEvent* candidate =  muel->ptrAt(i).get();
     if ( !best ) { 
//       std::cout << " muel " << std::endl;
      best = candidate;
      continue;
     } else {
      double bpt = (best->pt(0)+best->pt(1));
      double cpt = (candidate->pt(0)+candidate->pt(1));
      best = ( bpt > cpt ? best : candidate);
     }
    }
    for (size_t i=0; i<elmu->size(); ++i) {//---- em
     const reco::SkimEvent* candidate =  elmu->ptrAt(i).get();
     if ( !best ) { 
//       std::cout << " elmu " << std::endl;
      best = candidate;
      continue;
     } else {
      double bpt = (best->pt(0)+best->pt(1));
      double cpt = (candidate->pt(0)+candidate->pt(1));
      best = ( bpt > cpt ? best : candidate);
     }
    }
    for (size_t i=0; i<elel->size(); ++i) {//---- ee
     const reco::SkimEvent* candidate =  elel->ptrAt(i).get();
     if ( !best ) {
//       std::cout << " elel " << std::endl;
      best = candidate;
      continue;
     } else {
      double bpt = (best->pt(0)+best->pt(1));
      double cpt = (candidate->pt(0)+candidate->pt(1));
      best = ( bpt > cpt ? best : candidate);
     }
    }

//     std::cout << "->" << best << " " << mumu->size() << " " << elel->size() << " " << elmu->size() << " " << muel->size() << std::endl;
    if ( best ) skimEvent->push_back( *best );

    iEvent.put(skimEvent); 
}

// ------------ method called once each job just before starting event loop  ------------
void 
SkimEventBest2L2vProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimEventBest2L2vProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
SkimEventBest2L2vProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SkimEventBest2L2vProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SkimEventBest2L2vProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SkimEventBest2L2vProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SkimEventBest2L2vProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimEventBest2L2vProducer);
