// -*- C++ -*-
//
// Package:    HeavyFlavourFilter
// Class:      HeavyFlavourFilter
// 
/**\class HeavyFlavourFilter HeavyFlavourFilter.cc WWAnalysis/HeavyFlavourFilter/src/HeavyFlavourFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Adish Vartak,32 2-B09,+41227674752,
//         Created:  Sat Jun  9 23:20:48 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//
// class declaration
//

class HeavyFlavourFilter : public edm::EDFilter {
   public:
      explicit HeavyFlavourFilter(const edm::ParameterSet&);
      ~HeavyFlavourFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      edm::InputTag gensTag;
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
HeavyFlavourFilter::HeavyFlavourFilter(const edm::ParameterSet& iConfig):
    gensTag(iConfig.getParameter<edm::InputTag>("gensTag"))
{
   //now do what ever initialization is needed

}


HeavyFlavourFilter::~HeavyFlavourFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HeavyFlavourFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<GenParticleCollection> gensH;
    iEvent.getByLabel(gensTag, gensH);
    GenParticleCollection gens = *gensH;

    bool hasZDaughter = false;
    bool hasb = false;
    bool hasbbar = false;
    bool hasc = false;
    bool hascbar = false;

    for (size_t i = 0 ; i < gens.size(); i++) {
        if (gens[i].pdgId() == 5 && (gens[i].status() == 2 || gens[i].status() == 3)) {
            hasb = true;
            for( size_t j = 0; j < gens[i].numberOfDaughters(); j++) {
                if(gens[i].daughter(j) && gens[i].daughter(j)->pdgId() == 23) hasZDaughter = true;
            }    
        }
        else if (gens[i].pdgId() == -5 && (gens[i].status() == 2 || gens[i].status() == 3)) {
            hasbbar = true;
            for( size_t j = 0; j < gens[i].numberOfDaughters(); j++) {
                if(gens[i].daughter(j) && gens[i].daughter(j)->pdgId() == 23) hasZDaughter = true;
            }    
        }
        else if (gens[i].pdgId() == 4 && (gens[i].status() == 2 || gens[i].status() == 3)) {
            hasc = true;
            for( size_t j = 0; j < gens[i].numberOfDaughters(); j++) {
                if(gens[i].daughter(j) && gens[i].daughter(j)->pdgId() == 23) hasZDaughter = true;
            }    
        }
        else if (gens[i].pdgId() == -4 && (gens[i].status() == 2 || gens[i].status() == 3)) {
            hascbar = true;
            for( size_t j = 0; j < gens[i].numberOfDaughters(); j++) {
                if(gens[i].daughter(j) && gens[i].daughter(j)->pdgId() == 23) hasZDaughter = true;
            }    
        }
    }

    if (((hasb && hasbbar) || (hasc && hascbar)) && !hasZDaughter) return true;
    else return false;
   
}

// ------------ method called once each job just before starting event loop  ------------
void 
HeavyFlavourFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyFlavourFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
HeavyFlavourFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
HeavyFlavourFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
HeavyFlavourFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
HeavyFlavourFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HeavyFlavourFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(HeavyFlavourFilter);
