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

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

// from ROOT
#include "Math/VectorUtil.h"

//
// class declaration
//

class GenFilterDiBosons : public edm::EDFilter {
public:
  explicit GenFilterDiBosons(const edm::ParameterSet&);
  ~GenFilterDiBosons();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  const reco::Candidate* findMother(const reco::Candidate* cand);
  
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
GenFilterDiBosons::GenFilterDiBosons(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
}


GenFilterDiBosons::~GenFilterDiBosons()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GenFilterDiBosons::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   
   Handle<GenParticleCollection> particles;
   iEvent.getByLabel("genParticles",particles);

   float pt1,pt2;
   typedef vector<GenParticleCollection::const_iterator> vLeps;
   vLeps zLeptons;
   for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
     if(it->status()!=1) continue;
     if(fabs(it->pdgId())!=11 && fabs(it->pdgId())!=13) continue;
     const Candidate* mother = findMother(&*it);
     if(mother->pdgId()!=23) continue;
     zLeptons.push_back(it);
   }
   
   bool foundZpeakingEvent(false);
   for(vLeps::const_iterator it=zLeptons.begin(); it!=zLeptons.end(); ++it){
     pt1 = (*it)->pt();
     const Candidate* mother1 = findMother(&**it);
     for(vLeps::const_iterator it2=it+1; it2!=zLeptons.end(); ++it2){
       if(fabs((*it)->pdgId()) != fabs((*it2)->pdgId())) continue;
       pt2 = (*it2)->pt();
       if( max(pt1,pt2)<19) continue;
       if( min(pt1,pt2)<9) continue;
       const Candidate* mother2 = findMother(&**it2);
       if(mother1 != mother2) continue;
       foundZpeakingEvent = true;

       /*
       if(foundZpeakingEvent){
	 cout << "lep1 id,motherM,pt: " 
	      << (*it)->pdgId() << " , "
	      << mother1->mass() << " , "
	      << (*it)->pt() << endl;
	 cout << "lep2 id,motherM,pt: " 
	      << (*it2)->pdgId() << " , "
	      << mother2->mass() << " , "
	      << (*it2)->pt() << endl;
       }
       */
       break;
     }
   }


   return foundZpeakingEvent;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenFilterDiBosons::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
 void 
  GenFilterDiBosons::endJob() {
}


const reco::Candidate* GenFilterDiBosons::findMother(const reco::Candidate* cand){
  using namespace reco;
  //std::cout << "pdgId in findmother: " << cand->pdgId() << std::endl;

  int id = cand->pdgId();
  const Candidate* mother = cand->mother();
  int motherId = cand->mother()->pdgId();
  while(motherId ==id){
    mother = mother->mother();
    motherId = mother->pdgId();
  }
 
  return mother;

}
//define this as a plug-in
DEFINE_FWK_MODULE(GenFilterDiBosons);
