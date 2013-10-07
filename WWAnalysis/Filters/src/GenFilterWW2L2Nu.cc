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

class GenFilterWW2L2Nu : public edm::EDFilter {
   public:
      explicit GenFilterWW2L2Nu(const edm::ParameterSet&);
      ~GenFilterWW2L2Nu();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  // ----------member data ---------------------------
  //edm::InputTag inputMET_;


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
GenFilterWW2L2Nu::GenFilterWW2L2Nu(const edm::ParameterSet& iConfig)
  //: inputMET_(iConfig.getParameter<edm::InputTag>("srcMET")),
{
   //now do what ever initialization is needed

}


GenFilterWW2L2Nu::~GenFilterWW2L2Nu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GenFilterWW2L2Nu::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   
   Handle<GenParticleCollection> particles;
   //iEvent.getByLabel("prunedGen",particles);
   iEvent.getByLabel("genParticles",particles);


   bool foundWW(false), foundWZ(false), foundZZ(false);
   bool foundGoodLep1(false);
   bool foundGoodLep2(false);
   int countW(0);
   int countZ(0);

   for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
     if(it->status()!=3) continue;
     if(abs(it->pdgId())==24) countW++;
     if(abs(it->pdgId())==23) countZ++;

   }
   if(countW==2) {foundWW=true; 
     //cout << "found WW" << endl;
   }

   if(countW==1 && countZ==1) {foundWZ=true; 
     //cout << "found WZ" << endl;
   }

   if(countZ==2) {foundZZ=true; 
     //cout << "found ZZ" << endl;
   }


   //if(countW > 2) cout << "LOGIC ERROR: there are more than 2 Ws in the event." << endl; 
   //if(countZ > 2) cout << "LOGIC ERROR: there are more than 2 Zs in the event." << endl; 

   
   if(!(foundWW || foundWZ || foundZZ)){
     Candidate::LorentzVector vect4,vect4B;
     //cout << "WARNING: no simple.Trying WW" << endl;
     
     bool foundL1(false),foundNu1(false);
     bool foundL2(false),foundNu2(false);

     const Candidate* mother1; const Candidate* mother2;
     int pidL1; int pidL2;

     /*
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       //if(it->mother()) 
	 //cout << "pid,mother, mPid: " 
	 // << it->pdgId() << " , " 
	 //   << it->mother() << " , "
	 //   << it->mother()->pdgId() << endl;
     }
     */

     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid==11 || pid==13 || pid == 15){
	 //cout << "got L1. id,mId:" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundL1 = true;
	 mother1 = it->mother();
	 pidL1=pid;
	 break;
       }
     }//end loop particles
     if(!foundL1) goto secondSpecial;
     
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid==(pidL1+1) && it->mother() == mother1) {
	 //cout << "got Nu1" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundNu1 = true;
	 break;
       }
     }    
     if(!foundNu1) goto secondSpecial;

     // now search second pair
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if( (pid==11 || pid==13 || pid == 15) && it->mother() != mother1){
	 //cout << "got L2" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundL2 = true;
	 mother2 = it->mother();
	 pidL2=pid;
	 break;
       }
     }//end loop particles
     if(!foundL2) goto secondSpecial;
     
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid==(pidL2+1) && it->mother() == mother2) {
	 //cout << "got Nu2" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundNu2 = true;
	 break;
       }
     }    
     if(foundNu2) {
       foundWW = true;
       //cout << "found a special WW" << endl;
     }else{
       //cout << "WARNING: still a mystery" << endl;
     }


   }// end special case block


 secondSpecial:
   // --- special case for ZZ --
   if(!(foundWW || foundWZ || foundZZ)){
     Candidate::LorentzVector vect4,vect4B;
     //cout << "WARNING: no simple.Trying ZZ" << endl;
     
     bool foundL1(false),foundNu1(false);
     bool foundL2(false),foundNu2(false);

     const Candidate* mother1; const Candidate* mother2;
     int pidL1; int pidL2;

     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid>=11 && pid <=16){
	 //cout << "got L1A. id,mId:" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundL1 = true;
	 mother1 = it->mother();
	 pidL1=it->pdgId();
	 break;
       }
     }//end loop particles
     if(!foundL1) goto thirdSpecial;
     
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = it->pdgId();
       if(pid==(-pidL1) && it->mother() == mother1) {
	 //cout << "got L1B" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundNu1 = true;
	 break;
       }
     }    
     if(!foundNu1) goto thirdSpecial;


     // now search second pair
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid>=11 && pid <=16 && it->mother() != mother1){
	 //cout << "got L2A" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundL2 = true;
	 mother2 = it->mother();
	 pidL2=it->pdgId();
	 break;
       }
     }//end loop particles
     if(!foundL2) goto thirdSpecial;
     
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = it->pdgId();
       if(pid==(-pidL2) && it->mother() == mother2) {
	 //cout << "got L2B" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 foundNu2 = true;
	 break;
       }
     }    
     if(foundNu2) {
       foundZZ = true;
       //cout << "found a special ZZ" << endl;
     }else{
       //cout << "WARNING: still a mystery" << endl;

     }
   }// end special case block
   
 thirdSpecial:
   // --- special case for WZ --
   if(!(foundWW || foundWZ || foundZZ)){
     Candidate::LorentzVector vect4,vect4B;
     //cout << "WARNING: no simple.Trying WZ" << endl;
     
     bool foundW(false),foundZ(false),foundSecond(false);
     const Candidate* mother1; const Candidate* mother2;
     int pidL1; int pidL2;

     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid>=11 && pid <=16){
	 //cout << "got L1A. id,mId:" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 mother1 = it->mother();
	 pidL1=it->pdgId();
	 break;
       }
     }//end loop particles
     //if(!foundL1) goto thirdSpecial;
     
     
     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       if(it->mother()==mother1 &&
	  abs(it->pdgId())==(pidL1+1)) {
	 foundW=true;
	 break;
       }

       if(it->mother()==mother1 && it->pdgId() == (-pidL1) ) {
	 foundZ=true;       
	 break;
       }
     }
     
     //cout << "foundW: " << foundW << endl;
     //cout << "foundZ: " << foundZ << endl;

     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid>=11 && pid <=16 && it->mother()!=mother1){
	 //cout << "got L2. id,mId:" 
	 //   << pid << " , " << abs(it->mother()->pdgId()) << endl;
	 mother2 = it->mother();
	 pidL2=it->pdgId();
	 foundSecond = true;
	 break;
       }
     }//end loop particles
     //if(!foundL1) goto thirdSpecial;
     
     //cout << "foundSecond: " << foundSecond << endl;

     if(foundSecond){
       for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
	 if(it->status()!=3) continue;
	 if(it->mother()==mother2 &&
	    abs(it->pdgId())==(abs(pidL2)+1)) {foundW=true;break;}
	 if(it->mother()==mother2 &&
	    it->pdgId()==(-pidL2)) {foundZ=true; break;}       
       }
     }
     
     //cout << "foundW: " << foundW << endl;
     //cout << "foundZ: " << foundZ << endl;


     if(foundW && foundZ) {
       foundWZ = true;
       //cout << "found a special WZ" << endl;
     }else{
       //cout << "WARNING: still a mystery" << endl;
     }
   }// end special case block


   // --- very special case  --
   if(!(foundWW || foundWZ || foundZZ)){
     Candidate::LorentzVector vect4,vect4B;
     //cout << "WARNING: no simple.Trying all from one mother" << endl;

     bool first(false);
     int countTo4(0);
     const Candidate* mother1;
     int sumId(0);

     for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
       if(it->status()!=3) continue;
       int pid = abs(it->pdgId());
       if(pid>=11 && pid <=16){
	 if(!first){
	   mother1 = it->mother();
	   sumId += it->pdgId();
	   countTo4++;
	 }else{
	   if(it->mother()==mother1){
	     sumId += it->pdgId();
	     countTo4++;
	   }
	 }
       }
     }
     //cout << "countTo4,sumId: " << countTo4 << " , " << sumId << endl;
     if(countTo4==4 && sumId==0 && abs(mother1->pdgId())==23){
       foundZZ=true;
	 //cout << "found very special ZZ" << endl;
     }else if(countTo4==4 && sumId==0 && abs(mother1->pdgId())==24){
       foundWW=true;
       //cout << "found very special WW" << endl;
     }else if(countTo4==4 && abs(sumId)==1) {
       foundWZ=true;
       //cout << "found very special WZ" << endl;
     }else{
       cout << "still mystery" << endl;
     }
 
   }//end special block
   

   return foundWW;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenFilterWW2L2Nu::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenFilterWW2L2Nu::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenFilterWW2L2Nu);
