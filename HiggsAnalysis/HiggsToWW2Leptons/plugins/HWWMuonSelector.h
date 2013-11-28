#ifndef HWWMUONSELECTOR
#define HWWMUONSELECTOR


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class HWWMuonSelector  {
  public:
     explicit HWWMuonSelector(const edm::ParameterSet&);
     ~HWWMuonSelector();
       typedef reco::MuonCollection collection;
       typedef std::vector<reco::MuonRef>::const_iterator const_iterator;

       const_iterator begin() const{return selected_.begin();}
       const_iterator end() const{return selected_.end();}

       void select(edm::Handle<reco::MuonCollection>,const edm::Event&, 
                   const edm::EventSetup&);

  private:

 // ----------member data ---------------------------
 
  std::vector<reco::MuonRef> selected_;
 double muonPtMin_;
 double muonEtaMax_;

};

#endif
   

