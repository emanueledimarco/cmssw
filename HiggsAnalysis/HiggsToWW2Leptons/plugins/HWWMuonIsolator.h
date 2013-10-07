#ifndef HWWMUONISOLATOR
#define HWWMUONISOLATOR

#include "DataFormats/MuonReco/interface/MuonFwd.h"

class HWWMuonIsolator{
public:
 explicit HWWMuonIsolator(const edm::ParameterSet&);
 ~HWWMuonIsolator();

 // Collections to be selected
 typedef reco::MuonCollection collection;
 typedef std::vector<reco::MuonRef>::const_iterator const_iterator;

 //define iterators with above typedef
 const_iterator begin() const{return selected_.begin();}
 const_iterator end() const{return selected_.end();}


 void select(edm::Handle<reco::MuonCollection>,const edm::Event&, 
             const edm::EventSetup&);
private:
   std::vector<reco::MuonRef> selected_;
   double theTrackIsolCut_; 
   bool absolute_;
};


#endif
