#ifndef HWWELECTRONISOLATOR
#define HWWELECTRONISOLATOR

#include <memory>
#include <vector>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


class HWWElectronIsolator{
 public:
   explicit HWWElectronIsolator(const edm::ParameterSet&);
   ~HWWElectronIsolator();

   // Collections to be selected
   typedef reco::GsfElectronCollection collection;
   typedef std::vector<reco::GsfElectronRef> ::const_iterator const_iterator;

   //define iterators with above typedef
   const_iterator begin () const { return selected_.begin () ; }
   const_iterator end () const { return  selected_.end () ; }

   void select (edm::Handle<reco::GsfElectronCollection>,
                const edm::Event&, 
                const edm::EventSetup&) ;
	
 private:	
   std::vector<reco::GsfElectronRef> selected_;
   double theTrackIsolCut_; 
   bool absolute_;
};


#endif
