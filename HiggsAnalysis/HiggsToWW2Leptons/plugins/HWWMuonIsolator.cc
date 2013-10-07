// -*- C++ -*-
//
// Package:    MuonIsolatorProducer
// Class:      MuonIsolatorProducer
// 
/**\class MuonIsolatorProducer MuonIsolatorProducer.cc MyProducer/MuonIsolatorProducer/src/MuonIsolatorProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Thomas Punz
//         Created:  Fri Jan 25 17:54:06 CET 2008
// $Id: HWWMuonIsolator.cc,v 1.3 2008/09/24 13:26:17 emanuele Exp $
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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWMuonIsolator.h"


HWWMuonIsolator::HWWMuonIsolator(const edm::ParameterSet& iConfig)
{
  using namespace edm;  

  theTrackIsolCut_	     = iConfig.getParameter<double>("trackIsolCut");
  absolute_                  = iConfig.getParameter<bool>("absolute");

}

HWWMuonIsolator::~HWWMuonIsolator()
{
}

void
HWWMuonIsolator::select(edm::Handle<reco::MuonCollection> muons,
                        const edm::Event& iEvent,
		        const edm::EventSetup& iEventSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  selected_.clear();
  Handle<RefVector<MuonCollection> >muonsRef;

  for (unsigned i = 0;i<muons->size(); i++){ 
    Ref<MuonCollection>muonRAWRef(muons,i);

    double pt = (*muons)[i].pt();
    double sumPt = (*muons)[i].isolationR03().sumPt;
    
    if ( !absolute_ && pt != 0 ) sumPt = sumPt/pt;
    
    if ( sumPt < theTrackIsolCut_ ) selected_.push_back(muonRAWRef);

  } // loop over muons
}
//define this as a plug-in
//DEFINE_FWK_MODULE(HWWMuonIsolator);
