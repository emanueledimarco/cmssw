// -*- C++ -*-
//
// Package:    HWWElectronSelector
// Class:      HWWElectronSelector
// 
/**\class HWWElectronSelector HWWElectronSelector.cc HiggsAnalysis/HWWElectronSelector/src/HWWElectronSelector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Thomas Punz
// Revised by: Emanuele Di Marco
//         Created:  Wed Jan 30 11:13:18 CET 2008
// $Id: HWWElectronSelector.cc,v 1.6 2010/02/08 14:06:22 ceballos Exp $
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
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWElectronSelector.h"

#include <string.h>

HWWElectronSelector::HWWElectronSelector(const edm::ParameterSet& iConfig)
{
  hOverECut_	      = iConfig.getParameter<double>("hOverECut");
  sigmaiEtaiEtaEBCut_ = iConfig.getParameter<double>("sigmaiEtaiEtaEBCut");
  sigmaiEtaiEtaEECut_ = iConfig.getParameter<double>("sigmaiEtaiEtaEECut");
  deltaphiinCut_      = iConfig.getParameter<double>("deltaphiinCut");
  deltaetainCut_      = iConfig.getParameter<double>("deltaetainCut");

}


HWWElectronSelector::~HWWElectronSelector()
{}


void
HWWElectronSelector::select (edm::Handle<reco::GsfElectronCollection> electrons,
                             const edm::Event& iEvent,
		             const edm::EventSetup& iEventSetup)
{

  using namespace edm;
  using namespace reco;
  selected_.clear();

  // Loop over electrons
  for (unsigned int i = 0; i < electrons->size(); i++) {	
    Ref<reco::GsfElectronCollection> electronRef(electrons,i);

    bool eleid = false;

    if(electronRef->isEB() == true){
      if(electronRef->hcalOverEcal()  < hOverECut_ &&
         electronRef->sigmaIetaIeta() < sigmaiEtaiEtaEBCut_ &&
	 electronRef->deltaPhiSuperClusterTrackAtVtx() < deltaphiinCut_ &&
	 electronRef->deltaEtaSuperClusterTrackAtVtx() < deltaetainCut_)
	  eleid = true;
    } else {
      if(electronRef->hcalOverEcal()  < hOverECut_ &&
         electronRef->sigmaIetaIeta() < sigmaiEtaiEtaEECut_ &&
	 electronRef->deltaPhiSuperClusterTrackAtVtx() < deltaphiinCut_ &&
	 electronRef->deltaEtaSuperClusterTrackAtVtx() < deltaetainCut_)
	  eleid = true;
    }
    
    if (eleid==true) selected_.push_back (electronRef);
    
  }

}
