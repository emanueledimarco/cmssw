// -*- C++ -*-
//
// Package:    HWWKFactorProducer
// Class:      HWWKFactorProducer
//
/**\class HWWKFactorProducer HWWKFactorProducer.cc HiggsAnalysis/HWWKFactorProducer/src/HWWKFactorProducer.cc
 
Description: <one line class summary>
 
Implementation:
<Notes on implementation>
*/
//
// Original Author:  Joanna Weng
//         Created:  Fri Feb  1 15:30:42 CET 2008
// $Id: HWWKFactorProducer.cc,v 1.8 2011/05/17 09:15:33 mangano Exp $
//
//

// system include files
#include <memory>
#include <iostream>

// user include files
#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWKFactorProducer.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "HepMC/WeightContainer.h"
//#include "HepMC/GenEvent.h"
//#include "HepMC/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TH1D.h"
#include "TFile.h"



HWWKFactorProducer::HWWKFactorProducer(const edm::ParameterSet& iConfig): pt_histo_(0),debug_(0)
{
  produces<double>();
  genParticlesTag_ = iConfig.getParameter<edm::InputTag>("genParticlesTag");
  inputFilename_=iConfig.getUntrackedParameter<std::string>("inputFilename","dummy.root");  
  processID_ = iConfig.getUntrackedParameter<int>("ProcessID",0);  
  debug_ = iConfig.getUntrackedParameter<bool>("Debug",0);
  edm::FileInPath path_inputFilename(inputFilename_.c_str());
  pt_histo_ = new  HWWKfactorList("KFactorList", path_inputFilename.fullPath().c_str() );

}


HWWKFactorProducer::~HWWKFactorProducer()
{}

// ------------ method called to for each event  ------------
void HWWKFactorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)

{
  using namespace edm;

  // Event weight to be stored in event
  std::auto_ptr<double> pweight(new double);
  *pweight=1.;

  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticlesTag_, genParticles);


   // use event info product
  Handle<GenEventInfoProduct> hEvtInfo;
  iEvent.getByLabel("generator", hEvtInfo);

  if (debug_)std::cout << " Process Id " << hEvtInfo->signalProcessID()   << std::endl; 
  // Gluon Fusion found
  if ( processID_ == hEvtInfo->signalProcessID() ){
    // look for Higgs Boson and determine its pt
    reco::GenParticleCollection::const_iterator higgs;
    for(reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++){
      if(abs(mcIter->pdgId())==25)
	{higgs=mcIter; continue;}
    }
    //Print the complete table
    if (debug_) std::cout << (*pt_histo_);
    math::XYZTLorentzVector tot_momentum(higgs->p4());

    //old implemenation
    //calculate bin size
    //double binsize = (pt_histo_->GetXaxis()->GetXmax()-pt_histo_->GetXaxis()->GetXmin())/pt_histo_->GetNbinsX();
    //double higgspt = tot_momentum.pt();
    //int bin = 0;
    // underflow protection: use underflow entry
    //if(higgspt >= pt_histo_->GetXaxis()->GetXmin()){
    //  bin = Int_t((higgspt-pt_histo_->GetXaxis()->GetXmin())/binsize) + 1;
    //}
    // overflow protection: use overflow entry
    //if(bin > pt_histo_->GetNbinsX()) bin=pt_histo_->GetNbinsX()+1;
    //if (debug_){
    //std::cout <<" Bin Size "<< binsize <<std::endl;
    //std::cout <<" Higgs Pt "<< higgspt <<std::endl;
    //std::cout <<" Bin  "<< bin <<std::endl;
    //std::cout <<" KFactor "<<   pt_histo_->GetBinContent(bin) <<std::endl;
    //}
    // get KFactor
    //old implementation
    //*pweight=  pt_histo_->GetBinContent(bin);

    //new implementation
    if (debug_){
    std::cout <<" Bin Size "<< "it's variable" <<std::endl;
    std::cout <<" Higgs Pt "<< tot_momentum.pt() <<std::endl;
    std::cout <<" Bin  "<< pt_histo_->GetBinFromPt(tot_momentum.pt()) <<std::endl;
    std::cout <<" KFactor "<<   pt_histo_->GetWeightFromPt(tot_momentum.pt()) <<std::endl;
    }
    *pweight=  pt_histo_->GetWeightFromPt(tot_momentum.pt());
  }
  iEvent.put(pweight);
}

// ------------ method called once each job just before starting event loop  ------------
void
HWWKFactorProducer::beginJob(const edm::EventSetup&)
{}

// ------------ method called once each job just after ending the event loop  ------------
void
HWWKFactorProducer::endJob()
{
  if (pt_histo_) delete(pt_histo_);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HWWKFactorProducer);
