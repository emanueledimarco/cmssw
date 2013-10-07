#ifndef HWWKFACTORPRODUCER
#define HWWKFACTORPRODUCER

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "HepMC/WeightContainer.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWKFactorList.h"
#include "TH1D.h"
#include "TFile.h"

//
// class decleration
//

class HWWKFactorProducer : public edm::EDProducer {
   public:
      explicit HWWKFactorProducer(const edm::ParameterSet&);
      ~HWWKFactorProducer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
     
      edm::InputTag    genParticlesTag_;
      std::string inputFilename_;
      unsigned int  processID_;
      HWWKfactorList* pt_histo_;
      bool debug_;
};

#endif
