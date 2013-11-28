#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include<string>

class SkimEventJECInitializer : public edm::EDAnalyzer {
   public:
      explicit SkimEventJECInitializer(const edm::ParameterSet&);
      ~SkimEventJECInitializer();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
};

SkimEventJECInitializer::SkimEventJECInitializer(const edm::ParameterSet& iConfig) { 
    reco::SkimEvent::setupJEC(
        iConfig.getParameter<std::string>("l2File"),
        iConfig.getParameter<std::string>("l3File"),
        iConfig.getParameter<std::string>("resFile")
    );
}


void SkimEventJECInitializer::beginJob() { }
void SkimEventJECInitializer::endJob() { }
SkimEventJECInitializer::~SkimEventJECInitializer() { }
void SkimEventJECInitializer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { }
DEFINE_FWK_MODULE(SkimEventJECInitializer);
