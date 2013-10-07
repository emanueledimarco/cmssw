#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


class ChargedMetPrinter : public edm::EDAnalyzer {
   public:
      explicit ChargedMetPrinter(const edm::ParameterSet&);
      ~ChargedMetPrinter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
};

void ChargedMetPrinter::analyze(const edm::Event& evt, const edm::EventSetup& es) {

    edm::Handle<edm::ValueMap<reco::PFMET> > vmH;
    edm::Handle<reco::VertexCollection> vtxH;

    evt.getByLabel("chargedMetProducer",vmH);
    evt.getByLabel(edm::InputTag("goodPrimaryVertices","","Yield"),vtxH);

    //I am making an assumption that the first vertex was one of the vertices
    //used to build the charged mets.  But, actually, I only built the charged
    //met for the vertices closest to each of the leptons in the event, so you 
    //should find the closest vertex to your lepton, then dereference the 
    //value map.  It will return 0 if the vertex wasnt used for calculating the MET
    std::cout << (*vmH)[reco::VertexRef(vtxH,0)].pt() << std::endl;

}


ChargedMetPrinter::ChargedMetPrinter(const edm::ParameterSet& iConfig) { } 
ChargedMetPrinter::~ChargedMetPrinter() { }
void ChargedMetPrinter::beginJob() { }
void ChargedMetPrinter::endJob() { } 
void ChargedMetPrinter::beginRun(edm::Run const&, edm::EventSetup const&) { } 
void ChargedMetPrinter::endRun(edm::Run const&, edm::EventSetup const&) { }
void ChargedMetPrinter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { } 
void ChargedMetPrinter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }
void ChargedMetPrinter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ChargedMetPrinter);
