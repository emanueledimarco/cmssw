#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DPGAnalysis/EcalTools/interface/CmsEcalRecHitFiller.h"
#include "DPGAnalysis/EcalTools/plugins/RecHitDumper.h"

using namespace edm;
using namespace reco;

RecHitDumper::RecHitDumper(const edm::ParameterSet& iConfig)
{
  
  nameFile_      = iConfig.getUntrackedParameter<std::string>("nameFile", "RootOutput.root");
  nameTree_      = iConfig.getUntrackedParameter<std::string>("nameTree", "BaseTree");
  // ECAL rechits collections
  ecalBarrelRecHits_ = iConfig.getParameter<edm::InputTag>("EBRecHits");
  ecalEndcapRecHits_ = iConfig.getParameter<edm::InputTag>("EERecHits");
}

RecHitDumper::~RecHitDumper() { }

// ------------ method called to for each event  ------------
void RecHitDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 

  CmsEcalRecHitFiller ebFill(tree_);
  std::string prefix("");
  std::string suffix("EBRecHits");
  ebFill.writeCollectionToTree(ecalBarrelRecHits_, iEvent, iSetup, prefix, suffix, false);

  CmsEcalRecHitFiller eeFill(tree_);
  suffix = std::string("EERecHits");
  ebFill.writeCollectionToTree(ecalEndcapRecHits_, iEvent, iSetup, prefix, suffix, false);

  tree_->dumpData();
  

}

// ------------ method called once each job just before starting event loop  ------------
void RecHitDumper::beginJob() {
  
  fileOut_ = TFile::Open(nameFile_.c_str(), "RECREATE");

  tree_  = new  CmsTree(nameTree_.c_str(),nameTree_.c_str());
}

void RecHitDumper::beginRun( const Run & iRun, const EventSetup & iSetup ) { }



// ------------ method called once each job just after ending the event loop  ------------
void  RecHitDumper::endJob() {

  fileOut_->cd();

  TTree* treeEventsOut = tree_->getTree();
  treeEventsOut->Write();

  fileOut_->Close();

}

DEFINE_FWK_MODULE (RecHitDumper);

