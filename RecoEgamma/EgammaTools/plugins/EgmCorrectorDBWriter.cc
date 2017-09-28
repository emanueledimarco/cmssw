#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/EgammaObjects/interface/EgmCorrectorParameters.h"

class EgmCorrectorDBWriter : public edm::EDAnalyzer {
public:
  EgmCorrectorDBWriter(const edm::ParameterSet&);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override {}
  virtual void endJob() override {}
  ~EgmCorrectorDBWriter() {}
private:
  std::string txtfile;
  std::string payloadTag;
  std::string record;
  uint32_t startIOVRun;
};

// Constructor
EgmCorrectorDBWriter::EgmCorrectorDBWriter(const edm::ParameterSet& pSet)
{
  txtfile = pSet.getUntrackedParameter<std::string>("txtfile");
  payloadTag = pSet.getUntrackedParameter<std::string>("payloadTag");
  record = pSet.getUntrackedParameter<std::string>("record");
  startIOVRun = pSet.getUntrackedParameter<uint32_t>("startIOVRun");
}

// Begin Job
void EgmCorrectorDBWriter::beginJob() {
  std::cout << "Starting to import payload " << payloadTag << " from text files." << std::endl;
  try {
    std::cout << "Opened file " << txtfile << std::endl;
    // create the parameter object from file 
    EgmCorrectorParameters *payload = new EgmCorrectorParameters(txtfile);  

    // now write it into the DB
    std::cout << "Opening PoolDBOutputService" << std::endl;
    edm::Service<cond::service::PoolDBOutputService> s;
    if (s.isAvailable()) {
      if (s->isNewTagRequest(record)) 
        s->writeOne( payload, s->beginOfTime(), record);
      else
        s->writeOne( payload, startIOVRun, record);
      std::cout << "Wrote in CondDB payload label: " << record << std::endl;
    }
  }
  catch(edm::Exception ex) {
    std::cout << "Did not find the corrections file " << txtfile << std::endl;
  }
}

DEFINE_FWK_MODULE(EgmCorrectorDBWriter);
