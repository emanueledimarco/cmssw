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
  std::string name;
  std::string reco;
  std::string version;
  std::string object;
  std::string path;
  std::string payloadTag;
  std::string record;
};

// Constructor
EgmCorrectorDBWriter::EgmCorrectorDBWriter(const edm::ParameterSet& pSet)
{
  name = pSet.getUntrackedParameter<std::string>("name");
  reco = pSet.getUntrackedParameter<std::string>("reco");
  version = pSet.getUntrackedParameter<std::string>("version");
  object = pSet.getUntrackedParameter<std::string>("object");
  path   = pSet.getUntrackedParameter<std::string>("path");
  payloadTag = "EgmCorrectorParameters_"+object+"_"+name+"_"+reco+"_"+version;
  record = pSet.getUntrackedParameter<std::string>("record");
}

// Begin Job
void EgmCorrectorDBWriter::beginJob() {
  std::cout << "Starting to import payload " << payloadTag << " from text files." << std::endl;
  std::string inputTxtFile = path+payloadTag+".txt";
  try {
    edm::FileInPath fip(inputTxtFile);
    std::cout << "Opened file " << inputTxtFile << std::endl;
    // create the parameter object from file 
    EgmCorrectorParameters *payload = new EgmCorrectorParameters(fip.fullPath());  

    // now write it into the DB
    std::cout << "Opening PoolDBOutputService" << std::endl;
    edm::Service<cond::service::PoolDBOutputService> s;
    if (s.isAvailable()) {
      s->writeOne( payload, s->beginOfTime(), record);
      std::cout << "Wrote in CondDB payload label: " << record << std::endl;
    }
  }
  catch(edm::Exception ex) {
    std::cout << "Did not find the corrections file " << inputTxtFile << std::endl;
  }
}

DEFINE_FWK_MODULE(EgmCorrectorDBWriter);
