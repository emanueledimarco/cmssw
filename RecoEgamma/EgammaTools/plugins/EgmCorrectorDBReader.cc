// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/EgammaObjects/interface/EgmCorrectorParameters.h"
#include "CondFormats/DataRecord/interface/EgmCorrectorParametersRcd.h"

class EgmCorrectorDBReader : public edm::EDAnalyzer {
public:
  explicit EgmCorrectorDBReader(const edm::ParameterSet&);
  ~EgmCorrectorDBReader();
  
  
private:
  virtual void beginJob() override ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override ;
 
  std::string mRecord,mGlobalTag;
  bool mCreateTextFile,mPrintScreen;
};


EgmCorrectorDBReader::EgmCorrectorDBReader(const edm::ParameterSet& iConfig)
{
  mRecord         = iConfig.getUntrackedParameter<std::string>("record");
  mGlobalTag      = iConfig.getUntrackedParameter<std::string>("globalTag");  
  mPrintScreen    = iConfig.getUntrackedParameter<bool>("printScreen");
  mCreateTextFile = iConfig.getUntrackedParameter<bool>("createTextFile");
}


EgmCorrectorDBReader::~EgmCorrectorDBReader()
{
 
}

void EgmCorrectorDBReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<EgmCorrectorParameters> EgmCorParams;
  std::cout <<"Inspecting EGM corrections record with name: "<< mRecord <<std::endl;
  iSetup.get<EgmCorrectorParametersRcd>().get(EgmCorParams);
  std::cout<<"-------------------------------------------------" << std::endl;
  std::cout<<"Processing record = " << mRecord << std::endl;
  if (mCreateTextFile) {
    std::cout<<"Creating txt file: "<<mGlobalTag+"_"+mRecord+".txt"<<std::endl;
    EgmCorParams->printFile(mGlobalTag+"_"+mRecord+".txt");
  }
  if (mPrintScreen)
    EgmCorParams->printScreen();
}

void 
EgmCorrectorDBReader::beginJob()
{
}

void 
EgmCorrectorDBReader::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EgmCorrectorDBReader);
