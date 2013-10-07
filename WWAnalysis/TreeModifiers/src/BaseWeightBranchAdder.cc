#include "WWAnalysis/TreeModifiers/interface/BaseWeightBranchAdder.h"


BaseWeightBranchAdder::BaseWeightBranchAdder(const edm::ParameterSet& pset) : BranchAdder(pset){
  std::string datasetConf_     = pset.getParameter<std::string>("datasetConf"); 
  if(datasetConf_ == "7TeV") is7TeV_ = true;
  else if(datasetConf_ == "8TeV") is7TeV_ = false;
  else{
    std::cout << "WARNING: you didn't specify the \"datasetConf\" parameter. Using default: 8TeV" 
	      << std::endl;
    is7TeV_ = false;
  }

  xsecProvider_.initXsec();
}

float BaseWeightBranchAdder::calculateValue(TTree& tree,int k,  int id,
					    double& value)  {
  using namespace std;

  tree.GetEntry(k);
  
  /*
  map<int,double>::iterator it =  baseWeightMap_.find(id);
  if(it != baseWeightMap_.end())
    value = it->second;
  else{
    cout << "WARNING: baseWeight is not available for id=" << id << " and so it is set to 1." << endl;
  }
  */

  //TO DO: there is no check that the maps has been filled for that specific ID. Add a safety check.
  value = is7TeV_ ? xsecProvider_.getXSec7TeVByID(id)/xsecProvider_.getEvt7TeV(id) :
    xsecProvider_.getXSec8TeVByID(id)/xsecProvider_.getEvt8TeV(id) ;
  return 0;
}


