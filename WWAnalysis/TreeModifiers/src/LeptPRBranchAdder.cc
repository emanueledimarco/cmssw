#include "WWAnalysis/TreeModifiers/interface/LeptPRBranchAdder.h"


LeptPRBranchAdder::LeptPRBranchAdder(const edm::ParameterSet& pset) : BranchAdder(pset){
  std::string datasetConf_     = pset.getParameter<std::string>("datasetConf"); 
  bool is2011(false);
  if(datasetConf_ == "2011") is2011 = true;
  else if(datasetConf_ == "2012") is2011 = false;
  else{
    std::cout << "WARNING: you didn't specify the \"datasetConf\" parameter. Using default: 2012" 
	      << std::endl;
    is2011 = false;
  }
  sfProvider_.init(is2011);
  
}

void LeptPRBranchAdder::initialize(TTree& tree) {
  tree.SetBranchAddress("l1pt",&pt1_);tree.SetBranchAddress("l1eta",&eta1_);tree.SetBranchAddress("l1pdgId",&id1_);
  tree.SetBranchAddress("l2pt",&pt2_);tree.SetBranchAddress("l2eta",&eta2_);tree.SetBranchAddress("l2pdgId",&id2_);
  tree.SetBranchAddress("l3pt",&pt3_);tree.SetBranchAddress("l3eta",&eta3_);tree.SetBranchAddress("l3pdgId",&id3_);
  tree.SetBranchAddress("l4pt",&pt4_);tree.SetBranchAddress("l4eta",&eta4_);tree.SetBranchAddress("l4pdgId",&id4_);
}


float LeptPRBranchAdder::calculateValue(TTree& tree,int entry, int id, double& value)  {
  tree.GetEntry(entry);
  value = sfProvider_.getPR(pt1_,eta1_,id1_)*
    sfProvider_.getPR(pt2_,eta2_,id2_)*
    sfProvider_.getPR(pt3_,eta3_,id3_)*
    sfProvider_.getPR(pt4_,eta4_,id4_);


  return 0;
}


