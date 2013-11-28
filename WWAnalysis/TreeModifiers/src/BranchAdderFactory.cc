#include "WWAnalysis/TreeModifiers/interface/BranchAdderFactory.h"
#include <iostream>


//--- HERE INCLUDE THE CONCRETE BRANCHADDERS ---
#include "WWAnalysis/TreeModifiers/interface/BaseWeightBranchAdder.h"
#include "WWAnalysis/TreeModifiers/interface/PileupBranchAdder.h"
#include "WWAnalysis/TreeModifiers/interface/LeptSfBranchAdder.h"
#include "WWAnalysis/TreeModifiers/interface/LeptPRBranchAdder.h"
//#include "WWAnalysis/TreeModifiers/interface/myBranchAdder.h"


BranchAdderFactory::~BranchAdderFactory(){
  for(unsigned int i=0; i<branchAdders_.size();++i){
    delete branchAdders_[i];
  }

}

BranchAdder* 
BranchAdderFactory::getBranchAdder(std::string adderName, const edm::ParameterSet& pset){
  if(1) //here one should add a check to avoid to create the same branchAdder twice
    createBranchAdder(adderName, pset);    
  return branchAdders_.back();

}


void BranchAdderFactory::createBranchAdder(std::string adderName, const edm::ParameterSet& pset){
  using namespace std;
  if(adderName == "BaseWeightBranchAdder") branchAdders_.push_back(new BaseWeightBranchAdder(pset));
  else if(adderName == "PileupBranchAdder") branchAdders_.push_back(new PileupBranchAdder(pset));
  else if(adderName == "LeptSfBranchAdder") branchAdders_.push_back(new LeptSfBranchAdder(pset));
  else if(adderName == "LeptPRBranchAdder") branchAdders_.push_back(new LeptPRBranchAdder(pset));
  else {
    cout << "WARNING: branchAdderName \"" << adderName << "\" is not known. Instanciating default class."
	 << endl;
    branchAdders_.push_back(new BranchAdder(pset));
  }
    

}
