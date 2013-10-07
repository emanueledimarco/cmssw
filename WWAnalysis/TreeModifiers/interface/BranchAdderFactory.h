#ifndef TreeModifiers_BranchAdderFactory_h
#define TreeModifiers_BranchAdderFactory_h

#include "TTree.h"
#include <string>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"

#include <iostream>
#include <string>

/** 
 *
 */

class BranchAdderFactory {
 public:
 BranchAdderFactory()  {}
  virtual ~BranchAdderFactory();


  /// accessor to the branchAdder;
  BranchAdder* getBranchAdder(std::string adderName, const edm::ParameterSet& pset);
  
  
 private:
  void createBranchAdder(std::string adderName, const edm::ParameterSet& pset);

  //data member
  std::vector<BranchAdder*> branchAdders_;


};

#endif
