#ifndef TreeModifiers_BranchAdder_h
#define TreeModifiers_BranchAdder_h

#include "TTree.h"
#include <string>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>

/** This class defines the interface for a generic BranchAdder 
 * that add a new branch to a flat TTree provided in input.
 * The value stored in the new branch is calculated from 
 * the branches already present in the tree or from some
 * other external inputs like xml, text files, etc.
 */

class BranchAdder {
 public:
 BranchAdder(const edm::ParameterSet& pset)  {
    branchName_     = pset.getParameter<std::string>("newBranchName"); 
    treesToProcess_ = pset.getParameter<std::vector<std::string> >("treesToProcess");
    
  }
  virtual ~BranchAdder(){};


  /// add new branch to input tree
  int addBranch(TTree& tree, int id=-99);


  /// returns the names of the directories of the trees that have to be processed.
  std::vector<std::string> treesToProcess(){return treesToProcess_;}

  std::string branchName() const {return branchName_;}

  
  /// method were the value to be added to the tree branch is actually calculated
  virtual float calculateValue(TTree& tree, int entry, int id, double& value)  {return -999;}

  
 private:
  std::string branchName_;
  std::vector<std::string> treesToProcess_;

 protected:
  /// does initializations like mapping between tree branches and local variable of the branchAdder 
  virtual void initialize(TTree& tree) {
    //std::cout << "initialize in Base" << std::endl;
  }

};

#endif
