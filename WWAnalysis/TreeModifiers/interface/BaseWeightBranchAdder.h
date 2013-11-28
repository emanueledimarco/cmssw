#ifndef TreeModifiers_BaseWeightBranchAdder_h
#define TreeModifiers_BaseWeightBranchAdder_h


#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"
#include "WWAnalysis/TreeModifiers/interface/XSecProvider.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/** This class implements the calculateValue() method for 
 *  adding cross-section*BR*lumi weight.
 */


class BaseWeightBranchAdder : public BranchAdder {
 public:
  BaseWeightBranchAdder(const edm::ParameterSet& pset);
  virtual ~BaseWeightBranchAdder(){};


  /// add calculate the value to be added to the new branch
  virtual float calculateValue(TTree& tree,int entry, int id, double& value);

 private:
  //private methods here

  //datamembers
  XSecProvider xsecProvider_;
  bool is7TeV_;
};

#endif
