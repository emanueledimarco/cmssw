#ifndef TreeModifiers_KinematicBDTBranchAdder_h
#define TreeModifiers_KinematicBDTBranchAdder_h

#include "TH1F.h"

#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"
#include "WWAnalysis/TreeModifiers/interface/LeptSfProvider.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

/** This class implements the calculateValue() method for 
 *  adding lepton scale factors.
 */

#include "TMVA/Reader.h"


class KinematicBDTBranchAdder : public BranchAdder {
 public:
 KinematicBDTBranchAdder(const edm::ParameterSet& pset) ;
  
  virtual ~KinematicBDTBranchAdder(){
    //delete histoWeights_;
  };

  /// add calculate the value to be added to the new branch
  virtual float calculateValue(TTree& tree,int entry, int id, double& value) ;

  
  virtual void initialize(TTree& tree);


 private:
  //add datamembers and private methods here
  
    //leptons
    float l1pt, l1eta, l1phi, l1pdgId;
    float l2pt, l2eta, l2phi, l2pdgId;
    float pho1pt, pho1eta, pho1phi;
    float l3pt, l3eta, l3phi, l3pdgId;
    float l4pt, l4eta, l4phi, l4pdgId;
    float pho2pt, pho2eta, pho2phi;

    //BDT input variables
    TMVA::Reader *BDTReader;
    std::string weightfile;
    float MVAInputVar_costheta1;
    float MVAInputVar_costheta2;
    float MVAInputVar_costhetastar;
    float MVAInputVar_Phi;
    float MVAInputVar_Phi1;
    float MVAInputVar_mZ1;
    float MVAInputVar_mZ2;
    float MVAInputVar_ZZpt;
    float MVAInputVar_zzdotz1;
    float MVAInputVar_zzdotz2;
    float MVAInputVar_ZZptZ1ptCosDphi;
    float MVAInputVar_ZZptZ2ptCosDphi;
    float MVAInputVar_reduced_Z1pt;
    float MVAInputVar_reduced_Z2pt;
    float MVAInputVar_ZZy;
    float MVAInputVar_m4l;

};

#endif
