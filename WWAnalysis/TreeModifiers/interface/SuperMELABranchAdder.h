#ifndef TreeModifiers_SuperMELABranchAdder_h
#define TreeModifiers_SuperMELABranchAdder_h

#include "TH1F.h"

#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"
#include "WWAnalysis/TreeModifiers/interface/LeptSfProvider.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "ZZMatrixElement/SuperMELA/interface/SuperMELA.h"

#include "TMVA/Reader.h"


class SuperMELABranchAdder : public BranchAdder
{
public:
   SuperMELABranchAdder(const edm::ParameterSet& pset);
   virtual ~SuperMELABranchAdder();
   virtual float calculateValue(TTree& tree,int entry, int id, double& value);
   virtual void initialize(TTree& tree);

private:
   // Leptons
   float l1pt, l1eta, l1phi, l1pdgId;
   float l2pt, l2eta, l2phi, l2pdgId;
   float pho1pt, pho1eta, pho1phi;
   float l3pt, l3eta, l3phi, l3pdgId;
   float l4pt, l4eta, l4phi, l4pdgId;
   float pho2pt, pho2eta, pho2phi;
   float mass, z1mass, z2mass;

   // MELA variables
   float melaPhi, melaPhiStar1, melaCosTheta1, melaCosTheta2, melaCosThetaStar;
   float melaPt, melaY, melaPtY, melaLD;

   // SuperMELA stuff
   SuperMELA *CalculatorEE, *CalculatorEM, *CalculatorMM;
   double SuperMELAMass;
   double CMEnergy;
};

#endif
