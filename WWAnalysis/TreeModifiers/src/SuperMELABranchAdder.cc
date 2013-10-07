#include "TSystem.h"
#include "WWAnalysis/TreeModifiers/interface/SuperMELABranchAdder.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

SuperMELABranchAdder::SuperMELABranchAdder(const edm::ParameterSet& pset) : BranchAdder(pset)
{
   CalculatorEE = NULL;
   CalculatorEM = NULL;
   CalculatorMM = NULL;

   SuperMELAMass = pset.getParameter<double>("Mass");

   string DatasetConf = pset.getParameter<string>("datasetconf");
   if(DatasetConf == "2012")
      CMEnergy = 8;
   else
      CMEnergy = 7;
}

SuperMELABranchAdder::~SuperMELABranchAdder()
{
   if(CalculatorEE != NULL)
      delete CalculatorEE;
   if(CalculatorEM != NULL)
      delete CalculatorEM;
   if(CalculatorMM != NULL)
      delete CalculatorMM;

   CalculatorEE = NULL;
   CalculatorEM = NULL;
   CalculatorMM = NULL;
}

void SuperMELABranchAdder::initialize(TTree& tree)
{
   if(CalculatorEE != NULL)
      delete CalculatorEE;
   if(CalculatorEM != NULL)
      delete CalculatorEM;
   if(CalculatorMM != NULL)
      delete CalculatorMM;

   CalculatorEE = new SuperMELA(SuperMELAMass, "4e", CMEnergy);
   CalculatorEM = new SuperMELA(SuperMELAMass, "2e2mu", CMEnergy);
   CalculatorMM = new SuperMELA(SuperMELAMass, "4mu", CMEnergy);

   const char *CMSSWBase = gSystem->Getenv("CMSSW_BASE");
   char PathToCards[100] = "";
   sprintf(PathToCards, "%s/src/HiggsAnalysis/HZZ4L_CombinationPy/CreateDatacards/SM_inputs_%dTeV/",
      CMSSWBase, (int)CMEnergy);
   CalculatorEE->SetPathToCards(PathToCards);
   CalculatorEM->SetPathToCards(PathToCards);
   CalculatorMM->SetPathToCards(PathToCards);

   std::cout << PathToCards << std::endl;

   CalculatorEE->init();
   CalculatorEM->init();
   CalculatorMM->init();

   tree.SetBranchAddress("l1pt", &l1pt);
   tree.SetBranchAddress("l1eta", &l1eta);
   tree.SetBranchAddress("l1phi", &l1phi);
   tree.SetBranchAddress("l1pdgId", &l1pdgId);
   tree.SetBranchAddress("l2pt", &l2pt);
   tree.SetBranchAddress("l2eta", &l2eta);
   tree.SetBranchAddress("l2phi", &l2phi);
   tree.SetBranchAddress("l2pdgId", &l2pdgId);
   tree.SetBranchAddress("pho1pt", &pho1pt);
   tree.SetBranchAddress("pho1eta", &pho1eta);
   tree.SetBranchAddress("pho1phi", &pho1phi);

   tree.SetBranchAddress("l3pt", &l3pt);
   tree.SetBranchAddress("l3eta", &l3eta);
   tree.SetBranchAddress("l3phi", &l3phi);
   tree.SetBranchAddress("l3pdgId", &l3pdgId);
   tree.SetBranchAddress("l4pt", &l4pt);
   tree.SetBranchAddress("l4eta", &l4eta);
   tree.SetBranchAddress("l4phi", &l4phi);
   tree.SetBranchAddress("l4pdgId", &l4pdgId);
   tree.SetBranchAddress("pho2pt", &pho2pt);
   tree.SetBranchAddress("pho2eta", &pho2eta);
   tree.SetBranchAddress("pho2phi", &pho2phi);

   tree.SetBranchAddress("mass", &mass);
   tree.SetBranchAddress("z1mass", &z1mass);
   tree.SetBranchAddress("z2mass", &z2mass);

   tree.SetBranchAddress("melaPhi", &melaPhi);
   tree.SetBranchAddress("melaPhiStar1", &melaPhiStar1);
   tree.SetBranchAddress("melaCosTheta1", &melaCosTheta1);
   tree.SetBranchAddress("melaCosTheta2", &melaCosTheta2);
   tree.SetBranchAddress("melaCosThetaStar", &melaCosThetaStar);
   
   tree.SetBranchAddress("melaPt", &melaPt);
   tree.SetBranchAddress("melaY", &melaY);
   tree.SetBranchAddress("melaPtY", &melaPtY);
   tree.SetBranchAddress("melaLD", &melaLD);
}

float SuperMELABranchAdder::calculateValue(TTree &tree, int entry, int id, double &value)
{
   tree.GetEntry(entry);

   TLorentzVector lep1vec;
   TLorentzVector lep2vec;
   TLorentzVector lep3vec;
   TLorentzVector lep4vec;
   
   double mass1 = 0; 
   double mass2 = 0;
   double mass3 = 0;
   double mass4 = 0;
   
   if(abs(l1pdgId) == 11)
      mass1 = 0.51099892e-3;
   else if(abs(l1pdgId) == 13)
      mass1 = 105.658369e-3;
   if(abs(l2pdgId) == 11)
      mass2 = 0.51099892e-3;
   else if(abs(l2pdgId) == 13)
      mass2 = 105.658369e-3;
   if(abs(l3pdgId) == 11)
      mass3 = 0.51099892e-3;
   else if(abs(l3pdgId) == 13)
      mass3 = 105.658369e-3;
   if(abs(l4pdgId) == 11)
      mass4 = 0.51099892e-3;
   else if(abs(l4pdgId) == 13)
      mass4 = 105.658369e-3;

   if(l1pdgId > 0)
   {
      lep1vec.SetPtEtaPhiM(l1pt, l1eta, l1phi, mass1);
      lep2vec.SetPtEtaPhiM(l2pt, l2eta, l2phi, mass2);
   }
   else
   {
      lep1vec.SetPtEtaPhiM(l2pt, l2eta, l2phi, mass2);
      lep2vec.SetPtEtaPhiM(l1pt, l1eta, l1phi, mass1);
   }

   if(l3pdgId > 0)
   {
      lep3vec.SetPtEtaPhiM(l3pt, l3eta, l3phi, mass3);
      lep4vec.SetPtEtaPhiM(l4pt, l4eta, l4phi, mass4);
   }
   else
   {
      lep3vec.SetPtEtaPhiM(l4pt, l4eta, l4phi, mass4);
      lep4vec.SetPtEtaPhiM(l3pt, l3eta, l3phi, mass3);
   }

   // FSR recovery
   TLorentzVector pho1vec;
   TLorentzVector pho2vec;
   if(pho1pt > 0)
   {
      pho1vec.SetPtEtaPhiM(pho1pt, pho1eta, pho1phi, 0);
      if(reco::deltaR(lep1vec.Eta(), lep1vec.Phi(), pho1eta, pho1phi)
         < reco::deltaR(lep2vec.Eta(), lep2vec.Phi(), pho1eta, pho1phi))
         lep1vec += pho1vec;
      else
         lep2vec += pho1vec;
   }
   if(pho2pt > 0)
   {
      pho2vec.SetPtEtaPhiM(pho2pt, pho2eta, pho2phi, 0);
      if(reco::deltaR(lep3vec.Eta(), lep3vec.Phi(), pho2eta, pho2phi)
         < reco::deltaR(lep4vec.Eta(), lep4vec.Phi(), pho2eta, pho2phi))
         lep3vec += pho2vec;
      else
         lep4vec += pho2vec;
   }

   lep1vec[3] = sqrt(lep1vec.P() * lep1vec.P() + mass1 * mass1);
   lep2vec[3] = sqrt(lep2vec.P() * lep2vec.P() + mass2 * mass2);
   lep3vec[3] = sqrt(lep3vec.P() * lep3vec.P() + mass3 * mass3);
   lep4vec[3] = sqrt(lep4vec.P() * lep4vec.P() + mass4 * mass4);

   // Collecting info for super MELA
   double MZZ = (lep1vec + lep2vec + lep3vec + lep4vec).M();
   double M1 = (lep1vec + lep2vec).M();
   double M2 = (lep3vec + lep4vec).M();

   if(MZZ < 105 || MZZ > 140)
   {
      value = -1;
      return 0;
   }

   SuperMELA *Calculator = CalculatorEM;
   if(abs(l1pdgId) == 11 && abs(l3pdgId) == 11)
      Calculator = CalculatorEE;
   if(abs(l1pdgId) == 13 && abs(l3pdgId) == 13)
      Calculator = CalculatorMM;

   double PSuperSignal = 0;
   double PSuperBackground = 0;
   double PSuperMELA = 0;
   double PSignal = 0;
   double PBackground = 0;
   double PMELA = 0;

   Calculator->SetVerbosity(false);
   Calculator->RecalculateMELA(true);
   Calculator->SetDecayKinematics(z1mass, z2mass, melaCosThetaStar, melaCosTheta1, melaCosTheta2, melaPhi, melaPhiStar1);
   Calculator->computeKD(mass, false, PSuperMELA, PSuperSignal, PSuperBackground, PMELA, PSignal, PBackground);

   value = PSuperMELA;

   return 0;
}


