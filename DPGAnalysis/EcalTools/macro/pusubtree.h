//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug  1 10:46:55 2014 by ROOT version 5.34/10
// from TTree ntp1/ntp1
// found on file: root://eoscms//store/cmst3/user/emanuele/ecal/amplireco/ntuples/rereco_photongun_pu25_0.root
//////////////////////////////////////////////////////////

#ifndef pusubtree_h
#define pusubtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <math.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class pusubtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runNumber;
   ULong64_t       eventNumber;
   Int_t           lumiBlock;
   Int_t           bunchCrossing;
   Int_t           orbitNumber;
   Double_t        beamSpotX;
   Double_t        beamSpotY;
   Double_t        beamSpotZ;
   Int_t           nBX;
   Int_t           nObsPU[16];   //[nBX]
   Float_t         nTruePU[16];   //[nBX]
   Int_t           bxPU[16];   //[nBX]
   Int_t           nEBSuperClusters;
   Int_t           nBCEBSuperClusters[100];   //[nEBSuperClusters]
   Int_t           nCrystalsEBSuperClusters[100];   //[nEBSuperClusters]
   Float_t         rawEnergyEBSuperClusters[100];   //[nEBSuperClusters]
   Float_t         energyEBSuperClusters[100];   //[nEBSuperClusters]
   Float_t         etaEBSuperClusters[100];   //[nEBSuperClusters]
   Float_t         thetaEBSuperClusters[100];   //[nEBSuperClusters]
   Float_t         phiEBSuperClusters[100];   //[nEBSuperClusters]
   Int_t           nEESuperClusters;
   Int_t           nBCEESuperClusters[500];   //[nEESuperClusters]
   Int_t           nCrystalsEESuperClusters[500];   //[nEESuperClusters]
   Float_t         rawEnergyEESuperClusters[500];   //[nEESuperClusters]
   Float_t         energyEESuperClusters[500];   //[nEESuperClusters]
   Float_t         etaEESuperClusters[500];   //[nEESuperClusters]
   Float_t         thetaEESuperClusters[500];   //[nEESuperClusters]
   Float_t         phiEESuperClusters[500];   //[nEESuperClusters]
   Int_t           nEBSuperClusters2;
   Int_t           nBCEBSuperClusters2[100];   //[nEBSuperClusters2]
   Int_t           nCrystalsEBSuperClusters2[100];   //[nEBSuperClusters2]
   Float_t         rawEnergyEBSuperClusters2[100];   //[nEBSuperClusters2]
   Float_t         energyEBSuperClusters2[100];   //[nEBSuperClusters2]
   Float_t         etaEBSuperClusters2[100];   //[nEBSuperClusters2]
   Float_t         thetaEBSuperClusters2[100];   //[nEBSuperClusters2]
   Float_t         phiEBSuperClusters2[100];   //[nEBSuperClusters2]
   Int_t           nEESuperClusters2;
   Int_t           nBCEESuperClusters2[500];   //[nEESuperClusters2]
   Int_t           nCrystalsEESuperClusters2[500];   //[nEESuperClusters2]
   Float_t         rawEnergyEESuperClusters2[500];   //[nEESuperClusters2]
   Float_t         energyEESuperClusters2[500];   //[nEESuperClusters2]
   Float_t         etaEESuperClusters2[500];   //[nEESuperClusters2]
   Float_t         thetaEESuperClusters2[500];   //[nEESuperClusters2]
   Float_t         phiEESuperClusters2[500];   //[nEESuperClusters2]
   Int_t           nEBSuperClusters3;
   Int_t           nBCEBSuperClusters3[100];   //[nEBSuperClusters3]
   Int_t           nCrystalsEBSuperClusters3[100];   //[nEBSuperClusters3]
   Float_t         rawEnergyEBSuperClusters3[100];   //[nEBSuperClusters3]
   Float_t         energyEBSuperClusters3[100];   //[nEBSuperClusters3]
   Float_t         etaEBSuperClusters3[100];   //[nEBSuperClusters3]
   Float_t         thetaEBSuperClusters3[100];   //[nEBSuperClusters3]
   Float_t         phiEBSuperClusters3[100];   //[nEBSuperClusters3]
   Int_t           nEESuperClusters3;
   Int_t           nBCEESuperClusters3[500];   //[nEESuperClusters3]
   Int_t           nCrystalsEESuperClusters3[500];   //[nEESuperClusters3]
   Float_t         rawEnergyEESuperClusters3[500];   //[nEESuperClusters3]
   Float_t         energyEESuperClusters3[500];   //[nEESuperClusters3]
   Float_t         etaEESuperClusters3[500];   //[nEESuperClusters3]
   Float_t         thetaEESuperClusters3[500];   //[nEESuperClusters3]
   Float_t         phiEESuperClusters3[500];   //[nEESuperClusters3]
   Int_t           nMc;
   Float_t         pMc[1];   //[nMc]
   Float_t         thetaMc[1];   //[nMc]
   Float_t         etaMc[1];   //[nMc]
   Float_t         phiMc[1];   //[nMc]
   Float_t         energyMc[1];   //[nMc]
   Float_t         vxMc[1];   //[nMc]
   Float_t         vyMc[1];   //[nMc]
   Float_t         vzMc[1];   //[nMc]
   Int_t           idMc[1];   //[nMc]
   Int_t           mothMc[1];   //[nMc]
   Int_t           statusMc[1];   //[nMc]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_beamSpotX;   //!
   TBranch        *b_beamSpotY;   //!
   TBranch        *b_beamSpotZ;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nObsPU;   //!
   TBranch        *b_nTruePU;   //!
   TBranch        *b_bxPU;   //!
   TBranch        *b_nEBSuperClusters;   //!
   TBranch        *b_nBCEBSuperClusters;   //!
   TBranch        *b_nCrystalsEBSuperClusters;   //!
   TBranch        *b_rawEnergyEBSuperClusters;   //!
   TBranch        *b_energyEBSuperClusters;   //!
   TBranch        *b_etaEBSuperClusters;   //!
   TBranch        *b_thetaEBSuperClusters;   //!
   TBranch        *b_phiEBSuperClusters;   //!
   TBranch        *b_nEESuperClusters;   //!
   TBranch        *b_nBCEESuperClusters;   //!
   TBranch        *b_nCrystalsEESuperClusters;   //!
   TBranch        *b_rawEnergyEESuperClusters;   //!
   TBranch        *b_energyEESuperClusters;   //!
   TBranch        *b_etaEESuperClusters;   //!
   TBranch        *b_thetaEESuperClusters;   //!
   TBranch        *b_phiEESuperClusters;   //!
   TBranch        *b_nEBSuperClusters2;   //!
   TBranch        *b_nBCEBSuperClusters2;   //!
   TBranch        *b_nCrystalsEBSuperClusters2;   //!
   TBranch        *b_rawEnergyEBSuperClusters2;   //!
   TBranch        *b_energyEBSuperClusters2;   //!
   TBranch        *b_etaEBSuperClusters2;   //!
   TBranch        *b_thetaEBSuperClusters2;   //!
   TBranch        *b_phiEBSuperClusters2;   //!
   TBranch        *b_nEESuperClusters2;   //!
   TBranch        *b_nBCEESuperClusters2;   //!
   TBranch        *b_nCrystalsEESuperClusters2;   //!
   TBranch        *b_rawEnergyEESuperClusters2;   //!
   TBranch        *b_energyEESuperClusters2;   //!
   TBranch        *b_etaEESuperClusters2;   //!
   TBranch        *b_thetaEESuperClusters2;   //!
   TBranch        *b_phiEESuperClusters2;   //!
   TBranch        *b_nEBSuperClusters3;   //!
   TBranch        *b_nBCEBSuperClusters3;   //!
   TBranch        *b_nCrystalsEBSuperClusters3;   //!
   TBranch        *b_rawEnergyEBSuperClusters3;   //!
   TBranch        *b_energyEBSuperClusters3;   //!
   TBranch        *b_etaEBSuperClusters3;   //!
   TBranch        *b_thetaEBSuperClusters3;   //!
   TBranch        *b_phiEBSuperClusters3;   //!
   TBranch        *b_nEESuperClusters3;   //!
   TBranch        *b_nBCEESuperClusters3;   //!
   TBranch        *b_nCrystalsEESuperClusters3;   //!
   TBranch        *b_rawEnergyEESuperClusters3;   //!
   TBranch        *b_energyEESuperClusters3;   //!
   TBranch        *b_etaEESuperClusters3;   //!
   TBranch        *b_thetaEESuperClusters3;   //!
   TBranch        *b_phiEESuperClusters3;   //!
   TBranch        *b_nMc;   //!
   TBranch        *b_pMc;   //!
   TBranch        *b_thetaMc;   //!
   TBranch        *b_etaMc;   //!
   TBranch        *b_phiMc;   //!
   TBranch        *b_energyMc;   //!
   TBranch        *b_vxMc;   //!
   TBranch        *b_vyMc;   //!
   TBranch        *b_vzMc;   //!
   TBranch        *b_idMc;   //!
   TBranch        *b_mothMc;   //!
   TBranch        *b_statusMc;   //!

   pusubtree(TTree *tree=0);
   virtual ~pusubtree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char *outputfilename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool     mcmatch(TVector3 scdir, TVector3 truep, int charge);
};

#endif

#ifdef pusubtree_cxx
pusubtree::pusubtree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//store/cmst3/user/emanuele/ecal/amplireco/ntuples/rereco_photongun_pu25_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//store/cmst3/user/emanuele/ecal/amplireco/ntuples/rereco_photongun_pu25_0.root");
      }
      f->GetObject("ntp1",tree);

   }
   Init(tree);
}

pusubtree::~pusubtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pusubtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pusubtree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pusubtree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   fChain->SetBranchAddress("nObsPU", nObsPU, &b_nObsPU);
   fChain->SetBranchAddress("nTruePU", nTruePU, &b_nTruePU);
   fChain->SetBranchAddress("bxPU", bxPU, &b_bxPU);
   fChain->SetBranchAddress("nEBSuperClusters", &nEBSuperClusters, &b_nEBSuperClusters);
   fChain->SetBranchAddress("nBCEBSuperClusters", nBCEBSuperClusters, &b_nBCEBSuperClusters);
   fChain->SetBranchAddress("nCrystalsEBSuperClusters", nCrystalsEBSuperClusters, &b_nCrystalsEBSuperClusters);
   fChain->SetBranchAddress("rawEnergyEBSuperClusters", rawEnergyEBSuperClusters, &b_rawEnergyEBSuperClusters);
   fChain->SetBranchAddress("energyEBSuperClusters", energyEBSuperClusters, &b_energyEBSuperClusters);
   fChain->SetBranchAddress("etaEBSuperClusters", etaEBSuperClusters, &b_etaEBSuperClusters);
   fChain->SetBranchAddress("thetaEBSuperClusters", thetaEBSuperClusters, &b_thetaEBSuperClusters);
   fChain->SetBranchAddress("phiEBSuperClusters", phiEBSuperClusters, &b_phiEBSuperClusters);
   fChain->SetBranchAddress("nEESuperClusters", &nEESuperClusters, &b_nEESuperClusters);
   fChain->SetBranchAddress("nBCEESuperClusters", nBCEESuperClusters, &b_nBCEESuperClusters);
   fChain->SetBranchAddress("nCrystalsEESuperClusters", nCrystalsEESuperClusters, &b_nCrystalsEESuperClusters);
   fChain->SetBranchAddress("rawEnergyEESuperClusters", rawEnergyEESuperClusters, &b_rawEnergyEESuperClusters);
   fChain->SetBranchAddress("energyEESuperClusters", energyEESuperClusters, &b_energyEESuperClusters);
   fChain->SetBranchAddress("etaEESuperClusters", etaEESuperClusters, &b_etaEESuperClusters);
   fChain->SetBranchAddress("thetaEESuperClusters", thetaEESuperClusters, &b_thetaEESuperClusters);
   fChain->SetBranchAddress("phiEESuperClusters", phiEESuperClusters, &b_phiEESuperClusters);
   fChain->SetBranchAddress("nEBSuperClusters2", &nEBSuperClusters2, &b_nEBSuperClusters2);
   fChain->SetBranchAddress("nBCEBSuperClusters2", nBCEBSuperClusters2, &b_nBCEBSuperClusters2);
   fChain->SetBranchAddress("nCrystalsEBSuperClusters2", nCrystalsEBSuperClusters2, &b_nCrystalsEBSuperClusters2);
   fChain->SetBranchAddress("rawEnergyEBSuperClusters2", rawEnergyEBSuperClusters2, &b_rawEnergyEBSuperClusters2);
   fChain->SetBranchAddress("energyEBSuperClusters2", energyEBSuperClusters2, &b_energyEBSuperClusters2);
   fChain->SetBranchAddress("etaEBSuperClusters2", etaEBSuperClusters2, &b_etaEBSuperClusters2);
   fChain->SetBranchAddress("thetaEBSuperClusters2", thetaEBSuperClusters2, &b_thetaEBSuperClusters2);
   fChain->SetBranchAddress("phiEBSuperClusters2", phiEBSuperClusters2, &b_phiEBSuperClusters2);
   fChain->SetBranchAddress("nEESuperClusters2", &nEESuperClusters2, &b_nEESuperClusters2);
   fChain->SetBranchAddress("nBCEESuperClusters2", nBCEESuperClusters2, &b_nBCEESuperClusters2);
   fChain->SetBranchAddress("nCrystalsEESuperClusters2", nCrystalsEESuperClusters2, &b_nCrystalsEESuperClusters2);
   fChain->SetBranchAddress("rawEnergyEESuperClusters2", rawEnergyEESuperClusters2, &b_rawEnergyEESuperClusters2);
   fChain->SetBranchAddress("energyEESuperClusters2", energyEESuperClusters2, &b_energyEESuperClusters2);
   fChain->SetBranchAddress("etaEESuperClusters2", etaEESuperClusters2, &b_etaEESuperClusters2);
   fChain->SetBranchAddress("thetaEESuperClusters2", thetaEESuperClusters2, &b_thetaEESuperClusters2);
   fChain->SetBranchAddress("phiEESuperClusters2", phiEESuperClusters2, &b_phiEESuperClusters2);
   fChain->SetBranchAddress("nEBSuperClusters3", &nEBSuperClusters3, &b_nEBSuperClusters3);
   fChain->SetBranchAddress("nBCEBSuperClusters3", nBCEBSuperClusters3, &b_nBCEBSuperClusters3);
   fChain->SetBranchAddress("nCrystalsEBSuperClusters3", nCrystalsEBSuperClusters3, &b_nCrystalsEBSuperClusters3);
   fChain->SetBranchAddress("rawEnergyEBSuperClusters3", rawEnergyEBSuperClusters3, &b_rawEnergyEBSuperClusters3);
   fChain->SetBranchAddress("energyEBSuperClusters3", energyEBSuperClusters3, &b_energyEBSuperClusters3);
   fChain->SetBranchAddress("etaEBSuperClusters3", etaEBSuperClusters3, &b_etaEBSuperClusters3);
   fChain->SetBranchAddress("thetaEBSuperClusters3", thetaEBSuperClusters3, &b_thetaEBSuperClusters3);
   fChain->SetBranchAddress("phiEBSuperClusters3", phiEBSuperClusters3, &b_phiEBSuperClusters3);
   fChain->SetBranchAddress("nEESuperClusters3", &nEESuperClusters3, &b_nEESuperClusters3);
   fChain->SetBranchAddress("nBCEESuperClusters3", nBCEESuperClusters3, &b_nBCEESuperClusters3);
   fChain->SetBranchAddress("nCrystalsEESuperClusters3", nCrystalsEESuperClusters3, &b_nCrystalsEESuperClusters3);
   fChain->SetBranchAddress("rawEnergyEESuperClusters3", rawEnergyEESuperClusters3, &b_rawEnergyEESuperClusters3);
   fChain->SetBranchAddress("energyEESuperClusters3", energyEESuperClusters3, &b_energyEESuperClusters3);
   fChain->SetBranchAddress("etaEESuperClusters3", etaEESuperClusters3, &b_etaEESuperClusters3);
   fChain->SetBranchAddress("thetaEESuperClusters3", thetaEESuperClusters3, &b_thetaEESuperClusters3);
   fChain->SetBranchAddress("phiEESuperClusters3", phiEESuperClusters3, &b_phiEESuperClusters3);
   fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
   fChain->SetBranchAddress("pMc", pMc, &b_pMc);
   fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
   fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
   fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
   fChain->SetBranchAddress("energyMc", energyMc, &b_energyMc);
   fChain->SetBranchAddress("vxMc", vxMc, &b_vxMc);
   fChain->SetBranchAddress("vyMc", vyMc, &b_vyMc);
   fChain->SetBranchAddress("vzMc", vzMc, &b_vzMc);
   fChain->SetBranchAddress("idMc", idMc, &b_idMc);
   fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
   fChain->SetBranchAddress("statusMc", statusMc, &b_statusMc);
   Notify();
}

Bool_t pusubtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pusubtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pusubtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
bool pusubtree::mcmatch(TVector3 scdir, TVector3 truep, int charge) {
  // for the EB, back-propagate with the magnetic field, since it is simple.
  // for the EE, only match in deltaEta
  
  float abseta = fabs(scdir.Eta());
  float phicurv = atan2(1.24*3.8*charge,4.*scdir.Pt());
  
  float deltaeta = scdir.Eta() - truep.Eta();  
  if(abseta<1.479) {
    float deltaphi = scdir.DeltaPhi(truep) - phicurv;
    return sqrt(deltaeta*deltaeta + deltaphi*deltaphi) < 0.1;
  } else {
    return fabs(deltaeta) < 0.1;
  }
}
#endif // #ifdef pusubtree_cxx
