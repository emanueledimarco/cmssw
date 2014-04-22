//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 28 10:01:33 2014 by ROOT version 5.34/03
// from TTree ntp0/ntp0
// found on file: dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/emanuele/ECAL_612_V1/stdreco_nopu/tree_10_1_C6T.root
//////////////////////////////////////////////////////////

#ifndef ecaltree_h
#define ecaltree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ecaltree {
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
   Int_t           nEBRecHits;
   Float_t         etaEBRecHits[10000];   //[nEBRecHits]
   Float_t         phiEBRecHits[10000];   //[nEBRecHits]
   Float_t         ixEBRecHits[10000];   //[nEBRecHits]
   Float_t         iyEBRecHits[10000];   //[nEBRecHits]
   Float_t         energyEBRecHits[10000];   //[nEBRecHits]
   Float_t         timeEBRecHits[10000];   //[nEBRecHits]
   Float_t         swissXEBRecHits[10000];   //[nEBRecHits]
   Float_t         r9EBRecHits[10000];   //[nEBRecHits]
   Int_t           nEERecHits;
   Float_t         etaEERecHits[10000];   //[nEERecHits]
   Float_t         phiEERecHits[10000];   //[nEERecHits]
   Float_t         ixEERecHits[10000];   //[nEERecHits]
   Float_t         iyEERecHits[10000];   //[nEERecHits]
   Float_t         energyEERecHits[10000];   //[nEERecHits]
   Float_t         timeEERecHits[10000];   //[nEERecHits]
   Float_t         swissXEERecHits[10000];   //[nEERecHits]
   Float_t         r9EERecHits[10000];   //[nEERecHits]

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
   TBranch        *b_nEBRecHits;   //!
   TBranch        *b_etaEBRecHits;   //!
   TBranch        *b_phiEBRecHits;   //!
   TBranch        *b_ixEBRecHits;   //!
   TBranch        *b_iyEBRecHits;   //!
   TBranch        *b_energyEBRecHits;   //!
   TBranch        *b_timeEBRecHits;   //!
   TBranch        *b_swissXEBRecHits;   //!
   TBranch        *b_r9EBRecHits;   //!
   TBranch        *b_nEERecHits;   //!
   TBranch        *b_etaEERecHits;   //!
   TBranch        *b_phiEERecHits;   //!
   TBranch        *b_ixEERecHits;   //!
   TBranch        *b_iyEERecHits;   //!
   TBranch        *b_energyEERecHits;   //!
   TBranch        *b_timeEERecHits;   //!
   TBranch        *b_swissXEERecHits;   //!
   TBranch        *b_r9EERecHits;   //!

   ecaltree(TTree *tree=0);
   virtual ~ecaltree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char *outputfile);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ecaltree_cxx
ecaltree::ecaltree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/emanuele/ECAL_612_V1/stdreco_nopu/tree_10_1_C6T.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/emanuele/ECAL_612_V1/stdreco_nopu/tree_10_1_C6T.root");
      }
      f->GetObject("ntp0",tree);

   }
   Init(tree);
}

ecaltree::~ecaltree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ecaltree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ecaltree::LoadTree(Long64_t entry)
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

void ecaltree::Init(TTree *tree)
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
   fChain->SetBranchAddress("nEBRecHits", &nEBRecHits, &b_nEBRecHits);
   fChain->SetBranchAddress("etaEBRecHits", etaEBRecHits, &b_etaEBRecHits);
   fChain->SetBranchAddress("phiEBRecHits", phiEBRecHits, &b_phiEBRecHits);
   fChain->SetBranchAddress("ixEBRecHits", ixEBRecHits, &b_ixEBRecHits);
   fChain->SetBranchAddress("iyEBRecHits", iyEBRecHits, &b_iyEBRecHits);
   fChain->SetBranchAddress("energyEBRecHits", energyEBRecHits, &b_energyEBRecHits);
   fChain->SetBranchAddress("timeEBRecHits", timeEBRecHits, &b_timeEBRecHits);
   fChain->SetBranchAddress("swissXEBRecHits", swissXEBRecHits, &b_swissXEBRecHits);
   fChain->SetBranchAddress("r9EBRecHits", r9EBRecHits, &b_r9EBRecHits);
   fChain->SetBranchAddress("nEERecHits", &nEERecHits, &b_nEERecHits);
   fChain->SetBranchAddress("etaEERecHits", etaEERecHits, &b_etaEERecHits);
   fChain->SetBranchAddress("phiEERecHits", phiEERecHits, &b_phiEERecHits);
   fChain->SetBranchAddress("ixEERecHits", ixEERecHits, &b_ixEERecHits);
   fChain->SetBranchAddress("iyEERecHits", iyEERecHits, &b_iyEERecHits);
   fChain->SetBranchAddress("energyEERecHits", energyEERecHits, &b_energyEERecHits);
   fChain->SetBranchAddress("timeEERecHits", timeEERecHits, &b_timeEERecHits);
   fChain->SetBranchAddress("swissXEERecHits", swissXEERecHits, &b_swissXEERecHits);
   fChain->SetBranchAddress("r9EERecHits", r9EERecHits, &b_r9EERecHits);
   Notify();
}

Bool_t ecaltree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ecaltree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ecaltree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ecaltree_cxx
