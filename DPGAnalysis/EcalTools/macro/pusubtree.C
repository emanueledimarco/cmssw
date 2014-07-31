#define pusubtree_cxx
#include "pusubtree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void pusubtree::Loop(const char *outputfilename)
{
//   In a ROOT session, you can do:
//      Root > .L pusubtree.C
//      Root > pusubtree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile *fileo = TFile::Open(outputfilename,"recreate");

   TH1F *res = new TH1F("res","",100,-0.25,0.25);

   std::vector<TH1F*> resolutions_EB, resolutions_EE;
   float ptbins[7] = {1,10,20,30,50,100,300};
   for(int e=0;e<2;++e) {
     std::string suffix = (e==0) ? "EB" : "EE";
     for(int p=0;p<7;++p) {
       char namer[50];
       sprintf(namer,"res_%s_ptbin%d",suffix.c_str(),p);
       TH1F* ires = (TH1F*)res->Clone(namer);
       if(e==0) resolutions_EB.push_back(ires);
       else resolutions_EE.push_back(ires);
     }
   }

   Long64_t nentries = fChain->GetEntries();

   std::cout << "Total entries = " << nentries << std::endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry % 1000 == 0) std::cout << "Processing entry " << jentry << std::endl;

      if(idMc[0] != 22) continue;
      TVector3 pGamma;
      pGamma.SetMagThetaPhi(pMc[0],thetaMc[0],phiMc[0]);
      
      if(fabs(etaMc[0])<1.479) {
        for(int isc=0; isc<nEBSuperClusters; ++isc) {
          TVector3 scdir;
          scdir.SetMagThetaPhi(energyEBSuperClusters[isc],
                               thetaEBSuperClusters[isc],
                               phiEBSuperClusters[isc]);
          if(mcmatch(scdir,pGamma,0.3)) {
            int ptbin=0;
            for(int ipt=0; ipt<6; ++ipt) {
              if(pGamma.Pt() >= ptbins[ipt] && pGamma.Pt() < ptbins[ipt+1]) {
                ptbin=ipt; break;
              }
            }
            resolutions_EB[ptbin]->Fill((energyEBSuperClusters[isc]-energyMc[isc])/energyMc[isc]);
          }
        }
      } else {
        for(int isc=0; isc<nEESuperClusters; ++isc) {
          TVector3 scdir;
          scdir.SetMagThetaPhi(energyEESuperClusters[isc],
                               thetaEESuperClusters[isc],
                               phiEESuperClusters[isc]);
          if(mcmatch(scdir,pGamma,0.3)) {
            int ptbin=0;
            for(int ipt=0; ipt<6; ++ipt) {
              if(pGamma.Pt() >= ptbins[ipt] && pGamma.Pt() < ptbins[ipt+1]) {
                ptbin=ipt; break;
              }
            }
            resolutions_EE[ptbin]->Fill((energyEESuperClusters[isc]-energyMc[isc])/energyMc[isc]);
          }
        }
      }

   } // loop over the events

   for(int i=0; i<(int) resolutions_EB.size(); ++i) {
     resolutions_EB[i]->Write();
     resolutions_EE[i]->Write();
   }

   fileo->Close();
   
}
