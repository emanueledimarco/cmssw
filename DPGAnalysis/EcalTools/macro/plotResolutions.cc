// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TBinomialEfficiencyFitter.h>

#include "RooHZZStyle.C"

void plotResolutions(const char *file1, const char* file2) {

  vector<TH1F*> resolutions_EB_1, resolutions_EE_1, resolutions_EB_2, resolutions_EE_2;
  TFile *tfile1 = TFile::Open(file1);
  TFile *tfile2 = TFile::Open(file2);

  // get the histograms from the files
  for(int file=0;file<2;++file) {
    if(file==0) tfile1->cd();
    else tfile2->cd();

    for(int e=0;e<2;++e) {
      std::string suffix = (e==0) ? "EB" : "EE";
      for(int p=0;p<6;++p) {
        char namer[50];
        sprintf(namer,"res_%s_ptbin%d",suffix.c_str(),p);
        TH1F* ires = (TH1F*)gDirectory->Get(namer);
        if(file==0) {
          if(e==0) resolutions_EB_1.push_back(ires);
          else resolutions_EE_1.push_back(ires);
        } else {
          if(e==0) resolutions_EB_2.push_back(ires);
          else resolutions_EE_2.push_back(ires);
        }
      }
    }

  }

  // plot
  gStyle->SetOptStat(1100);

  int ptbins[7] = {1,10,20,30,50,100,300};

  vector< vector<TH1F*> > resolutions_1, resolutions_2;
  resolutions_1.push_back(resolutions_EB_1);
  resolutions_1.push_back(resolutions_EE_1);
  resolutions_2.push_back(resolutions_EB_2);
  resolutions_2.push_back(resolutions_EE_2);

  TCanvas *c1 = new TCanvas("c1","",600,600);
  TLatex* CP = new TLatex(0.15,0.96, "CMS Gun Simulation                                               #sqrt{s} = 8 TeV");
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.030);
  

  for(int idet=0;idet<2;++idet) {
    for(int p=0;p<6;++p) {
      
      (resolutions_1[idet])[p]->SetLineColor(kBlack);
      (resolutions_2[idet])[p]->SetLineColor(kRed+1);
      (resolutions_1[idet])[p]->SetLineWidth(2);
      (resolutions_2[idet])[p]->SetLineWidth(2);
      
      (resolutions_1[idet])[p]->Draw();
      (resolutions_2[idet])[p]->Draw("sames");

      char ptrange[50];
      sprintf(ptrange,"%d < p_{T} < %d GeV", ptbins[p], ptbins[p+1]);
      TLatex *ptr = new TLatex(0.15,0.70,ptrange);
      ptr->SetNDC(kTRUE);
      ptr->SetTextSize(0.040);
      ptr->Draw();

      c1->Update();

      TLegend *legend = new TLegend(0.15,0.75,0.40,0.85,NULL,"brNDC");
      legend->SetBorderSize(     0);
      legend->SetFillColor (     0);
      legend->SetTextAlign (    12);
      legend->SetTextFont  (    42);
  
      legend->AddEntry((resolutions_1[idet])[p], "std. reco");
      legend->AddEntry((resolutions_2[idet])[p], "PU sub.");

      TPaveStats *p1 = (TPaveStats*)(resolutions_1[idet])[p]->GetListOfFunctions()->FindObject("stats");
      p1->SetTextColor(kBlack);
      p1->SetX1NDC(0.7);
      p1->SetX2NDC(0.9);
      p1->SetY1NDC(0.7);
      p1->SetY2NDC(0.9);
      p1->Draw();

      TPaveStats *p2 = (TPaveStats*)(resolutions_2[idet])[p]->GetListOfFunctions()->FindObject("stats");
      p2->SetTextColor(kRed+1);
      p2->SetX1NDC(0.7);
      p2->SetX2NDC(0.9);
      p2->SetY1NDC(0.5);
      p2->SetY2NDC(0.7);
      p2->Draw();

      legend->Draw();
      CP->Draw();

      c1->SaveAs(TString("figures/")+(resolutions_1[idet])[p]->GetName()+TString(".pdf"));
      c1->SaveAs(TString("figures/")+(resolutions_1[idet])[p]->GetName()+TString(".png"));

    }
  }

}
