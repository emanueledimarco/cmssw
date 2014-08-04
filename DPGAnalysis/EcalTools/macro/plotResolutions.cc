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

void plotResolutions(const char *file) {

  vector<TH1F*> resolutions_EB_1, resolutions_EE_1, resolutions_EB_2, resolutions_EE_2;
  //               [type][eta][pt]
  TH1F *resolutions[3][2][6];
  TFile *tfile = TFile::Open(file);

  // get the histograms from the files
  for(int clustertype=0;clustertype<3;++clustertype) {
    for(int e=0;e<2;++e) {
      std::string suffix = (e==0) ? "EB" : "EE";
      for(int p=0;p<6;++p) {
        char namer[50];
        sprintf(namer,"cluster%d_res_%s_ptbin%d",clustertype,suffix.c_str(),p);
        TH1F* ires = (TH1F*)tfile->Get(namer);
        resolutions[clustertype][e][p] = ires;
        resolutions[clustertype][e][p]->Rebin(10);
      }
    }
  }

  // plot
  gStyle->SetOptStat(1100);

  int ptbins[7] = {1,10,20,30,50,100,300};

  TCanvas *c1 = new TCanvas("c1","",600,600);
  TLatex* CP = new TLatex(0.15,0.96, "CMS Gun Simulation                                               #sqrt{s} = 8 TeV");
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.030);
  
  for(int idet=0;idet<2;++idet) {
    for(int p=0;p<6;++p) {
      
      resolutions[0][idet][p]->SetLineColor(kBlack);
      resolutions[1][idet][p]->SetLineColor(kRed+1);
      resolutions[2][idet][p]->SetLineColor(kGreen+1);

      resolutions[0][idet][p]->Draw();
      resolutions[1][idet][p]->Draw("sames");
      resolutions[2][idet][p]->Draw("sames");

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
  
      legend->AddEntry(resolutions[0][idet][p], "weights");
      legend->AddEntry(resolutions[1][idet][p], "fit");
      legend->AddEntry(resolutions[2][idet][p], "fit + NoPU");

      TPaveStats *p1 = (TPaveStats*)resolutions[0][idet][p]->GetListOfFunctions()->FindObject("stats");
      p1->SetTextColor(kBlack);
      p1->SetX1NDC(0.7);
      p1->SetX2NDC(0.9);
      p1->SetY1NDC(0.7);
      p1->SetY2NDC(0.9);
      p1->Draw();

      TPaveStats *p2 = (TPaveStats*)resolutions[1][idet][p]->GetListOfFunctions()->FindObject("stats");
      p2->SetTextColor(kRed+1);
      p2->SetX1NDC(0.7);
      p2->SetX2NDC(0.9);
      p2->SetY1NDC(0.5);
      p2->SetY2NDC(0.7);
      p2->Draw();

      TPaveStats *p3 = (TPaveStats*)resolutions[2][idet][p]->GetListOfFunctions()->FindObject("stats");
      p3->SetTextColor(kGreen+1);
      p3->SetX1NDC(0.7);
      p3->SetX2NDC(0.9);
      p3->SetY1NDC(0.3);
      p3->SetY2NDC(0.5);
      p3->Draw();

      legend->Draw();
      CP->Draw();

      c1->SaveAs(TString("figures/")+resolutions[0][idet][p]->GetName()+TString(".pdf"));
      c1->SaveAs(TString("figures/")+resolutions[0][idet][p]->GetName()+TString(".png"));

    }
  }
  
}
