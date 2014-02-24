#include <iostream>

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

#include "RooDataSet.h"

using namespace std;

void makeplotstoysTime(TString dir = "./", TString file = "toyresults.root", 
		   TString suffix = "") {

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111111);

  TFile *resfile = TFile::Open(file);
  TTree *ntp = (TTree*)resfile->Get("dataset");

  TH1D* sigTime_pull = new TH1D("sigTime_pull", " ", 20,  -4., 4.  );
  TH1D* sigTime_err  = new TH1D("sigTime_err" , " ", 30,  0., 2.);

  ntp->Project("sigTime_pull","(sigIT_Time_mean_0 - sigIT_Time_meangen)/sigIT_Time_meanerr_0","covQual_0==3");
  ntp->Project("sigTime_err_0","sigIT_Time_mean_0","covQual_0==3");

  // pull distributions, fitted as Gaussians
  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  sigTime_pull->Fit("gaus");
  sigTime_pull->GetXaxis()->SetTitle("Pull N_{sig IT}");
  TString pullstring(dir+"/sigTime_pull");
  pullstring.Append(suffix);
  pullstring.Append(".pdf");
  c1->SaveAs(pullstring);

  // error distributions
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1111111);

}

void makeplotstoys(TString dir = "./", TString file = "toyresults.root", 
		   TString suffix = "") {

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111111);

  TFile *resfile = TFile::Open(file);
  TTree *ntp = (TTree*)resfile->Get("dataset");

  TH1D* N_sigIT_pull = new TH1D("N_sigIT_pull", " ", 20,  -4., 4.  );
  TH1D* N_sigIT_err  = new TH1D("N_sigIT_err" , " ", 200,  0., 200.);

  TH1D* N_bkgOoT_minus1bx_pull = new TH1D("N_bkgOoT_minus1bx_pull", " ", 20, -4., 4.  );
  TH1D* N_bkgOoT_minus1bx_err  = new TH1D("N_bkgOoT_minus1bx_err",  " ", 200, 0., 20. );

  TH1D* N_bkgOoT_minus2bx_pull = new TH1D("N_bkgOoT_minus2bx_pull", " ", 20, -4., 4.  );
  TH1D* N_bkgOoT_minus2bx_err  = new TH1D("N_bkgOoT_minus2bx_err",  " ", 200, 0., 20. );

  TH1D* N_bkgOoT_minus3bx_pull = new TH1D("N_bkgOoT_minus3bx_pull", " ", 20, -4., 4.  );
  TH1D* N_bkgOoT_minus3bx_err  = new TH1D("N_bkgOoT_minus3bx_err",  " ", 200, 0., 20. );

  TH1D* N_bkgOoT_minus4bx_pull = new TH1D("N_bkgOoT_minus4bx_pull", " ", 20, -4., 4.  );
  TH1D* N_bkgOoT_minus4bx_err  = new TH1D("N_bkgOoT_minus4bx_err",  " ", 200, 0., 20. );


  ntp->Project("N_sigIT_pull","(N_sigIT_0 - N_sigITgen)/N_sigITerr_0","covQual_0==3");
  ntp->Project("N_sigIT_err_0","N_sigITerr_0","covQual_0==3");

  ntp->Project("N_bkgOoT_minus1bx_pull","(N_bkgOoT_minus1bx_0 - N_bkgOoT_minus1bxgen)/N_bkgOoT_minus1bxerr_0","covQual_0==3");
  ntp->Project("N_bkgOoT_minus1bx_err","N_bkgOoT_minus1bxerr_0","covQual_0==3");

  ntp->Project("N_bkgOoT_minus2bx_pull","(N_bkgOoT_minus2bx_0 - N_bkgOoT_minus2bxgen)/N_bkgOoT_minus2bxerr_0","covQual_0==3");
  ntp->Project("N_bkgOoT_minus2bx_err","N_bkgOoT_minus2bxerr_0","covQual_0==3");

  ntp->Project("N_bkgOoT_minus3bx_pull","(N_bkgOoT_minus3bx_0 - N_bkgOoT_minus3bxgen)/N_bkgOoT_minus3bxerr_0","covQual_0==3");
  ntp->Project("N_bkgOoT_minus3bx_err","N_bkgOoT_minus3bxerr_0","covQual_0==3");

  ntp->Project("N_bkgOoT_minus4bx_pull","(N_bkgOoT_minus4bx_0 - N_bkgOoT_minus4bxgen)/N_bkgOoT_minus4bxerr_0","covQual_0==3");
  ntp->Project("N_bkgOoT_minus4bx_err","N_bkgOoT_minus4bxerr_0","covQual_0==3");

  // pull distributions, fitted as Gaussians
  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  N_sigIT_pull->Fit("gaus");
  N_sigIT_pull->GetXaxis()->SetTitle("Pull N_{sig IT}");
  TString pullstring(dir+"/N_sigIT_pull_");
  pullstring.Append(suffix);
  pullstring.Append(".pdf");
  c1->SaveAs(pullstring);

  N_bkgOoT_minus1bx_pull->Fit("gaus");
  N_bkgOoT_minus1bx_pull->GetXaxis()->SetTitle("Pull N_{bkg -1 bx}");
  pullstring = TString(dir+"/N_bkgOoT_minus1bx_pull_");
  pullstring.Append(suffix);
  pullstring.Append(".pdf");
  c1->SaveAs(pullstring);

  N_bkgOoT_minus2bx_pull->Fit("gaus");
  N_bkgOoT_minus2bx_pull->GetXaxis()->SetTitle("Pull N_{bkg -2 bx}");
  pullstring = TString(dir+"/N_bkgOoT_minus2bx_pull_");
  pullstring.Append(suffix);
  pullstring.Append(".pdf");
  c1->SaveAs(pullstring);

  N_bkgOoT_minus3bx_pull->Fit("gaus");
  N_bkgOoT_minus3bx_pull->GetXaxis()->SetTitle("Pull N_{bkg -3 bx}");
  pullstring = TString(dir+"/N_bkgOoT_minus3bx_pull_");
  pullstring.Append(suffix);
  pullstring.Append(".pdf");
  c1->SaveAs(pullstring);

  N_bkgOoT_minus4bx_pull->Fit("gaus");
  N_bkgOoT_minus4bx_pull->GetXaxis()->SetTitle("Pull N_{bkg -4 bx}");
  pullstring = TString(dir+"/N_bkgOoT_minus4bx_pull_");
  pullstring.Append(suffix);
  pullstring.Append(".pdf");
  c1->SaveAs(pullstring);

  // error distributions
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1111111);

  N_sigIT_err->GetXaxis()->SetTitle("#sigma(N_{sig IT})");
  N_sigIT_err->Draw();
  TString errstring(dir+"/N_sigIT_err_");
  errstring.Append(suffix);
  errstring.Append(".pdf");
  c2->SaveAs(errstring);

  // N_bkg_err->GetXaxis()->SetTitle("#sigma(N_{bkg})");
  // N_bkg_err->Draw();
  // errstring = TString(dir+"/N_bkg_err_");
  // errstring.Append(suffix);
  // errstring.Append(".eps");
  // c2->SaveAs(errstring);

}

