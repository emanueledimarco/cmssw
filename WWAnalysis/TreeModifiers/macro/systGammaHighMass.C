#include "TFile.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"

void getSyst(int channel) {
  
  TFile *fNom = TFile::Open("HighMass/8TeV/Nominal/parameters_channel2.root");
  TGraphErrors *gSigmaBWNom = (TGraphErrors*)fNom->Get("gSigmaBW");
  TFile *f1Up = TFile::Open("HighMass/8TeV/Sigma1Up/parameters_channel2.root");
  TGraphErrors *gSigmaBW1Up = (TGraphErrors*)f1Up->Get("gSigmaBW");
  TFile *f1Dn = TFile::Open("HighMass/8TeV/Sigma1Down/parameters_channel2.root");
  TGraphErrors *gSigmaBW1Dn = (TGraphErrors*)f1Dn->Get("gSigmaBW");

  gSigmaBWNom->SetMarkerColor(kBlack);
  gSigmaBWNom->SetLineColor(kBlack);

  gSigmaBW1Up->SetMarkerColor(kRed+2);
  gSigmaBW1Up->SetLineColor(kRed+2);

  gSigmaBW1Dn->SetMarkerColor(kGreen+2);
  gSigmaBW1Dn->SetLineColor(kGreen+2);

  gSigmaBWNom->GetYaxis()->SetRangeUser(0,200);

  Int_t n = gSigmaBWNom->GetN();
  Double_t *x = gSigmaBWNom->GetX();
  Double_t *yN = gSigmaBWNom->GetY();
  Double_t *yEN = gSigmaBWNom->GetEY();
  Double_t *yU = gSigmaBW1Up->GetY();
  Double_t *yEU = gSigmaBW1Up->GetEY();
  Double_t *yD = gSigmaBW1Dn->GetY();
  Double_t *yED = gSigmaBW1Dn->GetEY();

  TGraphErrors *gDiff = new TGraphErrors(n);
  gDiff->GetYaxis()->SetRangeUser(0,10);
  gDiff->SetMarkerStyle(8);
  gDiff->SetMarkerSize(1.5);

  for(int i=0; i<n; ++i) {
    float syst=std::max(fabs(yD[i]-yN[i]),fabs(yU[i]-yN[i]))/yN[i];
    if(syst<2) {
      gDiff->SetPoint(i,x[i],syst);
      gDiff->SetPointError(i,0,0); // for the moment
    } else gDiff->SetPoint(i,x[i],0);
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  gDiff->GetXaxis()->SetTitle("m_H [GeV]");
  gDiff->GetYaxis()->SetTitle("#Gamma_{BW} rel. unc.");
  gDiff->Draw("ape");
  c1->SaveAs("GammaBW_channel2.pdf");

}
