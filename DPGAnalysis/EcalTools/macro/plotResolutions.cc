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
#include <TString.h>
#include <TBinomialEfficiencyFitter.h>

#include "RooHZZStyle.C"

double effectiveSigma(TH1 *histo);

void plotResolutions(const char *file) {

  vector<TH1F*> resolutions_EB_1, resolutions_EE_1, resolutions_EB_2, resolutions_EE_2;
  //               [type][eta][pt]
  TH1F *resolutions[3][2][5];
  TFile *tfile = TFile::Open(file);

  // get the histograms from the files
  for(int clustertype=0;clustertype<3;++clustertype) {
    for(int e=0;e<2;++e) {
      std::string suffix = (e==0) ? "EB" : "EE";
      for(int p=0;p<5;++p) {
        char namer[50];
        sprintf(namer,"cluster%d_res_%s_ptbin%d",clustertype,suffix.c_str(),p);
        TH1F* ires = (TH1F*)tfile->Get(namer);
        resolutions[clustertype][e][p] = ires;
        resolutions[clustertype][e][p]->Rebin(10);
      }
    }
  }

  // plot
  gStyle->SetOptStat(0);

  // the sample is pT: [1-100] GeV
  int ptbins[7] = {1,10,20,30,50,100,300};

  TCanvas *c1 = new TCanvas("c1","",600,600);
  TLatex* CP = new TLatex(0.15,0.96, "CMS Gun Simulation                                               #sqrt{s} = 8 TeV");
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.030);
  
  for(int idet=0;idet<2;++idet) {
    for(int p=0;p<5;++p) {
      
      resolutions[0][idet][p]->SetLineColor(kBlack);
      resolutions[1][idet][p]->SetLineColor(kRed+1);
      resolutions[2][idet][p]->SetLineColor(kGreen+1);
      
      float maxy = std::max(resolutions[0][idet][p]->GetMaximum(),std::max(resolutions[1][idet][p]->GetMaximum(),resolutions[2][idet][p]->GetMaximum()));
      maxy*=1.2;

      resolutions[0][idet][p]->GetYaxis()->SetRangeUser(0,maxy);

      resolutions[0][idet][p]->Draw();
      resolutions[1][idet][p]->Draw("sames");
      resolutions[2][idet][p]->Draw("sames");

      char ptrange[50];
      sprintf(ptrange,"%d < p_{T} < %d GeV", ptbins[p], ptbins[p+1]);
      TLatex *ptr = new TLatex(0.15,0.70,ptrange);
      ptr->SetNDC(kTRUE);
      ptr->SetTextSize(0.040);
      ptr->Draw();

      TPaveText *results = new TPaveText(.15,.3,.60,.5,"NDC");
      results->SetBorderSize(0);
      results->SetFillColor (0);
      results->SetTextAlign(12);
      results->SetTextFont(42);
      //results->AddText(Form("weights: #sigma_{eff}=%.1f%%",100.0*effectiveSigma(resolutions[0][idet][p])));
      //results->AddText(Form("5^{th} sample: #sigma_{eff}=%.1f%%",100.0*effectiveSigma(resolutions[1][idet][p])));
      //results->AddText(Form("5^{th} sample NoPU: #sigma_{eff}=%.1f%%",100.0*effectiveSigma(resolutions[2][idet][p])));
      results->AddText(Form("weights: r.m.s.=%.2f%%",100.0*resolutions[0][idet][p]->GetRMS()));
      results->AddText(Form("5^{th} sample: r.m.s.=%.2f%%",100.0*resolutions[1][idet][p]->GetRMS()));
      results->AddText(Form("5^{th} sample NoPU: r.m.s.=%.2f%%",100.0*resolutions[2][idet][p]->GetRMS()));

      results->Draw();

      c1->Update();

      TLegend *legend = new TLegend(0.15,0.75,0.40,0.85,NULL,"brNDC");
      legend->SetBorderSize(     0);
      legend->SetFillColor (     0);
      legend->SetTextAlign (    12);
      legend->SetTextFont  (    42);
  
      legend->AddEntry(resolutions[0][idet][p], "weights");
      legend->AddEntry(resolutions[1][idet][p], "5^{th} sample");
      legend->AddEntry(resolutions[2][idet][p], "5^{th} sample NoPU");

      legend->Draw();
      CP->Draw();

      c1->SaveAs(TString("figures/fifthsampleNoPU_")+resolutions[0][idet][p]->GetName()+TString(".pdf"));
      c1->SaveAs(TString("figures/fifthsampleNoPU_")+resolutions[0][idet][p]->GetName()+TString(".png"));

    }
  }
  
}

double effectiveSigma(TH1 *histo) {
  
  TAxis *xaxis = histo->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  // Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = histo->GetMean();
  Double_t rms = histo->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=histo->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=histo->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=histo->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=histo->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;

}
