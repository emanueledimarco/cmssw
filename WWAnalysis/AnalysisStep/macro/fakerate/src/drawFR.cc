#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPaveText.h>

using namespace std;


void prep(TH1* hist, TF1* f, TPaveText* t, TCanvas* canvas, int pad, const char* leptoninfo, const char* isoinfo) {
    hist->SetTitle("");
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.035);
    hist->GetXaxis()->SetTitleSize(0.035);
    hist->GetYaxis()->SetTitleSize(0.035);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    hist->GetYaxis()->SetTitle("Fake Rate (Reco#rightarrow ID+Isolation)");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(0.7);
    hist->SetLineColor(4);
    hist->SetLineWidth(1);
    hist->SetStats(0);
    hist->GetYaxis()->SetRangeUser(0.0,0.4);

    f->SetLineWidth(1);
    f->SetLineStyle(2);
    f->SetLineColor(kBlue);

    t->AddText("#sqrt{s} = 7 TeV 2011, L = 4.9 fb^{-1}");
    t->AddText(leptoninfo);
    t->AddText(isoinfo);
    t->SetBorderSize(0);
    t->SetFillStyle(0);
    t->SetTextAlign(12);
    t->SetTextFont(42);
    t->SetTextSize(0.06);

    canvas->cd(pad);
    for (int i = 0; i < f->GetNumberFreeParameters(); ++i) f->SetParameter(i,0.1);
    hist->Fit(f, "", "", 5.0, 80.0);
    hist->Draw();
    t->Draw();

    for (int i = 0; i < f->GetNumberFreeParameters(); ++i) std::cerr << "\tfloat par" << i << " = " << f->GetParameter(i) << ";" << std::endl;
    std::cerr << std::endl;

}

void drawOneComparison(vector<TH1F*> histos, vector<TString> descr, TString xaxislabel, const char *filename) {

  if(histos.size()>3) {
    cout << "more than 3 histos not implemented." << endl;
    return;
  }

  if(histos.size()!=descr.size()) {
    cout << "description not complete!" << endl;
    return;
  }

  vector<int> colors;
  colors.push_back(kRed+1);
  colors.push_back(kAzure-6);
  colors.push_back(kTeal+3);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetGridx();
  c1->SetGridy();  
  TLegend* legend = new TLegend(0.20, 0.70, 0.43, 0.86);
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);

  for(int i=0;i<(int)histos.size();++i) {
    
    histos[i]->SetMinimum(0);
    histos[i]->SetMaximum(0.4);
    histos[i]->SetMarkerSize(2);
    histos[i]->SetMarkerStyle(20);
    histos[i]->SetMarkerColor(colors[i]);
    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetTitle("");
    if(TString(histos[i]->GetName()).Contains("PU")) {
      histos[i]->Fit("pol1","","same",4,35);
      histos[i]->GetFunction("pol1")->SetLineColor(colors[i]);
    }
    histos[i]->GetXaxis()->SetTitle(xaxislabel);
    histos[i]->GetYaxis()->SetTitleOffset(1.5);
    histos[i]->GetYaxis()->SetTitle("fake rate");

    legend->AddEntry(histos[i],descr[i]);

    if(i==0) { 
      histos[i]->Draw("pe1");
      for(int bin=1;bin<=histos[i]->GetNbinsX();bin++) {
	//cout << "bin i = " << bin << ": " <<  histos[i]->GetBinContent(bin) << " +/- " << histos[i]->GetBinError(bin) << endl;
	cout << "m35_fakeRate[" << bin << "] = " << histos[i]->GetBinContent(bin) << ";" << endl;
      }
      for(int bin=1;bin<=histos[i]->GetNbinsX();bin++) {
	cout << "m35_fakeRate_err[" << bin << "] = " << histos[i]->GetBinError(bin) << ";" << endl;
      }
    }
    else histos[i]->Draw("same pe1");
  }
  legend->Draw();

  TString basename(filename);
  basename.ReplaceAll("Eff","FR");
  c1->SaveAs(basename+TString(".png"));
  c1->SaveAs(basename+TString(".pdf"));

}

void drawOne(vector<TH1F*> set1, const char* desc1, const char* xaxislabel) {
  for(int i=0;i<(int)set1.size();++i) {
    vector<TH1F*> histos;
    vector<TString> desc;
    histos.push_back(set1[i]);
    desc.push_back(TString(desc1));
    if(set1[i]==0) {
      cout << "histogram not found!" << endl;
      continue;
    }
    drawOneComparison(histos,desc,TString(xaxislabel),set1[i]->GetName());
  }

}

void drawOneToOne(vector<TH1F*> set1, vector<TH1F*> set2, const char* desc1, const char* desc2, const char* xaxislabel) {
  if(set1.size()!=set2.size()) {
    cout << "first set and second set of histos have different sizes! ERROR! " << endl;
    return;
  }
  for(int i=0;i<(int)set1.size();++i) {
    vector<TH1F*> histos;
    vector<TString> desc;
    histos.push_back(set1[i]);
    histos.push_back(set2[i]);
    desc.push_back(TString(desc1));
    desc.push_back(TString(desc2));
    if(set1[i]==0 || set2[i]==0) {
      cout << "histogram not found!" << endl;
      continue;
    }
    drawOneComparison(histos,desc,TString(xaxislabel),set1[i]->GetName());
  }

}

void drawOneToTwo(vector<TH1F*> set1, vector<TH1F*> set2, vector<TH1F*> set3, const char* desc1, const char* desc2, const char* desc3, const char* xaxislabel) {
  if(set1.size()!=set2.size() || set1.size()!=set3.size()) {
    cout << "first set and second/third set of histos have different sizes! ERROR! " << endl;
    return;
  }
  for(int i=0;i<(int)set1.size();++i) {
    vector<TH1F*> histos;
    vector<TString> desc;
    histos.push_back(set1[i]);
    histos.push_back(set2[i]);
    histos.push_back(set3[i]);
    desc.push_back(TString(desc1));
    desc.push_back(TString(desc2));
    desc.push_back(TString(desc3));
    if(set1[i]==0 || set2[i]==0 || set3[i]==0) {
      cout << "histogram not found!" << endl;
      continue;
    }
    drawOneComparison(histos,desc,TString(xaxislabel),set1[i]->GetName());
  }

}


// ------------------------------------------
// just draw FR for ONE scenario (reference)
// ------------------------------------------
void drawElectronIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zll1e.root");
  // pt
  TH1F *helhzzPfIsoPtBarrel = (TH1F*)file->Get("helhzzPfIsoPtBarrel_Eff");
  TH1F *helhzzPfIsoPtEndcap = (TH1F*)file->Get("helhzzPfIsoPtEndcap_Eff");
  vector<TH1F*> ptSet1,ptSet2;
  ptSet1.push_back(helhzzPfIsoPtBarrel);
  ptSet1.push_back(helhzzPfIsoPtEndcap);
  drawOne(ptSet1,"Electrons","p_{T} [GeV]");
}

void drawMuonIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zll1m.root");
  // pt
  TH1F *hmuhzzPfIsoPtBarrel = (TH1F*)file->Get("hmuhzzPfIsoPtBarrel_Eff");
  TH1F *hmuhzzPfIsoPtEndcap = (TH1F*)file->Get("hmuhzzPfIsoPtEndcap_Eff");
  vector<TH1F*> ptSet1,ptSet2;
  ptSet1.push_back(hmuhzzPfIsoPtBarrel);
  ptSet1.push_back(hmuhzzPfIsoPtEndcap);
  drawOne(ptSet1,"Muons","p_{T} [GeV]");
}








// compare two scenarios

void compareElectronIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zll1e.root");

  // eta
  // std iso
  TH1F *helhzzPfIsoEtaHighPt = (TH1F*)file->Get("helhzzPfIsoEtaHighPt_Eff");
  TH1F *helhzzPfIsoEtaLowPt = (TH1F*)file->Get("helhzzPfIsoEtaLowPt_Eff");
  // MVA iso
  TH1F *helhzzMvaLooseEtaHighPt = (TH1F*)file->Get("helhzzMvaLooseEtaHighPt_Eff");
  TH1F *helhzzMvaLooseEtaLowPt = (TH1F*)file->Get("helhzzMvaLooseEtaLowPt_Eff");
  TH1F *helhzzMvaTightEtaHighPt = (TH1F*)file->Get("helhzzMvaTightEtaHighPt_Eff");
  TH1F *helhzzMvaTightEtaLowPt = (TH1F*)file->Get("helhzzMvaTightEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2, etaSet3;
  etaSet1.push_back(helhzzPfIsoEtaHighPt);
  etaSet1.push_back(helhzzPfIsoEtaLowPt);
  etaSet2.push_back(helhzzMvaLooseEtaHighPt);
  etaSet2.push_back(helhzzMvaLooseEtaLowPt);
  etaSet3.push_back(helhzzMvaTightEtaHighPt);
  etaSet3.push_back(helhzzMvaTightEtaLowPt);
 
  drawOneToTwo(etaSet1,etaSet2,etaSet3,"H#rightarrow ZZ MVA ref.", "H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","#eta");

  // pt
  // std iso
  TH1F *helhzzPfIsoPtBarrel = (TH1F*)file->Get("helhzzPfIsoPtBarrel_Eff");
  TH1F *helhzzPfIsoPtEndcap = (TH1F*)file->Get("helhzzPfIsoPtEndcap_Eff");
  // MVA iso
  TH1F *helhzzMvaLoosePtBarrel = (TH1F*)file->Get("helhzzMvaLoosePtBarrel_Eff");
  TH1F *helhzzMvaLoosePtEndcap = (TH1F*)file->Get("helhzzMvaLoosePtEndcap_Eff");
  TH1F *helhzzMvaTightPtBarrel = (TH1F*)file->Get("helhzzMvaTightPtBarrel_Eff");
  TH1F *helhzzMvaTightPtEndcap = (TH1F*)file->Get("helhzzMvaTightPtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2, ptSet3;
  ptSet1.push_back(helhzzPfIsoPtBarrel);
  ptSet1.push_back(helhzzPfIsoPtEndcap);
  ptSet2.push_back(helhzzMvaLoosePtBarrel);
  ptSet2.push_back(helhzzMvaLoosePtEndcap);
  ptSet3.push_back(helhzzMvaTightPtBarrel);
  ptSet3.push_back(helhzzMvaTightPtEndcap);

  drawOneToTwo(ptSet1,ptSet2,ptSet3,"H#rightarrow ZZ MVA ref","H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","p_{T} [GeV]");

  // PU
  TH1F *helhzzPfIsoPUBarrel = (TH1F*)file->Get("helhzzPfIsoPUBarrel_Eff");
  TH1F *helhzzPfIsoPUEndcap = (TH1F*)file->Get("helhzzPfIsoPUEndcap_Eff");
  // MVA iso
  TH1F *helhzzMvaLoosePUBarrel = (TH1F*)file->Get("helhzzMvaLoosePUBarrel_Eff");
  TH1F *helhzzMvaLoosePUEndcap = (TH1F*)file->Get("helhzzMvaLoosePUEndcap_Eff");
  TH1F *helhzzMvaTightPUBarrel = (TH1F*)file->Get("helhzzMvaTightPUBarrel_Eff");
  TH1F *helhzzMvaTightPUEndcap = (TH1F*)file->Get("helhzzMvaTightPUEndcap_Eff");

  vector<TH1F*> puSet1, puSet2, puSet3;
  puSet1.push_back(helhzzPfIsoPUBarrel);
  puSet1.push_back(helhzzPfIsoPUEndcap);
  puSet2.push_back(helhzzMvaLoosePUBarrel);
  puSet2.push_back(helhzzMvaLoosePUEndcap);
  puSet3.push_back(helhzzMvaTightPUBarrel);
  puSet3.push_back(helhzzMvaTightPUEndcap);

  drawOneToTwo(puSet1,puSet2,puSet3,"H#rightarrow ZZ MVA ref","H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","# vertices");


}


void compareMuonIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zll1m.root");

  // eta
  TH1F *hmuhzzPfIsoEtaHighPt = (TH1F*)file->Get("hmuhzzPfIsoEtaHighPt_Eff");
  TH1F *hmuhzzPfIsoEtaLowPt = (TH1F*)file->Get("hmuhzzPfIsoEtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *hmuhzzMvaPfIsoEtaHighPt = (TH1F*)file->Get("hmuhzzMvaPfIsoEtaHighPt_Eff");
  TH1F *hmuhzzMvaPfIsoEtaLowPt = (TH1F*)file->Get("hmuhzzMvaPfIsoEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(hmuhzzPfIsoEtaHighPt);
  etaSet1.push_back(hmuhzzPfIsoEtaLowPt);
  etaSet2.push_back(hmuhzzMvaPfIsoEtaHighPt);
  etaSet2.push_back(hmuhzzMvaPfIsoEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"H#rightarrow ZZ PF iso","H#rightarrow ZZ MVA PF Iso","#eta");

  // pt
  TH1F *hmuhzzPfIsoPtBarrel = (TH1F*)file->Get("hmuhzzPfIsoPtBarrel_Eff");
  TH1F *hmuhzzPfIsoPtEndcap = (TH1F*)file->Get("hmuhzzPfIsoPtEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *hmuhzzMvaPfIsoPtBarrel = (TH1F*)file->Get("hmuhzzMvaPfIsoPtBarrel_Eff");
  TH1F *hmuhzzMvaPfIsoPtEndcap = (TH1F*)file->Get("hmuhzzMvaPfIsoPtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(hmuhzzPfIsoPtBarrel);
  ptSet1.push_back(hmuhzzPfIsoPtEndcap);
  ptSet2.push_back(hmuhzzMvaPfIsoPtBarrel);
  ptSet2.push_back(hmuhzzMvaPfIsoPtEndcap);

  drawOneToOne(ptSet1,ptSet2,"H#rightarrow ZZ PF iso","H#rightarrow ZZ MVA PF iso","p_{T} [GeV]");

  // PU
  TH1F *hmuhzzPfIsoPUBarrel = (TH1F*)file->Get("hmuhzzPfIsoPUBarrel_Eff");
  TH1F *hmuhzzPfIsoPUEndcap = (TH1F*)file->Get("hmuhzzPfIsoPUEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *hmuhzzMvaPfIsoPUBarrel = (TH1F*)file->Get("hmuhzzMvaPfIsoPUBarrel_Eff");
  TH1F *hmuhzzMvaPfIsoPUEndcap = (TH1F*)file->Get("hmuhzzMvaPfIsoPUEndcap_Eff");

  vector<TH1F*> puSet1, puSet2;
  puSet1.push_back(hmuhzzPfIsoPUBarrel);
  puSet1.push_back(hmuhzzPfIsoPUEndcap);
  puSet2.push_back(hmuhzzMvaPfIsoPUBarrel);
  puSet2.push_back(hmuhzzMvaPfIsoPUEndcap);

  drawOneToOne(puSet1,puSet2,"H#rightarrow ZZ PF iso","H#rightarrow ZZ MVA PF iso","# vertices");


}

enum {
  kMuPfIso = 0,
  kMuMvaPfIso
};

void fitMuonIdsUnbiased(int wp) {

  if(wp>1) {
    cout << "wp must be kPfIso / kMvaPfIso" << endl;
    return;
  }
  TString wpstr = (wp==kMuPfIso) ? TString("PfIso") : TString("MvaPfIso");

  TFile *file = TFile::Open("fakerates_zll1m.root");

  TH1F *barrel = (TH1F*)file->Get(TString("hmuhzz")+wpstr+TString("PtBarrel_Eff"));
  TH1F *endcap = (TH1F*)file->Get(TString("hmuhzz")+wpstr+TString("PtEndcap_Eff"));

  std::vector<TH1F*> fakerates1D;
  fakerates1D.push_back(barrel);
  fakerates1D.push_back(endcap);

  TCanvas* canvas = new TCanvas("fakerate","Fake Rates",1000,1000);
  canvas->Divide(1,2);

  TPaveText* text1  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
  TPaveText* text2  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");

  TF1 * f1 = new TF1("f1","[0]+[1]*(exp(-[2]*x) - exp(-[3]*x))");
  for (int i = 0; i <=f1->GetNumberFreeParameters(); ++i) f1->SetParameter(i,1);
  f1->SetParLimits(1,0,5);

  prep(barrel, f1, text1, canvas, 1, "Barrel", "");
  prep(endcap, f1, text2, canvas, 2, "Endcap", "");
  
  TString pdffile = (wp==kMuPfIso) ? TString("fakerates_zll1m_PfIso.pdf") : TString("fakerates_zll1m_MvaPfIso.pdf");
  canvas->SaveAs(pdffile);

  TString pngfile = (wp==kMuPfIso) ? TString("fakerates_zll1m_PfIso.png") : TString("fakerates_zll1m_MvaPfIso.png");
  canvas->SaveAs(pngfile);

}

enum {
  kElePfIso = 0,
  kEleMvaLoose,
  kEleMvaTight
};

void fitElectronIdsUnbiased(int wp) {

  if(wp>2) {
    cout << "wp must be kElePfIso / kMvaLoose / kMvaTight" << endl;
    return;
  }
  TString wpstr;
  if(wp==kElePfIso) wpstr = TString("PfIso");
  else if(wp==kEleMvaLoose) wpstr = TString("MvaLoose");
  else  wpstr = TString("MvaTight");

  TFile *file = TFile::Open("fakerates_zll1e.root");

  TH1F *barrel = (TH1F*)file->Get(TString("helhzz")+wpstr+TString("PtBarrel_Eff"));
  TH1F *endcap = (TH1F*)file->Get(TString("helhzz")+wpstr+TString("PtEndcap_Eff"));

  std::vector<TH1F*> fakerates1D;
  fakerates1D.push_back(barrel);
  fakerates1D.push_back(endcap);

  TCanvas* canvas = new TCanvas("fakerate","Fake Rates",1000,1000);
  canvas->Divide(1,2);

  TPaveText* text1  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
  TPaveText* text2  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");

  TF1 * f1 = new TF1("f1","[0]+[1]*(exp(-[2]*x) - exp(-[3]*x))");
  for (int i = 0; i <=f1->GetNumberFreeParameters(); ++i) f1->SetParameter(i,1);

  prep(barrel, f1, text1, canvas, 1, "Barrel", "");
  prep(endcap, f1, text2, canvas, 2, "Endcap", "");
  
  TString pdffile;
  if(wp==kElePfIso) pdffile=TString("fakerates_zll1e_PfIso.pdf");
  else if(wp==kEleMvaLoose) pdffile=TString("fakerates_zll1e_MvaLoosePfIso.pdf");
  else pdffile=TString("fakerates_zll1e_MvaTightPfIso.pdf");

  canvas->SaveAs(pdffile);

  TString pngfile;
  if(wp==kElePfIso) pngfile=TString("fakerates_zll1e_PfIso.png");
  else if(wp==kEleMvaLoose) pngfile=TString("fakerates_zll1e_MvaLoosePfIso.png");
  else pngfile=TString("fakerates_zll1e_MvaTightPfIso.png");

  canvas->SaveAs(pngfile);

}
