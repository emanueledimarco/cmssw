#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "RooHZZStyle.C"

#include <sstream>
#include <iostream>

using namespace std;

TLatex *CMSPreliminary(float lumi7TeV=5.1, float lumi8TeV=19.6) {

  stringstream line;
  line << "CMS Preliminary  #sqrt{s}=7 TeV, L=" << lumi7TeV << " fb^{-1}  #sqrt{s}=8 TeV, L=" << lumi8TeV << " fb^{-1}";
  TLatex* CP = new TLatex(0.15,0.96, line.str().c_str());
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.032);

  return CP;

}

TPaveText *text(const char *txt, float x1, float y1, float x2, float y2) {
  TPaveText *text = new TPaveText(x1,y1,x2,y2,"brNDC");
  text->AddText(txt);
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.05);
  return text;
}

void analyzeMassToys(int fitdimensions) {

  // fit limits
  float massMin=120;
  float massMax=130;

  stringstream oss;
  if(fitdimensions==1) oss << "histos1D.root";
  else if(fitdimensions==2) oss << "histos2D.root";
  else if(fitdimensions==3) oss << "histos3D.root";
  
  TFile *histos = TFile::Open(oss.str().c_str(),"recreate");
        
  stringstream fss;
  if(fitdimensions==1) fss << "results/c1DNoMassErr/toys.root";
  else if(fitdimensions==2) fss << "results/c1DMassErr/toys.root";
  else if(fitdimensions==3) fss << "results/c2DMassErr/toys.root";
  else {
    cout << "fitdimensions should be 1,2,3" << endl;
    return;
  }

  TFile *tfile = TFile::Open(fss.str().c_str());
  TTree *toys = (TTree*)tfile->Get("limit");
  
  TH1F *massh = new TH1F("massh","",25,123.,129.);
  TH1F *pullh = new TH1F("pullh","",25,-5.,5.);
  TH1F *errh = new TH1F("errh","",25,0.,2.);
  TH1F *errl = new TH1F("errl","",25,0.,2.);

  errh->GetXaxis()->SetTitle("mass error (GeV)");
  errl->GetXaxis()->SetTitle("mass error (GeV)");
  pullh->GetXaxis()->SetTitle("pull");

  errh->GetYaxis()->SetTitle("n. toys");
  errl->GetYaxis()->SetTitle("n. toys");
  pullh->GetYaxis()->SetTitle("n. toys");

  double genmh;
  float MH;
  float quantileExpected;
  toys->SetBranchAddress("mh", &genmh);
  toys->SetBranchAddress("MH", &MH);
  toys->SetBranchAddress("quantileExpected", &quantileExpected);

  float mass, lowere, uppere;
  mass=lowere=uppere=0;
  for(int i=0; i<(int)toys->GetEntries(); ++i) {
    toys->GetEntry(i);
    // this to distinguish -1 sigma / +1 sigma
    if(quantileExpected==1) {
      massh->Fill(MH);
      mass = MH;
    } else {
      if((i-1)%4==0) {
        lowere = mass-MH;
        if(mass!=massMin && mass!=massMax && lowere>0.1) errl->Fill(lowere);
      }
      if((i-2)%4==0) {
        uppere = MH-mass;
        if(mass!=massMin && mass!=massMax && uppere>0.1) errh->Fill(uppere);
      }
      if((i-3)%4==0) {
        float avgerr = (lowere + uppere)/2.;
        if(mass!=massMin && mass!=massMax && lowere>0.1 && uppere>0.1) pullh->Fill((mass-genmh)/avgerr);
      }
    }
  }

  histos->cd();
  massh->Write();
  pullh->Write();
  errh->Write();
  errl->Write();
  histos->Close();

}

void finalizeMassToys() {

  // data fits (uppere)
  //  float dm1h=0.572;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *fit1D = TFile::Open("histos1D.root");
  TH1F *pull1D = (TH1F*)fit1D->Get("pullh");
  TH1F *errh1D = (TH1F*)fit1D->Get("errh");

  TFile *fit2D = TFile::Open("histos2D.root");
  TH1F *pull2D = (TH1F*)fit2D->Get("pullh");
  TH1F *errh2D = (TH1F*)fit2D->Get("errh");

  TFile *fit3D = TFile::Open("histos3D.root");
  TH1F *pull3D = (TH1F*)fit3D->Get("pullh");
  TH1F *errh3D = (TH1F*)fit3D->Get("errh");

  float m1=errh1D->GetMean();
  float m2=errh2D->GetMean();
  float m3=errh3D->GetMean();

  float me1=errh1D->GetMeanError();
  float me2=errh2D->GetMeanError();
  float me3=errh3D->GetMeanError();

  stringstream m1s,m2s,m3s;
  m1s.precision(2); m2s.precision(2); m3s.precision(2);
  m1s << std::fixed << "m_{4l}. Mean = " << m1 << " #pm " << me1;
  m2s << std::fixed << "m_{4l}, #sigma_{m}. Mean = " << m2 << " #pm " << me2;
  m3s << std::fixed << "m_{4l}, #sigma_{m}, KD. Mean = " << m3 << " #pm " << me3;

  std::vector<TH1F*> pulls,errs;
  pulls.push_back(pull1D);
  pulls.push_back(pull2D);
  pulls.push_back(pull3D);
  errs.push_back(errh1D);
  errs.push_back(errh2D);
  errs.push_back(errh3D);

  
  TPaveText *text = new TPaveText(0.50,0.50,0.85,0.60,"brNDC");
  text->AddText(m1s.str().c_str());
  text->AddText(m2s.str().c_str());
  text->AddText(m3s.str().c_str());
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.03);


  // draw the legend
  TLegend *legend = new TLegend(0.50,0.65,0.85,0.85,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);
  legend->AddEntry(errh1D, "m_{4l}","l");
  legend->AddEntry(errh2D, "m_{4l},#sigma_{m}","l");
  legend->AddEntry(errh3D, "m_{4l},#sigma_{m},KD","l");

  TCanvas *c1 = new TCanvas("c1","",750,750);
  errs[0]->SetMaximum(3000);
  for(int i=0;i<3;++i) {
    errs[i]->SetLineColor(i+1);
    if(i==0) errs[i]->Draw();
    else {
      errs[i]->Scale(errs[0]->Integral()/errs[i]->Integral());
      errs[i]->Draw("same");
    }
  }
  legend->Draw();
  text->Draw();
  c1->SaveAs("errors.png");

  for(int i=0;i<3;++i) {
    pulls[i]->Draw();
    pulls[i]->Fit("gaus");
    stringstream pullc;
    pullc << "pull_" << i+1 << "D.png";
    c1->SaveAs(pullc.str().c_str());
  }  
 
}



void plotScanByChannel(int ndim) {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  
  TStyle *mystyle = RooHZZStyle("ZZ");
  mystyle->cd();

  stringstream file4mu,file4e,file2e2mu,filecomb;
  file4mu   << "higgsCombineSCAN" << ndim << "D4muFast.MultiDimFit.mH125.8.root";
  file4e    << "higgsCombineSCAN" << ndim << "D4eFast.MultiDimFit.mH125.8.root";
  file2e2mu << "higgsCombineSCAN" << ndim << "D2e2muFast.MultiDimFit.mH125.8.root";
  filecomb  << "higgsCombineSCAN" << ndim << "DcombFast.MultiDimFit.mH125.8.root";

  TFile *fit4mu = TFile::Open(file4mu.str().c_str());
  TTree *tree4mu = (TTree*)fit4mu->Get("limit");

  TFile *fit4e = TFile::Open(file4e.str().c_str());
  TTree *tree4e = (TTree*)fit4e->Get("limit");

  TFile *fit2e2mu = TFile::Open(file2e2mu.str().c_str());
  TTree *tree2e2mu = (TTree*)fit2e2mu->Get("limit");

  TFile *fitcomb = TFile::Open(filecomb.str().c_str());
  TTree *treecomb = (TTree*)fitcomb->Get("limit");

  TGraph *g4mu = new TGraph();
  TGraph *g4e = new TGraph();
  TGraph *g2e2mu = new TGraph();
  TGraph *gcomb = new TGraph();

  vector<TTree*> trees;
  trees.push_back(tree4mu);
  trees.push_back(tree4e);
  trees.push_back(tree2e2mu);
  trees.push_back(treecomb);

  vector<TGraph*> graphs;
  graphs.push_back(g4mu);
  graphs.push_back(g4e);
  graphs.push_back(g2e2mu);
  graphs.push_back(gcomb);


  for(int cha=0; cha<(int)trees.size(); ++cha) {
    
    cout << "Analyzing scan for channel = " << cha << endl;

    float MH;
    float deltaNLL;
    trees[cha]->SetBranchAddress("MH", &MH);
    trees[cha]->SetBranchAddress("deltaNLL", &deltaNLL);

    if(cha<3) graphs[cha]->SetLineWidth(2);
    else graphs[cha]->SetLineWidth(4);

    for(int i=1; i<(int)trees[cha]->GetEntries();++i) {
       trees[cha]->GetEntry(i);
       graphs[cha]->SetPoint(i-1,MH,2*deltaNLL);
    }

  }

  graphs[0]->SetLineColor(kRed);
  graphs[1]->SetLineColor(kGreen+2);
  graphs[2]->SetLineColor(kBlue);
  graphs[3]->SetLineColor(kBlack);

  TCanvas *c1 = new TCanvas("c1","",750,750);
  graphs[0]->GetXaxis()->SetRangeUser(122,130);
  graphs[0]->GetYaxis()->SetRangeUser(0,10);
  graphs[0]->GetXaxis()->SetTitle("m_{X} (GeV)");
  graphs[0]->GetYaxis()->SetTitle("-2 #Delta lnL");

  graphs[0]->Draw("al");
  for(int i=1;i<4;++i) graphs[i]->Draw("l");

  // draw the legend
  TLegend *legend = new TLegend(0.55,0.65,0.85,0.90,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);
  legend->AddEntry(graphs[3], "Combined","l");
  legend->AddEntry(graphs[1], "H #rightarrow ZZ #rightarrow 4e ","l");
  legend->AddEntry(graphs[0], "H #rightarrow ZZ #rightarrow 4#mu ","l");
  legend->AddEntry(graphs[2], "H #rightarrow ZZ #rightarrow 2e2#mu ","l");
  legend->Draw();

  TLine *line1 = new TLine(122,1,130,1);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(3.0);
  line1->Draw();

  TLine *line2 = new TLine(122,4,130,4);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(1.5);
  line2->Draw();

  TLatex *CP = CMSPreliminary();
  CP->Draw();

  TPaveText *comment;
  if(ndim==1) comment = text("H #rightarrow ZZ, 1D",0.20,0.90,0.40,0.90);
  if(ndim==2) comment = text("H #rightarrow ZZ, 2D",0.20,0.90,0.40,0.90);
  if(ndim==3) comment = text("H #rightarrow ZZ, 3D",0.20,0.90,0.40,0.90);
  comment->Draw();


  stringstream outnamepdf;
  outnamepdf << "scanmass_bychannel_" << ndim << "D.pdf";
  c1->SaveAs(outnamepdf.str().c_str());
  stringstream outnamepng;
  outnamepng << "scanmass_bychannel_" << ndim << "D.png";
  c1->SaveAs(outnamepng.str().c_str());
}

void plotScanByDim(int channel) {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  
  TStyle *mystyle = RooHZZStyle("ZZ");
  mystyle->cd();

  string channelstr;
  if (channel==0) channelstr = "4mu";
  if (channel==1) channelstr = "4e";
  if (channel==2) channelstr = "2e2mu";

  stringstream file1D,file2D,file3D;
  file1D << "higgsCombineSCAN1D" << channelstr << "Fast.MultiDimFit.mH125.8.root";
  file2D << "higgsCombineSCAN2D" << channelstr << "Fast.MultiDimFit.mH125.8.root";
  file3D << "higgsCombineSCAN3D" << channelstr << "Fast.MultiDimFit.mH125.8.root";

  TFile *fit1D = TFile::Open(file1D.str().c_str());
  TTree *tree1D = (TTree*)fit1D->Get("limit");

  TFile *fit2D = TFile::Open(file2D.str().c_str());
  TTree *tree2D = (TTree*)fit2D->Get("limit");

  TFile *fit3D = TFile::Open(file3D.str().c_str());
  TTree *tree3D = (TTree*)fit3D->Get("limit");

  TGraph *g1D = new TGraph();
  TGraph *g2D = new TGraph();
  TGraph *g3D = new TGraph();

  vector<TTree*> trees;
  trees.push_back(tree1D);
  trees.push_back(tree2D);
  trees.push_back(tree3D);

  vector<TGraph*> graphs;
  graphs.push_back(g1D);
  graphs.push_back(g2D);
  graphs.push_back(g3D);

  for(int dim=0; dim<(int)trees.size(); ++dim) {
    
    cout << "Analyzing scan for n. dim. fit = " << dim << endl;

    float MH;
    float deltaNLL;
    trees[dim]->SetBranchAddress("MH", &MH);
    trees[dim]->SetBranchAddress("deltaNLL", &deltaNLL);

    graphs[dim]->SetLineWidth(2);

    for(int i=1; i<(int)trees[dim]->GetEntries();++i) {
       trees[dim]->GetEntry(i);
       graphs[dim]->SetPoint(i-1,MH,2*deltaNLL);
    }

  }

  graphs[0]->SetLineColor(kRed);
  graphs[1]->SetLineColor(kGreen+2);
  graphs[2]->SetLineColor(kBlue);

  TCanvas *c1 = new TCanvas("c1","",750,750);
  graphs[0]->GetXaxis()->SetRangeUser(122,130);
  graphs[0]->GetYaxis()->SetRangeUser(0,10);
  graphs[0]->GetXaxis()->SetTitle("m_{X} (GeV)");
  graphs[0]->GetYaxis()->SetTitle("-2 #Delta lnL");

  graphs[0]->Draw("al");
  for(int i=1;i<3;++i) graphs[i]->Draw("l");

  // draw the legend
  TLegend *legend = new TLegend(0.75,0.65,0.90,0.90,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);
  legend->AddEntry(graphs[0], "1D fit","l");
  legend->AddEntry(graphs[1], "2D fit","l");
  legend->AddEntry(graphs[2], "3D fit","l");
  legend->Draw();

  TLine *line1 = new TLine(122,1,130,1);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(3.0);
  line1->Draw();

  TLine *line2 = new TLine(122,4,130,4);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(1.5);
  line2->Draw();

  TLatex *CP = CMSPreliminary();
  CP->Draw();

  TPaveText *comment;
  if(channel==0) comment = text("H #rightarrow ZZ #rightarrow 4#mu",0.20,0.90,0.40,0.90);
  if(channel==1) comment = text("H #rightarrow ZZ #rightarrow 4e",0.20,0.90,0.40,0.90);
  if(channel==2) comment = text("H #rightarrow ZZ #rightarrow 2e2#mu",0.20,0.90,0.40,0.90);
  comment->Draw();


  stringstream outnamepdf;
  outnamepdf << "scanmass_bydim_" << channelstr << ".pdf";
  c1->SaveAs(outnamepdf.str().c_str());

  stringstream outnamepng;
  outnamepng << "scanmass_bydim_" << channelstr << ".png";
  c1->SaveAs(outnamepng.str().c_str());

}


    
void cccPlot(int ndim) {

  float fitval[4], fiterrl[4], fiterrh[4];
  for(int cha=0; cha<4; ++cha) {
    fiterrl[cha]=0;
    fiterrh[cha]=0;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  
  TStyle *mystyle = RooHZZStyle("ZZ");
  mystyle->cd();

  stringstream file4mu,file4e,file2e2mu,filecomb;
  file4mu   << "higgsCombineSCAN" << ndim << "D4muFast.MultiDimFit.mH125.8.root";
  file4e    << "higgsCombineSCAN" << ndim << "D4eFast.MultiDimFit.mH125.8.root";
  file2e2mu << "higgsCombineSCAN" << ndim << "D2e2muFast.MultiDimFit.mH125.8.root";
  filecomb  << "higgsCombineSCAN" << ndim << "DcombFast.MultiDimFit.mH125.8.root";

  TFile *fit4mu = TFile::Open(file4mu.str().c_str());
  TTree *tree4mu = (TTree*)fit4mu->Get("limit");

  TFile *fit4e = TFile::Open(file4e.str().c_str());
  TTree *tree4e = (TTree*)fit4e->Get("limit");

  TFile *fit2e2mu = TFile::Open(file2e2mu.str().c_str());
  TTree *tree2e2mu = (TTree*)fit2e2mu->Get("limit");

  TFile *fitcomb = TFile::Open(filecomb.str().c_str());
  TTree *treecomb = (TTree*)fitcomb->Get("limit");

  vector<TTree*> trees;
  trees.push_back(tree4mu);
  trees.push_back(tree4e);
  trees.push_back(tree2e2mu);
  trees.push_back(treecomb);

  for(int cha=0; cha<(int)trees.size(); ++cha) {
    
    cout << "Analyzing scan for channel = " << cha << endl;

    float MH;
    float deltaNLL;
    trees[cha]->SetBranchAddress("MH", &MH);
    trees[cha]->SetBranchAddress("deltaNLL", &deltaNLL);
    
    bool firstcross=false;
    bool secondcross=false;
    float prevStepMH=-1;
    for(int i=0; i<(int)trees[cha]->GetEntries();++i) {
       trees[cha]->GetEntry(i);
       if(i==0) fitval[cha]=MH;
       else {
	 if(2*deltaNLL<1 && firstcross==false) {
	   fiterrl[cha]=MH;
	   firstcross=true;
	 }
	 if(2*deltaNLL>1 && firstcross==true && secondcross==false) {
	   fiterrh[cha]=prevStepMH;
	   secondcross=true;
	 }
	 prevStepMH=MH;
       }
    }

  }

  for(int cha=0; cha<4; ++cha) {
    fiterrl[cha]=fitval[cha]-fiterrl[cha];
    fiterrh[cha]=fiterrh[cha]-fitval[cha];
    // patch if the scan arrested too early
    if(fiterrh[cha]==(-fitval[cha])) fiterrh[cha]=fiterrl[cha];
  }


    TLatex l; l.SetTextFont(43); l.SetNDC(); l.SetTextSize(25);

    TCanvas *c1 = new TCanvas("c1","",750,750);
    c1->SetLeftMargin(0.4);
    c1->SetGridx(1);

    int nChann = 3;
    TH2F frame("frame",";best fit m_{X} (GeV);",1,122,132,nChann,0,nChann);

    TGraphAsymmErrors points(nChann);
    for (int cha=0; cha<3; ++cha) {
      TString channame("");
      if (cha==0) channame+=" 4#mu";
      if (cha==1) channame+=" 4e";
      if (cha==2) channame+=" 2e2#mu";
      points.SetPoint(cha,       fitval[cha],  cha+0.5);
      points.SetPointError(cha,  fiterrl[cha], fiterrh[cha], 0, 0);
      frame.GetYaxis()->SetBinLabel(cha+1, channame);
    }
    points.SetLineColor(kRed);
    points.SetLineWidth(3);
    points.SetMarkerStyle(21);
    frame.GetXaxis()->SetNdivisions(5,kFALSE);
    frame.GetXaxis()->SetTitleSize(0.05);
    frame.GetXaxis()->SetLabelSize(0.04);
    frame.GetYaxis()->SetLabelSize(0.06);
    frame.Draw(); gStyle->SetOptStat(0);
    TBox globalFitBand(fitval[3]-fiterrl[3], 0, fitval[3]+fiterrh[3], nChann);
    globalFitBand.SetFillStyle(3013);
    globalFitBand.SetFillColor(65);
    globalFitBand.SetLineStyle(0);
    globalFitBand.DrawClone();
    TLine globalFitLine(fitval[3], 0, fitval[3], nChann);
    globalFitLine.SetLineWidth(4);
    globalFitLine.SetLineColor(214);
    globalFitLine.DrawClone();
    points.Draw("P SAME");
    if(ndim==1) l.DrawLatex(0.45, 0.90, Form("H #rightarrow ZZ, 1D"));
    if(ndim==2) l.DrawLatex(0.45, 0.90, Form("H #rightarrow ZZ, 2D"));
    if(ndim==3) l.DrawLatex(0.45, 0.90, Form("H #rightarrow ZZ, 3D"));

    stringstream outnamepdf;
    outnamepdf << "bestfit_bychannel_" << ndim << "D.pdf";
    c1->SaveAs(outnamepdf.str().c_str());
    stringstream outnamepng;
    outnamepng << "bestfit_bychannel_" << ndim << "D.png";
    c1->SaveAs(outnamepng.str().c_str());

}
