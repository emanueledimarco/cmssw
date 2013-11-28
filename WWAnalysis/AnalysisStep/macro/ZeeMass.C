#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;


int plotMass(const char* cut, char *filename) {
  
  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  gStyle->SetOptStat("rme");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  float min=75;
  float max=105;
  int nbins=70;
  float step=(max-min)/float(nbins);
  char stepstr[20];
  sprintf(stepstr,"Events/%4.2f GeV",step);

  TFile *fileData = TFile::Open("results/hzzTree.root");
  TFile *fileDY = TFile::Open("results/hzzTree_id111.root");

  TTree *treeData = (TTree*)fileData->Get("zeetree/probe_tree");
  TTree *treeDY = (TTree*)fileDY->Get("zeetree/probe_tree");
  treeDY->AddFriend( "zeetree/pu_probe_tree = zeetree/probe_tree", "results/hzzTree_id111_puw.root" );

  TH1D *massDYsmear = new TH1D("massDYsmear","",nbins,min,max);
  TH1D *massData = new TH1D("massData","",nbins,min,max);

  // cosmetics
  massDYsmear->Sumw2();
  massDYsmear->SetFillColor(kCyan-9);
  massDYsmear->SetFillStyle(3003);
  massDYsmear->SetLineColor(kAzure+3);
  massDYsmear->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
  massDYsmear->GetYaxis()->SetTitle(stepstr);
  massData->SetMarkerSize(0.8);
  massData->SetMarkerStyle(21);
  massData->SetMaximum(massData->GetMaximum()+0.1*massData->GetMaximum());
  
  TString finalcutData = TString("(")+TString(cut)+TString(")*1");
  TString finalcutMC = TString("(")+TString(cut)+TString(")*puW");
  
  treeData->Project("massData","zmll",finalcutData);
  treeDY->Project("massDYsmear","zmll",finalcutMC);

  massDYsmear->Scale(massData->Integral()/massDYsmear->Integral());

  TLegend* leg = new TLegend(0.2,0.7,0.4,0.89);
  leg->SetBorderSize(     0);
  leg->SetFillColor (     0);
  leg->SetTextAlign (    12);
  leg->SetTextFont  (    42);
  leg->SetTextSize  (0.05);
  leg->SetFillColor(0);
  leg->AddEntry(massDYsmear,"MC smeared","F");
  leg->AddEntry(massData,"Data","pe");

  TCanvas *c1 = new TCanvas("c1","c1");
  massDYsmear->Draw("hist");
  massData->Draw("same pe1");
  leg->Draw();

  c1->SaveAs(filename);

  fileData->Close();
  fileDY->Close();

  return 0;
 
}

int plotAll() {
  
  vector<string> cutsele1;
  cutsele1.push_back("abs(l1eta)<1.479");
  cutsele1.push_back("abs(l1eta)>1.479");
  cutsele1.push_back("abs(l1eta)<1 && l1bdtID>0.9");
  cutsele1.push_back("abs(l1eta)<1 && l1bdtID<0.9");
  cutsele1.push_back("abs(l1eta)>1 && abs(l1eta)<1.479 && l1bdtID>0.9");
  cutsele1.push_back("abs(l1eta)>1 && abs(l1eta)<1.479 && l1bdtID<0.9");
  cutsele1.push_back("abs(l1eta)>1.479 && abs(l1eta)<2 && l1bdtID>0.9");
  cutsele1.push_back("abs(l1eta)>1.479 && abs(l1eta)<2 && l1bdtID<0.9");
  cutsele1.push_back("abs(l1eta)>2 && l1bdtID>0.9");
  cutsele1.push_back("abs(l1eta)>2 && l1bdtID<0.9");

  vector<string> cutsele2;
  cutsele2.push_back("abs(l2eta)<1.479");
  cutsele2.push_back("abs(l2eta)>1.479");
  cutsele2.push_back("abs(l2eta)<1 && l2bdtID>0.9");
  cutsele2.push_back("abs(l2eta)<1 && l2bdtID<0.9");
  cutsele2.push_back("abs(l2eta)>1 && abs(l2eta)<1.479 && l2bdtID>0.9");
  cutsele2.push_back("abs(l2eta)>1 && abs(l2eta)<1.479 && l2bdtID<0.9");
  cutsele2.push_back("abs(l2eta)>1.479 && abs(l2eta)<2 && l2bdtID>0.9");
  cutsele2.push_back("abs(l2eta)>1.479 && abs(l2eta)<2 && l2bdtID<0.9");
  cutsele2.push_back("abs(l2eta)>2 && l2bdtID>0.9");
  cutsele2.push_back("abs(l2eta)>2 && l2bdtID<0.9");

  for(int l1cut=0;l1cut<(int)cutsele1.size();++l1cut) {
    for(int l2cut=0;l2cut<(int)cutsele2.size();++l2cut) {
      string cut=cutsele1[l1cut]+string(" && ")+cutsele2[l2cut];
      cout << "plotting cut = " << cut << endl;
      char filename[100];
      sprintf(filename,"mass_zee_%d_%d.png",l1cut,l2cut);
      plotMass(cut.c_str(),filename);
    }
  }

  return 0;

}
