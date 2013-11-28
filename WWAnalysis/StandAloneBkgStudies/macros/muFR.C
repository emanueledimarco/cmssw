#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCut.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <sstream>
#include <math.h>


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
   }


void makeFakeRatePlot(TTree* treeZ, TTree* treeW, TTree* treeQCD10, TTree* treeQCD15, TTree* treeData,
		      double scale1,double scale2,double scale3, double scale4, 
		      float* ptBin,int ptBinsSize,
		      string selection, string variable,
		      string labelX);
	
void makeFakeRatePlot2D(TTree* tree1,TTree* tree2,TTree* tree3a,TTree* tree3b,TTree* tree4,
			double scale1,double scale2,double scale3a,double scale3b,
			float* binsX,int binsSizeX,float* binsY,int binsSizeY,
			string selection,string variables);

void makeStackPlot(TTree* tree1,TTree* tree2,TTree* tree3a,TTree* tree3b,TTree* tree4,
		   double scale1,double scale2,double scale3a,double scale3b,
		   int nbins,float binMin,float binMax,
		   string selection, string variable,
		   string labelX,bool isLog=true);


void setHistoAttributes(TH1D* h,
			int lCol,int lWidth,float max, float min,string xTitle, string yTitle,string Title);	      

void muFRtest2D() {

  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleOffset(0.80,"X");
  gStyle->SetTitleOffset(0.80,"Y");
  // ------------------------------ 


  // yes MC triggers
  TFile* f1 = new TFile("new/tree_004.root");
  TFile* f2 = new TFile("new/tree_026.root");
  TFile* f3a = new TFile("new/tree_101.root");
  TFile* f3b = new TFile("new/tree_102.root");
  TFile* f4 = new TFile("new/tree_151_data.root");

  double scale1 = 0.8396525673;
  double scale2 = 2.0722688028;
  double scale3a = 43.5174874180;
  double scale3b = 2.8768663179;
    
  //scale to 1 /pb
  scale1 *= 1./1000.;
  scale2 *= 1./1000.;
  scale3a *= 1./1000.;
  scale3b *= 1./1000.;
  
  
  string commonSelection = "dR>1.0 && jetPt >20 && mt < 20 && pfmet < 20";
  string zVeto = "(!(mass>0 && mass<12) && !(mass>0 && abs(mass-91.1876) < 15) )";
  stringstream selectionStream;
  selectionStream << commonSelection << " && " << zVeto;
  
  string variableToPlot("");
  int nbins; float binMin; float binMax;

  TTree* tree1 = (TTree*)f1->Get("treeL/probe_tree");
  TTree* tree2 = (TTree*)f2->Get("treeL/probe_tree");
  TTree* tree3a = (TTree*)f3a->Get("treeL/probe_tree");
  TTree* tree3b = (TTree*)f3b->Get("treeL/probe_tree");
  TTree* tree4 = (TTree*)f4->Get("treeL/probe_tree");


  //====================================
  const int ptBinsSize(5);
  float ptBins[ptBinsSize+1];
  ptBins[0] = 10;
  ptBins[1] = 14;
  ptBins[2] = 20;
  ptBins[3] = 26;
  ptBins[4] = 50;
  ptBins[5] = 100;

  const int etaBinsSize(3);
  float etaBins[etaBinsSize+1];
  etaBins[0] = -2.5;
  etaBins[1] = -1.5;
  etaBins[2] =  1.5;
  etaBins[3] =  2.5;
  
  
  const int phiBinsSize(10);
  float phiBins[phiBinsSize+1];
  float rangePhi = 2*M_PI;
  float binSizePhi = rangePhi/10.;
  for(unsigned int i=0;i<=10;i++){
    phiBins[i] = -M_PI + i*binSizePhi;
  }
  
  
  //====================================
  string labelX("");
  string ptLabelX("pT (GeV)");
  string etaLabelX("eta");
  string phiLabelX("phi");

  //------  
  string varPt("muPt");
  makeFakeRatePlot(tree1,tree2,tree3a,tree3b,tree4,
		   scale1,scale2,scale3a,scale3b,
		   ptBins,ptBinsSize,
		   selectionStream.str(), varPt,
		   ptLabelX);
  
 
  //------  
  string varEta("muEta");
  makeFakeRatePlot(tree1,tree2,tree3a,tree3b,tree4,
		   scale1,scale2,scale3a,scale3b,
		   etaBins,etaBinsSize,
		   selectionStream.str(), varEta,
		   etaLabelX);


  //------
  string varPhi("muPhi");
  makeFakeRatePlot(tree1,tree2,tree3a,tree3b,tree4,
		   scale1,scale2,scale3a,scale3b,
		   phiBins,phiBinsSize,
		   selectionStream.str(), varPhi,
		   phiLabelX);
  
		   
  string variables("muEta:muPt");
  makeFakeRatePlot2D(tree1,tree2,tree3a,tree3b,tree4,
		     scale1,scale2,scale3a,scale3b,
		     ptBins,ptBinsSize,etaBins,etaBinsSize,
		     selectionStream.str(), variables);


  
  // ---------  
  variableToPlot = "muPt";
  nbins = 90;
  binMin = 10.;
  binMax = 100.;
  makeStackPlot(tree1,tree2,tree3a,tree3b,tree4,
		scale1,scale2,scale3a,scale3b,
		nbins,binMin,binMax,
		selectionStream.str(),variableToPlot,
		ptLabelX);  



  
  stringstream plusMuPtSelection;
  plusMuPtSelection << commonSelection << " && lMuPt<70";
  
  variableToPlot = "jetPt";
  labelX = "Jet pT (GeV)";
  nbins = 10;
  binMin = 10.;
  binMax = 160.;
  makeStackPlot(tree1,tree2,tree3a,tree3b,tree4,		   
		scale1,scale2,scale3a,scale3b,
		nbins,binMin,binMax,
		plusMuPtSelection.str(),variableToPlot,
		labelX);


  variableToPlot = "mass";
  labelX = "mass_ll (GeV)";
  nbins = 35;
  binMin = 70.;
  binMax = 140.;
  makeStackPlot(tree1,tree2,tree3a,tree3b,tree4,		   
		scale1,scale2,scale3a,scale3b,
		nbins,binMin,binMax,
		plusMuPtSelection.str(),variableToPlot,
		labelX,false);  
  
}


   
void makeFakeRatePlot(TTree* tree1,TTree* tree2,TTree* tree3a,TTree* tree3b,TTree* tree4,
		      double scale1,double scale2,double scale3a,double scale3b,
		      float* bins,int binsSize,
		      string selection, string variable,
		      string labelX){

  TCanvas* canv = new TCanvas("canv","myCanv for FakeRate",500,750);
  canv->Divide(2,3);
  canv->cd(1);


  TH1D* hRatio1 = new TH1D("hRatio1","hRatio1",binsSize,bins);
  TH1D* hRatio2 = new TH1D("hRatio2","hRatio2",binsSize,bins);
  TH1D* hRatio3 = new TH1D("hRatio3","hRatio3",binsSize,bins);
  TH1D* hRatio4 = new TH1D("hRatio4","hRatio4",binsSize,bins);

  //mixing of all MCs
  TH1D* hRatio5 = new TH1D("hRatio5","hRatio5",binsSize,bins);

  //data, corrected for W/Z MCs 
  TH1D* hRatio6 = new TH1D("hRatio6","hRatio6",binsSize,bins);


  TH1D* hD1 = new TH1D("hD1","histo den1",binsSize,bins);
  TH1D* hD2 = new TH1D("hD2","histo den2",binsSize,bins);
  TH1D* hD3a = new TH1D("hD3a","histo den3a",binsSize,bins);
  TH1D* hD3b = new TH1D("hD3b","histo den3b",binsSize,bins);
  TH1D* hD4 = new TH1D("hD4","histo den4",binsSize,bins);
  TH1D* hD5 = new TH1D("hD5","histo den5",binsSize,bins);
  TH1D* hD6 = new TH1D("hD6","histo den6",binsSize,bins);

  TH1D* hN1 = new TH1D("hN1","histo num1",binsSize,bins);
  TH1D* hN2 = new TH1D("hN2","histo num2",binsSize,bins);
  TH1D* hN3a = new TH1D("hN3a","histo num3a",binsSize,bins);
  TH1D* hN3b = new TH1D("hN3b","histo num3b",binsSize,bins);
  TH1D* hN4 = new TH1D("hN4","histo num4",binsSize,bins);
  TH1D* hN5 = new TH1D("hN5","histo num5",binsSize,bins);
  TH1D* hN6 = new TH1D("hN6","histo num6",binsSize,bins);


    
  stringstream streamD1; streamD1 << variable << ">>hD1";
  stringstream streamD2; streamD2 << variable << ">>hD2";
  stringstream streamD3a; streamD3a << variable << ">>hD3a";
  stringstream streamD3b; streamD3b << variable << ">>hD3b";
  stringstream streamD4; streamD4 << variable << ">>hD4";

  stringstream streamN1; streamN1 << variable << ">>hN1";
  stringstream streamN2; streamN2 << variable << ">>hN2";
  stringstream streamN3a; streamN3a << variable << ">>hN3a";
  stringstream streamN3b; streamN3b << variable << ">>hN3b";
  stringstream streamN4; streamN4 << variable << ">>hN4";

  stringstream selectionLoose;
  stringstream selectionTight;

  selectionLoose << "weight*(" << selection << " && isLoose)" ;
  selectionTight << "weight*(" << selection << " && isLoose && isTight)" ;

  hD1->Sumw2();
  hD2->Sumw2();
  hD3a->Sumw2();
  hD3b->Sumw2();
  hD4->Sumw2();

  tree1->Draw(streamD1.str().c_str(),selectionLoose.str().c_str());
  tree2->Draw(streamD2.str().c_str(),selectionLoose.str().c_str());
  tree3a->Draw(streamD3a.str().c_str(),selectionLoose.str().c_str());
  tree3b->Draw(streamD3b.str().c_str(),selectionLoose.str().c_str());
  tree4->Draw(streamD4.str().c_str(),selectionLoose.str().c_str());

  hN1->Sumw2();
  hN2->Sumw2();
  hN3a->Sumw2();
  hN3b->Sumw2();
  hN4->Sumw2();
  tree1->Draw(streamN1.str().c_str(),selectionTight.str().c_str());
  tree2->Draw(streamN2.str().c_str(),selectionTight.str().c_str());
  tree3a->Draw(streamN3a.str().c_str(),selectionTight.str().c_str());
  tree3b->Draw(streamN3b.str().c_str(),selectionTight.str().c_str());
  tree4->Draw(streamN4.str().c_str(),selectionTight.str().c_str());
  
  hN1->Scale(scale1); 
  hN2->Scale(scale2); 

  hN3a->Scale(scale3a); 
  hN3b->Scale(scale3b); 
  hN3a->Add(hN3b);

  hN5->Add(hN1);
  hN5->Add(hN2);
  hN5->Add(hN3a);

  hN6->Add(hN4);//adding data
  hN6->Add(hN1,-1.);//subtracting MC
  hN6->Add(hN2,-1.);//subtracting MC

  hD1->Scale(scale1);   
  hD2->Scale(scale2); 

  hD3a->Scale(scale3a); 
  hD3b->Scale(scale3b); 
  hD3a->Add(hD3b);

  hD5->Add(hD1);
  hD5->Add(hD2);
  hD5->Add(hD3a);

  hD6->Add(hD4);//adding data
  hD6->Add(hD1,-1.);//subtracting MC
  hD6->Add(hD2,-1.);//subtracting MC

  /*
  canv->cd(1); hN1->Draw();   canv->cd(2); hN2->Draw();
  canv->cd(3); hN3a->Draw();   canv->cd(4); hN4->Draw();
  canv->cd(5); hN5->Draw();   canv->cd(6); hN6->Draw();
  gPad->Update(); Wait();


  canv->cd(1); hD1->Draw();   canv->cd(2); hD2->Draw();
  canv->cd(3); hD3a->Draw();   canv->cd(4); hD4->Draw();
  canv->cd(5); hD5->Draw();   canv->cd(6); hD6->Draw();
  gPad->Update(); Wait();
  */


  hRatio1->Divide(hN1,hD1,1.,1.,"B");
  hRatio2->Divide(hN2,hD2,1.,1.,"B");
  hRatio3->Divide(hN3a,hD3a,1.,1.,"B");
  hRatio4->Divide(hN4,hD4,1.,1.,"B");
  hRatio5->Divide(hN5,hD5,1.,1.,"B");
  hRatio6->Divide(hN6,hD6,1.,1.,"B");

  string yTitle("#LooseAndTight / #Loose");
  string xTitle(labelX);

  setHistoAttributes(hRatio1,2,2,1.01,0.0,xTitle,yTitle,"ZJets MC");
  setHistoAttributes(hRatio2,2,2,1.01,0.0,xTitle,yTitle,"WJets MC");
  setHistoAttributes(hRatio3,2,2,1.01,0.0,xTitle,yTitle,"QCD MC");
  setHistoAttributes(hRatio4,1,2,1.01,0.0,xTitle,yTitle,"");
  setHistoAttributes(hRatio5,2,2,1.01,0.0,xTitle,yTitle,"");
  setHistoAttributes(hRatio6,1,2,1.01,0.0,xTitle,yTitle,"");

  hRatio5->SetFillColor(2);
  hRatio3->SetFillColor(2);

  //canv->Divide(2,3);
  canv->cd(1);  hRatio1->Draw();
  canv->cd(2);  hRatio2->Draw();
  canv->cd(3);  hRatio3->Draw();


  canv->cd(5);  hRatio5->Draw("hist"); hRatio4->Draw("same P");  
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(hRatio5,"Z+W+QCD MCs","F");
  leg->AddEntry(hRatio4,"May10ReReco data","L");
  leg->Draw();


  canv->cd(6);  
  TH1D* hRatio3Clone = (TH1D*) hRatio3->Clone(); 
  hRatio3Clone->SetTitle("");
  hRatio3Clone->SetFillColor(2);
  hRatio3Clone->Draw("hist"); hRatio6->Draw("same P");  
  TLegend* leg2 = new TLegend(0.5,0.7,0.9,0.89);
  leg2->SetFillColor(0);
  leg2->AddEntry(hRatio3,"QCD MC","F");
  leg2->AddEntry(hRatio4,"Data, W/Z corrected","L");
  leg2->Draw();


  gPad->Update();
  stringstream fileNameOut; fileNameOut << "fakeRate_" << variable << ".pdf";
  //c1->Print("fakeRate.pdf"); Wait();
  canv->Print(fileNameOut.str().c_str()); Wait();
  delete canv;

}


void setHistoAttributes(TH1D* h,
			int lCol,int lWidth,float max, float min,string xTitle, string yTitle,string title){
  
  h->SetLineColor(lCol);
  h->SetLineWidth(lWidth);
  h->SetMaximum(max);
  h->SetMinimum(min);

  h->GetXaxis()->SetTitle(xTitle.c_str());
  h->GetYaxis()->SetTitle(yTitle.c_str());
  h->SetTitle(title.c_str());



}	      





void makeStackPlot(TTree* tree1,TTree* tree2,TTree* tree3a,TTree* tree3b,TTree* tree4,
		   double scale1,double scale2,double scale3a,double scale3b,
		   int nbins,float binMin,float binMax,
		   string selection, string variable,
		   string labelX,bool isLog){

  TCanvas* canv = new TCanvas("canv","myCanv for Stack Plot",500,500);

  
  canv->SetLogy(isLog);
  canv->cd();

  THStack *hs2 = new THStack("hs2","Stacked 1D histograms");

  TH1F* hpt1 = new TH1F("hpt1","histo pt1",nbins,binMin,binMax);
  TH1F* hpt2 = new TH1F("hpt2","histo pt2",nbins,binMin,binMax);
  TH1F* hpt3a = new TH1F("hpt3a","histo pt3a",nbins,binMin,binMax);
  TH1F* hpt3b = new TH1F("hpt3b","histo pt3b",nbins,binMin,binMax);
  TH1F* hpt4 = new TH1F("hpt4","histo pt4",nbins,binMin,binMax);

  stringstream stream1; stream1 << variable << ">>hpt1";
  stringstream stream2; stream2 << variable << ">>hpt2";
  stringstream stream3a; stream3a << variable << ">>hpt3a";
  stringstream stream3b; stream3b << variable << ">>hpt3b";
  stringstream stream4; stream4 << variable << ">>hpt4";

  stringstream selectionLoose;
  stringstream selectionTight;

  selectionLoose << "weight*(" << selection << " && isLoose)" ;
  selectionTight << "weight*(" << selection << " && isLoose && isTight)" ;

  hpt4->Sumw2();

  tree1->Draw(stream1.str().c_str(),selectionLoose.str().c_str());
  tree2->Draw(stream2.str().c_str(),selectionLoose.str().c_str());
  tree3a->Draw(stream3a.str().c_str(),selectionLoose.str().c_str());
  tree3b->Draw(stream3b.str().c_str(),selectionLoose.str().c_str());
  tree4->Draw(stream4.str().c_str(),selectionLoose.str().c_str());


  hpt1->SetFillColor(kBlue);
  hpt2->SetFillColor(kRed);
  hpt3a->SetFillColor(kGreen);
  hpt3b->SetFillColor(kGreen);
  hpt4->SetLineWidth(2);

  hpt1->Scale(scale1);
  hpt2->Scale(scale2);

  hpt3a->Scale(scale3a);
  hpt3b->Scale(scale3b);
  hpt3a->Add(hpt3b);


  hs2->Add(hpt1);
  hs2->Add(hpt2);
  hs2->Add(hpt3a);

 
  float maximum = max(hpt4->GetMaximum(),hs2->GetMaximum());
 
  hpt4->SetMaximum(maximum*1.10);
  hs2->SetMaximum(maximum*1.10);
  hs2->Draw("");
  
  hpt4->Draw("same p");


  hs2->GetXaxis()->SetTitle(labelX.c_str());
  
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(hpt1,"ZJets","F");
  leg->AddEntry(hpt2,"WJets","F");
  leg->AddEntry(hpt3a,"QCD","F");
  leg->AddEntry(hpt4,"Data","L");
  leg->Draw();
  //hpt4->Scale(0.5);
  //hpt4->Draw("same");
  gPad->Update();
  stringstream nameFileOut; 
  nameFileOut << "stack_" << variable << ".pdf";
  gPad->Print(nameFileOut.str().c_str()); Wait();
  delete canv;
  
  
}


void makeFakeRatePlot2D(TTree* tree1,TTree* tree2,TTree* tree3a,TTree* tree3b,TTree* tree4,
			double scale1,double scale2,double scale3a,double scale3b,
			float* binsX,int binsSizeX,float* binsY,int binsSizeY,
			string selection,string variables)
{
  TCanvas* canv = new TCanvas("canv","myCanv for FakeRate",750,500);
  canv->Divide(1,2);
  //canv->SetLogx(1);
  canv->cd(1);


  TH2D* hRatio1 = new TH2D("hRatio1","hRatio1",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hRatio2 = new TH2D("hRatio2","hRatio2",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hRatio3 = new TH2D("hRatio3","hRatio3",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hRatio4 = new TH2D("hRatio4","hRatio4",binsSizeX,binsX,binsSizeY,binsY);

  //mixing of all MCs
  TH2D* hRatio5 = new TH2D("hRatio5","hRatio5",binsSizeX,binsX,binsSizeY,binsY);

  //data, corrected for W/Z MCs 
  TH2D* hRatio6 = new TH2D("hRatio6","hRatio6",binsSizeX,binsX,binsSizeY,binsY);


  TH2D* hD1 = new TH2D("hD1","histo den1",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hD2 = new TH2D("hD2","histo den2",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hD3a = new TH2D("hD3a","histo den3a",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hD3b = new TH2D("hD3b","histo den3b",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hD4 = new TH2D("hD4","histo den4",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hD5 = new TH2D("hD5","histo den5",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hD6 = new TH2D("hD6","histo den6",binsSizeX,binsX,binsSizeY,binsY);

  TH2D* hN1 = new TH2D("hN1","histo num1",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hN2 = new TH2D("hN2","histo num2",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hN3a = new TH2D("hN3a","histo num3a",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hN3b = new TH2D("hN3b","histo num3b",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hN4 = new TH2D("hN4","histo num4",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hN5 = new TH2D("hN5","histo num5",binsSizeX,binsX,binsSizeY,binsY);
  TH2D* hN6 = new TH2D("hN6","histo num6",binsSizeX,binsX,binsSizeY,binsY);


    
  stringstream streamD1; streamD1 << variables << ">>hD1";
  stringstream streamD2; streamD2 << variables << ">>hD2";
  stringstream streamD3a; streamD3a << variables << ">>hD3a";
  stringstream streamD3b; streamD3b << variables << ">>hD3b";
  stringstream streamD4; streamD4 << variables << ">>hD4";

  stringstream streamN1; streamN1 << variables << ">>hN1";
  stringstream streamN2; streamN2 << variables << ">>hN2";
  stringstream streamN3a; streamN3a << variables << ">>hN3a";
  stringstream streamN3b; streamN3b << variables << ">>hN3b";
  stringstream streamN4; streamN4 << variables << ">>hN4";

  stringstream selectionLoose;
  stringstream selectionTight;

  selectionLoose << "weight*(" << selection << " && isLoose)" ;
  selectionTight << "weight*(" << selection << " && isLoose && isTight)" ;

  hD1->Sumw2();
  hD2->Sumw2();
  hD3a->Sumw2();
  hD3b->Sumw2();
  hD4->Sumw2();

  tree1->Draw(streamD1.str().c_str(),selectionLoose.str().c_str());
  tree2->Draw(streamD2.str().c_str(),selectionLoose.str().c_str());
  tree3a->Draw(streamD3a.str().c_str(),selectionLoose.str().c_str());
  tree3b->Draw(streamD3b.str().c_str(),selectionLoose.str().c_str());
  tree4->Draw(streamD4.str().c_str(),selectionLoose.str().c_str());

  hN1->Sumw2();
  hN2->Sumw2();
  hN3a->Sumw2();
  hN3b->Sumw2();
  hN4->Sumw2();
  tree1->Draw(streamN1.str().c_str(),selectionTight.str().c_str());
  tree2->Draw(streamN2.str().c_str(),selectionTight.str().c_str());
  tree3a->Draw(streamN3a.str().c_str(),selectionTight.str().c_str());
  tree3b->Draw(streamN3b.str().c_str(),selectionTight.str().c_str());
  tree4->Draw(streamN4.str().c_str(),selectionTight.str().c_str());
  
  hN1->Scale(scale1); 
  hN2->Scale(scale2); 

  hN3a->Scale(scale3a); 
  hN3b->Scale(scale3b); 
  hN3a->Add(hN3b);

  hN5->Add(hN1);
  hN5->Add(hN2);
  hN5->Add(hN3a);

  hN6->Add(hN4);//adding data
  hN6->Add(hN1,-1.);//subtracting MC
  hN6->Add(hN2,-1.);//subtracting MC

  hD1->Scale(scale1);   
  hD2->Scale(scale2); 

  hD3a->Scale(scale3a); 
  hD3b->Scale(scale3b); 
  hD3a->Add(hD3b);

  hD5->Add(hD1);
  hD5->Add(hD2);
  hD5->Add(hD3a);

  hD6->Add(hD4);//adding data
  hD6->Add(hD1,-1.);//subtracting MC
  hD6->Add(hD2,-1.);//subtracting MC

  /*
  canv->cd(1); hN1->Draw("text");   canv->cd(2); hN2->Draw("text");
  canv->cd(3); hN3a->Draw("text");   canv->cd(4); hN4->Draw("text");
  canv->cd(5); hN5->Draw("text");   canv->cd(6); hN6->Draw("text");
  hN1->Print();
  gPad->Update(); Wait();


  canv->cd(1); hD1->Draw("text");   canv->cd(2); hD2->Draw("text");
  canv->cd(3); hD3a->Draw("text");   canv->cd(4); hD4->Draw("text");
  canv->cd(5); hD5->Draw("text");   canv->cd(6); hD6->Draw("text");
  gPad->Update(); Wait();
  */


  hRatio1->Divide(hN1,hD1,1.,1.,"B");
  hRatio2->Divide(hN2,hD2,1.,1.,"B");
  hRatio3->Divide(hN3a,hD3a,1.,1.,"B");
  hRatio4->Divide(hN4,hD4,1.,1.,"B");
  hRatio5->Divide(hN5,hD5,1.,1.,"B");
  hRatio6->Divide(hN6,hD6,1.,1.,"B");

  //string yTitle("#LooseAndTight / #Loose");
  //string xTitle(labelX);

  /*
  setHistoAttributes(hRatio1,2,2,1.01,0.0,xTitle,yTitle,"ZJets MC");
  setHistoAttributes(hRatio2,2,2,1.01,0.0,xTitle,yTitle,"WJets MC");
  setHistoAttributes(hRatio3,2,2,1.01,0.0,xTitle,yTitle,"QCD MC");
  setHistoAttributes(hRatio4,1,2,1.01,0.0,xTitle,yTitle,"");
  setHistoAttributes(hRatio5,2,2,1.01,0.0,xTitle,yTitle,"");
  setHistoAttributes(hRatio6,1,2,1.01,0.0,xTitle,yTitle,"");
  */

  hRatio5->SetFillColor(2);
  hRatio3->SetFillColor(2);

  //canv->Divide(2,3);
  //canv->cd(1);  hRatio1->Draw("text");
  //canv->cd(2);  hRatio2->Draw("text");
  //canv->cd(3);  hRatio3->Draw("text");

  //canv->cd(5);  hRatio5->Draw("hist"); hRatio4->Draw("same P");  
  hRatio4->SetMarkerSize(2.8);
  canv->cd(1);  hRatio4->Draw("text45");  
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(hRatio5,"Z+W+QCD MCs","F");
  leg->AddEntry(hRatio4,"May10ReReco data","L");
  //leg->Draw();


  canv->cd(2);  
  TH2D* hRatio3Clone = (TH2D*) hRatio3->Clone(); 
  hRatio3Clone->SetTitle("");
  hRatio3Clone->SetFillColor(2);
  //hRatio3Clone->Draw("hist"); hRatio6->Draw("same P");  
  hRatio6->SetMarkerSize(2.8);
  hRatio6->Draw("text45");  
  TLegend* leg2 = new TLegend(0.5,0.7,0.9,0.89);
  leg2->SetFillColor(0);
  leg2->AddEntry(hRatio3,"QCD MC","F");
  leg2->AddEntry(hRatio4,"Data, W/Z corrected","L");
  //leg2->Draw();


  gPad->Update();
  stringstream fileNameOut; fileNameOut << "fakeRate_2D.pdf";
  //c1->Print("fakeRate.pdf"); Wait();
  canv->Print(fileNameOut.str().c_str()); Wait();

}
