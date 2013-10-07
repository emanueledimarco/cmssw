#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSpline.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"

using namespace RooFit;



#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

void interpolateAgain(int channel, int year=2012, bool highmass=false) {

  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(1.3);
  gStyle->SetMarkerColor(kGreen+3);
  gStyle->SetLineColor(kGreen+3);
  gStyle->SetLineWidth(2);
  gStyle->SetTitle("H#rightarrow ZZ Signal Lineshape Interpolation");

  stringstream filename;
  filename << "signalPdfs/" << ((year==2012) ? "8TeV" : "7TeV" ) << "/YesRegrYesCalibV14/" << ((highmass) ? "HighMass" : "LowMass") << "/parameters_channel" << channel << ".root";
  TFile *resultfile = TFile::Open(filename.str().c_str());
  TGraphErrors *gA1 = (TGraphErrors*)resultfile->Get("gA1");
  TGraphErrors *gA2 = (TGraphErrors*)resultfile->Get("gA2");
  TGraphErrors *gN1 = (TGraphErrors*)resultfile->Get("gN1");
  TGraphErrors *gN2 = (TGraphErrors*)resultfile->Get("gN2");
  TGraphErrors *gMeanCB = (TGraphErrors*)resultfile->Get("gMeanCB");
  TGraphErrors *gSigmaCB = (TGraphErrors*)resultfile->Get("gSigmaCB");
  TGraphErrors *gSigmaBW = 0;
  if(highmass) gSigmaBW = (TGraphErrors*)resultfile->Get("gSigmaBW");

  Int_t n = gA1->GetN();

  Double_t *x = gA1->GetX();  
  Double_t *yA1 = gA1->GetY();  Double_t *yEA1 = gA1->GetEY();
  Double_t *yA2 = gA2->GetY();  Double_t *yEA2 = gA2->GetEY();
  Double_t *yN1 = gN1->GetY();  Double_t *yEN1 = gN1->GetEY();
  Double_t *yN2 = gN2->GetY();  Double_t *yEN2 = gN2->GetEY();
  Double_t *yMeanCB = gMeanCB->GetY();  Double_t *yEMeanCB = gMeanCB->GetEY();
  Double_t *ySigmaCB = gSigmaCB->GetY();  Double_t *yESigmaCB = gSigmaCB->GetEY();
  Double_t *ySigmaBW, *yESigmaBW;
  if(highmass) { ySigmaBW = gSigmaBW->GetY();  yESigmaBW = gSigmaBW->GetEY(); }

  // skip some mass points where we did not have MC                                                   
  std::vector<int> pointstoskip;
  if(year==2012) {
    if(!highmass) {
      if(channel==0) {
	pointstoskip.push_back(325);
	pointstoskip.push_back(375);
      } 
      if(channel==1) {
	pointstoskip.push_back(180);
	pointstoskip.push_back(190);
	pointstoskip.push_back(275);
	pointstoskip.push_back(325);
      }
      if(channel==2) {
	pointstoskip.push_back(325);
	pointstoskip.push_back(350);
      }
    } else {
      if(channel==0) { 
	pointstoskip.push_back(650);
	pointstoskip.push_back(700);
	pointstoskip.push_back(900);
      } 
      if(channel==1) {
	pointstoskip.push_back(700);
	//      pointstoskip.push_back(750);
	pointstoskip.push_back(400);
	pointstoskip.push_back(800);
      }
      if(channel==2) {
	//      pointstoskip.push_back(400);
	pointstoskip.push_back(900);
	pointstoskip.push_back(1000);
      }
    }
  } else { // 7 TeV
    if(!highmass) {
      if(channel==0) {
	pointstoskip.push_back(300);
	pointstoskip.push_back(375);
      } 
      if(channel==1) {
	pointstoskip.push_back(150);
	pointstoskip.push_back(190);
	pointstoskip.push_back(220);
	pointstoskip.push_back(375);
      }
      if(channel==2) {
	pointstoskip.push_back(170);
	pointstoskip.push_back(190);
	pointstoskip.push_back(230);
	pointstoskip.push_back(375);
      }
    } else {
      if(channel==0) { 
	pointstoskip.push_back(550);
	pointstoskip.push_back(575);
	pointstoskip.push_back(600);
	pointstoskip.push_back(850);
      } 
      if(channel==1) {
      }
      if(channel==2) {
      }
    }
  }

  Int_t m = n-pointstoskip.size();
  TGraphErrors *cA1 = new TGraphErrors(m);
  TGraphErrors *cA2 = new TGraphErrors(m);
  TGraphErrors *cN1 = new TGraphErrors(m);
  TGraphErrors *cN2 = new TGraphErrors(m);
  TGraphErrors *cMeanCB = new TGraphErrors(m);
  TGraphErrors *cSigmaCB = new TGraphErrors(m);
  TGraphErrors *cSigmaBW = 0;
  if(highmass) cSigmaBW = new TGraphErrors(m);


  cA1->SetName("cA1");
  cA2->SetName("cA2");
  cN1->SetName("cN1");
  cN2->SetName("cN2");
  cMeanCB->SetName("cMeanCB");
  cSigmaCB->SetName("cSigmaCB");
  if(highmass) cSigmaBW->SetName("cSigmaBW");

  if(!highmass) { cA1->SetMinimum(0);  cA1->SetMaximum(3); }
  else          { cA1->SetMinimum(0);  cA1->SetMaximum(7); }
  if(!highmass) { cA2->SetMinimum(0);  cA2->SetMaximum(12); }
  else          { cA2->SetMinimum(0);  cA2->SetMaximum(20); }
  if(!highmass) { cN1->SetMinimum(0);  cN1->SetMaximum(20); }
  else          { cN1->SetMinimum(0);  cN1->SetMaximum(50); }
  cN2->SetMinimum(10); cN2->SetMaximum(30);
  if(!highmass) { cMeanCB->SetMinimum(-2);   cMeanCB->SetMaximum(3);    }
  else          { cMeanCB->SetMinimum(-30);  cMeanCB->SetMaximum(30);   }
  if(!highmass) { cSigmaCB->SetMinimum(0);   cSigmaCB->SetMaximum(5);   }
  else          { cSigmaCB->SetMinimum(0);   cSigmaCB->SetMaximum(150); }
  if(highmass)  { cSigmaBW->SetMinimum(0);   cSigmaBW->SetMaximum(200); }

  int k=0; 
  float xMin=7000; float xMax=-100;
  for(int i=0;i<n;++i) {
    bool skip=false;
    for(int j=0;j<(int)pointstoskip.size();++j) {
      if(x[i]==int(pointstoskip[j])) { cout << "skipping point " << i << " m(H) = " << x[i] << endl; skip=true; break; }
    }
    if(skip) continue;
    if(x[i]<xMin) xMin=x[i];
    if(x[i]>xMax) xMax=x[i];
    float penalty= (highmass) ? 10. : 1.; // for highmass, penalty term of 100 added to constrain sigmaCB<sigmaBW
    if(x[i]>=350) penalty=0.1;
    cA1->SetPoint(k,x[i],yA1[i]);            cA1->SetPointError(k,0,penalty*yEA1[i]);
    cA2->SetPoint(k,x[i],yA2[i]);            cA2->SetPointError(k,0,penalty*yEA2[i]);
    cN1->SetPoint(k,x[i],yN1[i]);            cN1->SetPointError(k,0,penalty*yEN1[i]);
    cN2->SetPoint(k,x[i],yN2[i]);            cN2->SetPointError(k,0,penalty*yEN2[i]);
    cMeanCB->SetPoint(k,x[i],yMeanCB[i]);    cMeanCB->SetPointError(k,0,penalty*yEMeanCB[i]);
    cSigmaCB->SetPoint(k,x[i],ySigmaCB[i]);  cSigmaCB->SetPointError(k,0,penalty*yESigmaCB[i]);
    if(highmass) { cSigmaBW->SetPoint(k,x[i],ySigmaBW[i]);  cSigmaBW->SetPointError(k,0,penalty*yESigmaBW[i]); }
    k++;
  }

  cout << k << "   " << m << endl;
  cout << "xMin = " << xMin << "  xMax = " << xMax << endl;

  std::stringstream tfile;
  tfile << "splines_channel" << channel << "_8TeV.root";
  TFile *outputfile = TFile::Open(tfile.str().c_str(),"recreate");

  TSpline3 *sA1 = new TSpline3("sA1",cA1);
  TSpline3 *sA2 = new TSpline3("sA2",cA2);
  TSpline3 *sN1 = new TSpline3("sN1",cN1);
  TSpline3 *sN2 = new TSpline3("sN2",cN2);
  TSpline3 *sMeanCB = new TSpline3("sMeanCB",cMeanCB);
  TSpline3 *sSigmaCB = new TSpline3("sSigmaCB",cSigmaCB);

  sA1->SetLineColor(kRed);
  sA2->SetLineColor(kRed);
  sN1->SetLineColor(kRed);
  sN2->SetLineColor(kRed);
  sMeanCB->SetLineColor(kRed);
  sSigmaCB->SetLineColor(kRed);

  sA1->SetName("sA1");
  sA2->SetName("sA2");
  sN1->SetName("sN1");
  sN2->SetName("sN2");
  sMeanCB->SetName("sMeanCB");
  sSigmaCB->SetName("sSigmaCB");

  outputfile->cd();
  sA1->Write();
  sA2->Write();
  sN1->Write();
  sN2->Write();
  sMeanCB->Write();
  sSigmaCB->Write();
  outputfile->Close();

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  cA1->GetYaxis()->SetTitle("alpha_{L}");
  cA2->GetYaxis()->SetTitle("alpha_{R}");
  cN1->GetYaxis()->SetTitle("N_{L}");
  cN2->GetYaxis()->SetTitle("N_{R}");
  cMeanCB->GetYaxis()->SetTitle("mean");
  cSigmaCB->GetYaxis()->SetTitle("#sigma");
  if(highmass) cSigmaBW->GetYaxis()->SetTitle("#Gamma BW");

  cA1->GetXaxis()->SetTitle("m_{H} [GeV]");
  cA2->GetXaxis()->SetTitle("m_{H} [GeV]");
  cN1->GetXaxis()->SetTitle("m_{H} [GeV]");
  cN2->GetXaxis()->SetTitle("m_{H} [GeV]");
  cMeanCB->GetXaxis()->SetTitle("m_{H} [GeV]");
  cSigmaCB->GetXaxis()->SetTitle("m_{H} [GeV]");
  if(highmass) cSigmaBW->GetXaxis()->SetTitle("m_{H} [GeV]");

  cA1->SetTitle("H #rightarrow ZZ Signal Lineshape Interpolation");
  cA2->SetTitle("H #rightarrow ZZ Signal Lineshape Interpolation");
  cN1->SetTitle("H #rightarrow ZZ Signal Lineshape Interpolation");
  cN2->SetTitle("H #rightarrow ZZ Signal Lineshape Interpolation");
  cMeanCB->SetTitle("H #rightarrow ZZ Signal Lineshape Interpolation");
  cSigmaCB->SetTitle("H #rightarrow ZZ Signal Lineshape Interpolation");
  if(highmass) cSigmaBW->SetTitle("H #rightarrow ZZ Signal Lineshape Interpolation");

  stringstream channame;
  channame << "_channel" << channel << ".pdf";
  if(!highmass) {
    cA1->Fit("pol5","","",xMin,xMax); cA1->Draw("Ap"); gPad->Update(); gPad->Print((string("gA1")+channame.str()).c_str());
    cA2->Fit("pol5","","",xMin,xMax); cA2->Draw("Ap"); gPad->Update(); gPad->Print((string("gA2")+channame.str()).c_str());
    cN1->Fit("pol5","","",xMin,xMax); cN1->Draw("Ap"); gPad->Update(); gPad->Print((string("gN1")+channame.str()).c_str());
    cN2->Fit("pol0","","",xMin,xMax); cN2->Draw("Ap"); gPad->Update(); gPad->Print((string("gN2")+channame.str()).c_str());
    cMeanCB->Fit("pol5","","",xMin,xMax); cMeanCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gMeanCB")+channame.str()).c_str());
    cSigmaCB->Fit("pol5","","",xMin,xMax); cSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaCB")+channame.str()).c_str());
  } else {
    cA1->Fit("pol1","","",xMin,xMax); cA1->Draw("Ap"); gPad->Update(); gPad->Print((string("gA1")+channame.str()).c_str());
    cA2->Fit("pol1","","",xMin,xMax); cA2->Draw("Ap"); gPad->Update(); gPad->Print((string("gA2")+channame.str()).c_str());
    cN1->Fit("pol3","","",xMin,xMax); cN1->Draw("Ap"); gPad->Update(); gPad->Print((string("gN1")+channame.str()).c_str());
    cN2->Fit("pol0","","",xMin,xMax); cN2->Draw("Ap"); gPad->Update(); gPad->Print((string("gN2")+channame.str()).c_str());
    cMeanCB->Fit("pol3","","",xMin,xMax); cMeanCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gMeanCB")+channame.str()).c_str());
    cSigmaCB->Fit("pol3","","",xMin,xMax); cSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaCB")+channame.str()).c_str());
    cSigmaBW->Fit("pol3","","",xMin,xMax); cSigmaBW->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaBW")+channame.str()).c_str()); 
  }


  TF1 *fA1 = (highmass) ? (TF1*)cA1->GetFunction("pol1") : (TF1*)cA1->GetFunction("pol5");
  TF1 *fA2 = (highmass) ? (TF1*)cA2->GetFunction("pol1") : (TF1*)cA2->GetFunction("pol5");
  TF1 *fN1 = (highmass) ? (TF1*)cN1->GetFunction("pol3") : (TF1*)cN1->GetFunction("pol5");
  TF1 *fN2 = (highmass) ? (TF1*)cN2->GetFunction("pol0") : (TF1*)cN2->GetFunction("pol0");
  TF1 *fMeanCB  = (highmass) ? (TF1*)cMeanCB->GetFunction("pol3")  : (TF1*)cMeanCB->GetFunction("pol5");
  TF1 *fSigmaCB = (highmass) ? (TF1*)cSigmaCB->GetFunction("pol3") : (TF1*)cSigmaCB->GetFunction("pol5");
  TF1 *fSigmaBW = (highmass) ? (TF1*)cSigmaBW->GetFunction("pol3") : 0;

  std::vector<string> names;
  std::vector<TF1*> fcns;
  fcns.push_back(fA1);       names.push_back("fA1");
  fcns.push_back(fA2);       names.push_back("fA2");
  fcns.push_back(fN1);       names.push_back("fN1");
  fcns.push_back(fN2);       names.push_back("fN2");
  fcns.push_back(fMeanCB);   names.push_back("fMeanCB");
  fcns.push_back(fSigmaCB);  names.push_back("fSigmaCB");
  if(highmass) fcns.push_back(fSigmaBW);  names.push_back("fSigmaBW");

  for(int fcn=0;fcn<(int)fcns.size();++fcn) {

    std::stringstream ss;
    
    for (int i = 0; i < fcns[fcn]->GetNumberFreeParameters(); i++) {
      if (i != 0) ss << " + (";
      else ss << "(";
      ss << fcns[fcn]->GetParameter(i);
      for (int j = 0; j < i; j++) {
        ss << "*@0";
      }
      ss << ")";
    }
    cout << "RooFormulaVar for " << names[fcn] << " = " << ss.str() << endl;
  }


}


void interpolateAgainMassErr(int channel, int year=2012) {

  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(1.3);
  gStyle->SetMarkerColor(kGreen+3);
  gStyle->SetLineColor(kGreen+3);
  gStyle->SetLineWidth(2);
  gStyle->SetTitle("H#rightarrow ZZ Signal Lineshape Interpolation");

  stringstream filename;
  filename << "parameters_masserr_channel" << channel << ((year==2012) ? "_8TeV" : "_7TeV") << ".root";

  TFile *resultfile = TFile::Open(filename.str().c_str());
  TGraphErrors *gLdM = (TGraphErrors*)resultfile->Get("gLdM");
  TGraphErrors *gLdS = (TGraphErrors*)resultfile->Get("gLdS");
  TGraphErrors *gLnM = (TGraphErrors*)resultfile->Get("gLnM");
  TGraphErrors *gLnK = (TGraphErrors*)resultfile->Get("gLnK");
  TGraphErrors *gF1 = (TGraphErrors*)resultfile->Get("gF1");

  TCanvas *c1 = new TCanvas("c1","c1",750,750);
  c1->cd();

  stringstream channame;
  channame << (year==2011 ? "_7TeV" : "_8TeV") << "_channel" << channel << ".pdf";

  gLdM->Fit("pol2"); gLdM->Draw("Ap"); gPad->Update(); gPad->Print((string("gLdM")+channame.str()).c_str());
  gLdS->Fit("pol2"); gLdS->Draw("Ap"); gPad->Update(); gPad->Print((string("gLdS")+channame.str()).c_str());
  gLnM->Fit("pol2"); gLnM->Draw("Ap"); gPad->Update(); gPad->Print((string((channel==0) ? "gLnM" : "gGaM")+channame.str()).c_str());
  gLnK->Fit("pol2"); gLnK->Draw("Ap"); gPad->Update(); gPad->Print((string((channel==0) ? "gLnK" : "gGaS")+channame.str()).c_str());
  gF1->Fit("pol2"); gF1->Draw("Ap"); gPad->Update(); gPad->Print((string("gF1")+channame.str()).c_str());

  TF1 *fLdM = (TF1*)gLdM->GetFunction("pol2");
  TF1 *fLdS = (TF1*)gLdS->GetFunction("pol2");
  TF1 *fLnM = (TF1*)gLnM->GetFunction("pol2");
  TF1 *fLnK = (TF1*)gLnK->GetFunction("pol2");
  TF1 *fF1 = (TF1*)gF1->GetFunction("pol2");

  std::vector<string> names;
  std::vector<TF1*> fcns;
  fcns.push_back(fLdM);       names.push_back("fLdM");
  fcns.push_back(fLdS);       names.push_back("fLdS");
  fcns.push_back(fLnM);       names.push_back("fLnM");
  fcns.push_back(fLnK);       names.push_back("fLnK");
  fcns.push_back(fF1);        names.push_back("fF1");

  for(int fcn=0;fcn<(int)fcns.size();++fcn) {

    std::stringstream ss;
    
    for (int i = 0; i < fcns[fcn]->GetNumberFreeParameters(); i++) {
      if (i != 0) ss << " + (";
      else ss << "(";
      ss << fcns[fcn]->GetParameter(i);
      for (int j = 0; j < i; j++) {
        ss << "*@0";
      }
      ss << ")";
    }
    cout << "RooFormulaVar for " << names[fcn] << " = " << ss.str() << endl;
  }

}
