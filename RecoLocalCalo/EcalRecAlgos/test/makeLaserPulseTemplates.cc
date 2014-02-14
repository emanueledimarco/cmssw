#include <string>
#include <sstream>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TProfile2D.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooLandau.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooPlot.h"

using namespace std;
using namespace RooFit;

void templatesFromLaserAllCrystals(const char *dqmfile) {

  TFile *tfile = TFile::Open(dqmfile);
  string dirEB = "/DQMData/Run 202299/EcalBarrel/Run summary/EBLaserTask/Laser1";
  string dirEE = "/DQMData/Run 202299/EcalEndcap/Run summary/EELaserTask/Laser1";

  tfile->cd(dirEB.c_str());
  vector<TProfile2D> templates;
  for(int ism=1; ism<=36; ++ism) {
    stringstream ebname;
    if(ism<10) ebname << "EBLT shape EB+0" << ism << " L1";
    else if(ism<19) ebname << "EBLT shape EB+" << ism << " L1";
    else if(ism<28) ebname << "EBLT shape EB-0" << ism-18 << " L1";
    else ebname << "EBLT shape EB-" << ism-18 << " L1";
    cout << "ism = " << ism << "  SM name = " << ebname.str() << endl;
    TProfile2D *shape = (TProfile2D*)gDirectory->Get(ebname.str().c_str());
    templates.push_back(*shape);
  }

  tfile->cd(dirEE.c_str());
  for(int ism=-9; ism<=9; ++ism) {
    if(ism==0) continue;
    stringstream eename;
    if(ism<0) eename << "EELT shape EE-0" << abs(ism) << " L1";
    else eename << "EELT shape EE+0" << ism << " L1";
    cout << "ism = " << ism << "  SM name = " << eename.str() << endl;
    TProfile2D *shape = (TProfile2D*)gDirectory->Get(eename.str().c_str());
    templates.push_back(*shape);
  }

  TFile *strippedfile = TFile::Open("templates.root","recreate");
  strippedfile->cd();
  
  for(int t=0; t<(int)templates.size(); ++t) templates[t].Write();
  
  strippedfile->Close();

}

void templatesFromLaserSingleCrystal(const char *dqmfile) {

  TFile *tfile = TFile::Open(dqmfile);
  string dirEB = "/DQMData/Run 202299/EcalBarrel/Run summary/EBLaserClient";
  string dirEE = "/DQMData/Run 202299/EcalEndcap/Run summary/EELaserClient";

  tfile->cd(dirEB.c_str());
  vector<TH1F> templates;
  for(int ism=1; ism<=36; ++ism) {
    stringstream ebname;
    if(ism<10) ebname << "EBLT laser shape L1 EB+0" << ism;
    else if(ism<19) ebname << "EBLT laser shape L1 EB+" << ism;
    else if(ism<28) ebname << "EBLT laser shape L1 EB-0" << ism-18;
    else ebname << "EBLT laser shape L1 EB-" << ism-18;
    cout << "ism = " << ism << "  SM name = " << ebname.str() << endl;
    TH1F *shape = (TH1F*)gDirectory->Get(ebname.str().c_str());
    string newname = ebname.str().replace(5,5,"single");
    shape->SetName(newname.c_str());
    templates.push_back(*shape);
  }

  tfile->cd(dirEE.c_str());
  for(int ism=-9; ism<=9; ++ism) {
    if(ism==0) continue;
    stringstream eename;
    if(ism<0) eename << "EELT laser shape L1 EE-0" << abs(ism);
    else eename << "EELT laser shape L1 EE+0" << ism;
    cout << "ism = " << ism << "  SM name = " << eename.str() << endl;
    TH1F *shape = (TH1F*)gDirectory->Get(eename.str().c_str());
    string newname = eename.str().replace(5,5,"single");
    shape->SetName(newname.c_str());
    templates.push_back(*shape);
  }

  TFile *strippedfile = TFile::Open("templates.root","update");
  strippedfile->cd();
  
  for(int t=0; t<(int)templates.size(); ++t) templates[t].Write();
  
  strippedfile->Close();

}


void makePileupTemplates() {

  TFile *tfile = TFile::Open("data/templates.root");
  // take one representative template, for the moment
  TH1F *templ = (TH1F*)tfile->Get("EBLT single shape L1 EB+08");

  RooRealVar *amplitude = new RooRealVar("amplitude","time",0,10,"samples");
  RooArgList vars(*amplitude);
  RooDataHist *hist = new RooDataHist("pileupTemplate","pileupTemplate",vars,templ);

  //--- Pedestal
  RooRealVar pedSlope("pedSlope","pedSlope",7,0,100);
  RooPolynomial pedestal("pedestal","pedestal",*amplitude,pedSlope);

  //--- Landau
  RooRealVar mean("mean","mean",3,7,"samples") ;
  RooRealVar sigma("#sigma","width",1,5,"samples"); 
  RooLandau landau("landau","landau",*amplitude,mean,sigma);
  
  //--- sum of the Landau with poly
  RooRealVar frac("frac","ped. fraction",0,1);
  RooAddPdf *pdf = new RooAddPdf("pulse","pulse",landau,pedestal,frac);

  pdf->fitTo(*hist,SumW2Error(1),RooFit::Range(0.,10.),Strategy(2),NumCPU(8));

  TCanvas *canv = new TCanvas("canv","",600,600);

  RooPlot* xframe = amplitude->frame(0,10,10) ;
  xframe->SetTitle("Pulse shape");
  hist->plotOn(xframe,DataError(RooAbsData::SumW2) );
  pdf->plotOn(xframe,LineColor(kBlack),Range(0.,10.));
  pdf->plotOn(xframe, Components("pedestal"), LineStyle(kDashed), LineColor(kBlack));

  xframe->Draw(); gPad->Update(); canv->SaveAs("shapeAnalyticFit.png");


  RooAddPdf *pdfExt = new RooAddPdf("pulseExt","pulseExt",landau,pedestal,frac);
  amplitude->setRange(0,30);

  RooPlot* xframeExt = amplitude->frame(0,30,30) ;
  pdfExt->plotOn(xframeExt,LineColor(kBlack));
  xframeExt->Draw(); gPad->Update(); canv->SaveAs("shapeAnalyticFitExt.png");

}

float getNPerGeV(bool iseb) {

  float adcToGeV_EB = 0.035;
  float adcToGeV_EE = 0.060;

  TFile *tfile = TFile::Open("data/templates.root");
  // take one representative template, for the moment
  TH1F *templ = (TH1F*)tfile->Get("EBLT single shape L1 EB+08");
  
  int ped = (templ->GetBinContent(1)+templ->GetBinContent(2))/2.;
  float integral=0;
  for(int i=1;i<11;++i) {
    integral += templ->GetBinContent(i)-ped;
  }

  float ampli_1GeV = (iseb) ? 1./adcToGeV_EB : 1./adcToGeV_EE;
  float scale = ampli_1GeV/templ->GetBinContent(5);
  float frac_in25ns = 0.80;

  return integral * scale / frac_in25ns;

}

void makeTemplates(const char *dqmfile="~/Work/data/ecalreco/DQM_V0013_EcalBarrel_R000202299.root") {
  templatesFromLaserAllCrystals(dqmfile);
  templatesFromLaserSingleCrystal(dqmfile);
}
