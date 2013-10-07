// ----- here are the ROOT includes -----
#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"

// ----- here are the rooFit includes -----
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
#include "Math/MinimizerOptions.h"

// ----- here is the dependency from CMSSW package -----
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"


// ----- standar library includes -----
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

/*
  This macro fits the shape of signal ggH MC for all the mass bins available. Then
  fits the parameters of the signal model as a function of mH.
  Plots for both the fit in each mass bins and for the signal shape parameters 
  are produced.

  Function fitSignalShapeW() is used to fit a single mass bin.
  Function all() iterates over all the available mass bins, calling
  fitSignalShapeW() for each of them. 

  When calling all(channel,year), user should provide year (2011 or 2012) 
  to define the set of input files and channel. For the latter the convention is:

    channels = 0 --> 4mu
    channels = 1 --> 4el
    channels = 2 --> 2mu2el


  Prerequisites to run the code:
  - compiled version of HiggsAnalysis/CombinedLimit containing the relativistic-BW 
  function definition. The corresponding library has to be loaded in root.
 
  - input step3 trees containing updated puW and lepW weights.
  

 */




using namespace RooFit;

int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}


void fitSignalShapeW(int massBin=120, int id=1120, int ch=0, int year = 2011,
		     float lumi=10,double rangeLow=105., double rangeHigh=130.,
		     double bwSigma=1.,
		     double fitValues[5]=0, double fitErrors[5]=0);

void all(int channels=0,int year=2011,double lumi=1.){

  double bwSigma[30];
  int mass[30]; int id[30]; double xLow[30]; double xHigh[30];  
  int maxMassBin;
  if(year==2011){
    mass[0] = 115; id[0]=200; xLow[0] = 100; xHigh[0] = 125; bwSigma[0] = 3.1/1000.;
    mass[1] = 120; id[1]=201; xLow[1] = 105; xHigh[1] = 130; bwSigma[1] = 3.5/1000.;
    mass[2] = 130; id[2]=202; xLow[2] = 115; xHigh[2] = 140; bwSigma[2] = 4.9/1000.;
    mass[3] = 140; id[3]=203; xLow[3] = 125; xHigh[3] = 150; bwSigma[3] = 8.1/1000.;
    mass[4] = 150; id[4]=204; xLow[4] = 135; xHigh[4] = 160; bwSigma[4] = 1.7/100.;
    mass[5] = 160; id[5]=205; xLow[5] = 145; xHigh[5] = 170; bwSigma[5] = 8.3/100.;
    mass[6] = 170; id[6]=206; xLow[6] = 155; xHigh[6] = 180; bwSigma[6] = 3.8/10.;
    mass[7] = 180; id[7]=207; xLow[7] = 165; xHigh[7] = 190; bwSigma[7] = 6.3/10.;
    mass[8] = 190; id[8]=208; xLow[8] = 170; xHigh[8] = 200; bwSigma[8] = 1.04;
    mass[9] = 200; id[9]=209; xLow[9] = 175; xHigh[9] = 210; bwSigma[9] = 1.43;
    mass[10] = 210; id[10]=210; xLow[10] = 185; xHigh[10] = 220; bwSigma[10] = 1.85;
    mass[11] = 220; id[11]=211; xLow[11] = 190; xHigh[11] = 230; bwSigma[11] = 2.31;
    maxMassBin = 12;
  }


  if(year==2012){
    mass[0] = 115; id[0]=1115; xLow[0] = 100; xHigh[0] = 125; bwSigma[0] = 3.1/1000.;
    mass[1] = 117; id[1]=1117; xLow[1] = 102; xHigh[1] = 127; bwSigma[1] = 3.2/1000.;
    mass[2] = 119; id[2]=1119; xLow[2] = 104; xHigh[2] = 129; bwSigma[2] = 3.4/1000.;
    mass[3] = 120; id[3]=1120; xLow[3] = 105; xHigh[3] = 130; bwSigma[3] = 3.5/1000.;
    mass[4] = 121; id[4]=1121; xLow[4] = 106; xHigh[4] = 131; bwSigma[4] = 3.6/1000.;
    mass[5] = 123; id[5]=1123; xLow[5] = 108; xHigh[5] = 133; bwSigma[5] = 3.8/1000.;
    mass[6] = 124; id[6]=1124; xLow[6] = 109; xHigh[6] = 134; bwSigma[6] = 3.9/1000.;
    mass[7] = 125; id[7]=1125; xLow[7] = 110; xHigh[7] = 135; bwSigma[7] = 4.0/1000.;
    mass[8] = 126; id[8]=1126; xLow[8] = 111; xHigh[8] = 136; bwSigma[8] = 4.2/1000.;
    mass[9] = 127; id[9]=1127; xLow[9] = 112; xHigh[9] = 137; bwSigma[9] = 4.3/1000.;
    mass[10] = 145; id[10]=1145; xLow[10] = 130; xHigh[10] = 155; bwSigma[10] = 1.1/100.;
    //
    mass[11] = 150; id[11]=1150; xLow[11] = 135; xHigh[11] = 160; bwSigma[11] = 1.7/100.;
    mass[12] = 180; id[12]=1180; xLow[12] = 150; xHigh[12] = 195; bwSigma[12] = 6.3/10.;
    mass[13] = 200; id[13]=1200; xLow[13] = 165; xHigh[13] = 220; bwSigma[13] = 1.4;
    mass[14] = 220; id[14]=1220; xLow[14] = 175; xHigh[14] = 250; bwSigma[14] = 2.3;
    maxMassBin = 15;
  }
  // -----------------------



  double massV[30],massE[30];
  for(int i=0; i<maxMassBin;++i){
    massV[i]=mass[i];
    massE[i]=0;
  }

  double aVal[30],aErr[30];
  double nVal[30],nErr[30];
  double meanCBVal[30],meanCBErr[30];
  double sigmaCBVal[30],sigmaCBErr[30];
  double meanBWVal[30],meanBWErr[30];

  double fitValues[5];
  double fitErrors[5];

  double extendL(1),extendH(1);

  if(channels==0) {extendL=1.0;extendH=1.0;}
  if(channels==1) {extendL=0.90;extendH=1.05;}
  if(channels==2) {extendL=0.95;extendH=1.04;}



  //fitSignalShapeW(mass[0],channels,10.,doSfLepton,max(95.,xLow[0]*0.90),xHigh[0]*1.1,bwSigma[12],
  // 		  fitValues,fitErrors);  Wait();
  
  for(int i=0; i<maxMassBin;++i){

    fitSignalShapeW(mass[i],id[i],channels,year,lumi,max(95.,xLow[i]*extendL),xHigh[i]*extendH,bwSigma[i],
		    fitValues,fitErrors);  
  
    cout << "a value,error: " << fitValues[0] << " , " << fitErrors[0] << endl; 
    aVal[i]=fitValues[0]; aErr[i]=fitErrors[0];

    cout << "n value,error: " << fitValues[3] << " , " << fitErrors[3] << endl; 
    nVal[i]=fitValues[3]; nErr[i]=fitErrors[3];

    cout << "meanCB value,error: " << fitValues[1] << " , " << fitErrors[1] << endl;
    meanCBVal[i]=fitValues[1]; meanCBErr[i]=fitErrors[1];
    
    cout << "sigmaCB value,error: " << fitValues[4] << " , " << fitErrors[4] << endl; 
    sigmaCBVal[i]=fitValues[4]; sigmaCBErr[i]=fitErrors[4];

    cout << "meanBW value,error: " << fitValues[2] << " , " << fitErrors[2] << endl; 
    meanBWVal[i]=fitValues[2]; meanBWErr[i]=fitErrors[2];

    //Wait();
  }
  

  TGraphErrors* gA = new TGraphErrors(maxMassBin,massV,aVal,massE,aErr);
  TGraphErrors* gN = new TGraphErrors(maxMassBin,massV,nVal,massE,nErr);
  TGraphErrors* gMeanCB = new TGraphErrors(maxMassBin,massV,meanCBVal,massE,meanCBErr);
  TGraphErrors* gSigmaCB = new TGraphErrors(maxMassBin,massV,sigmaCBVal,massE,sigmaCBErr);
  TGraphErrors* gMeanBW = new TGraphErrors(maxMassBin,massV,meanBWVal,massE,meanBWErr);

  gA->SetMarkerStyle(20);   gA->SetMarkerSize(1);
  gN->SetMarkerStyle(20);   gN->SetMarkerSize(1);
  gMeanCB->SetMarkerStyle(20);   gMeanCB->SetMarkerSize(1);
  gSigmaCB->SetMarkerStyle(20);   gSigmaCB->SetMarkerSize(1);
  gMeanBW->SetMarkerStyle(20);   gMeanBW->SetMarkerSize(1);
  

  gA->SetTitle("");
  gA->GetXaxis()->SetTitle("mass (GeV)");
  gA->GetYaxis()->SetTitle("CB a-parameter");

  gN->SetTitle("");
  gN->GetXaxis()->SetTitle("mass (GeV)");
  gN->GetYaxis()->SetTitle("CB n-parameter");

  gMeanCB->SetTitle("");
  gMeanCB->GetXaxis()->SetTitle("mass (GeV)");
  gMeanCB->GetYaxis()->SetTitle("CB mean");

  gSigmaCB->SetTitle("");
  gSigmaCB->GetXaxis()->SetTitle("mass (GeV)");
  gSigmaCB->GetYaxis()->SetTitle("CB sigma");

  gMeanBW->SetTitle("");
  gMeanBW->GetXaxis()->SetTitle("mass (GeV)");
  gMeanBW->GetYaxis()->SetTitle("BW mean");



  gA->Fit("pol0"); gA->Draw("Ap"); gPad->Update(); gPad->Print("aFit.pdf"); Wait();
  gN->Fit("pol1"); gN->Draw("Ap"); gPad->Update(); gPad->Print("nFit.pdf"); Wait();
  gMeanCB->Fit("pol1"); gMeanCB->Draw("Ap"); gPad->Update(); gPad->Print("meanFit.pdf"); Wait();
  gSigmaCB->Fit("pol1"); gSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print("sigmaFit.pdf"); Wait();


}



void fitSignalShapeW(int massBin,int id, int channels, int year, 
		     float lumi,double rangeLow, double rangeHigh,
		     double bwSigma,
		     double fitValues[5], double fitErrors[5]){
 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  stringstream ggFileName,vbfFileName;
  ggFileName << "/data/hzz4l/step3/" << year << "/MC/HZZ4L_S3_V02/hzzTree_id" << id << ".root"; 


  TFile* ggFile = new TFile(ggFileName.str().c_str()); 
  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  float channel;
  float nTrueInt,mass;
  double lepW(1.);
  double puW(1.);

  int  nentries = ggTree->GetEntries();
 
 
  //--- ggTree part
  ggTree->SetBranchAddress("channel",&channel);
  ggTree->SetBranchAddress("numTrueInteractions",&nTrueInt);
  ggTree->SetBranchAddress("mass",&mass);
  ggTree->SetBranchAddress("puW",&puW);
  ggTree->SetBranchAddress("lepW",&lepW);


  //--- rooFit part
  double xMin,xMax,xInit;
  xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;
  //xMin = 114.;
  //xMax = 122 ;

  //---------  
  TCut cut = "channel == 0";
  RooRealVar x("mass","mass (GeV)",xInit,xMin,xMax);
  //RooRealVar x("m4l","m4l",xInit,xMin,xMax);
  RooRealVar w("myW","myW",1.0,0.,1000.);
  RooArgSet ntupleVarSet(x,w);
  RooDataSet dataset("mass4l","mass4l",ntupleVarSet,WeightVar("myW"));


  for(int k=0; k<nentries; k++){
    ggTree->GetEvent(k);

    if(channels==0){if(channel != 0) continue;}
    if(channels==1){if(channel != 1) continue;}
    if(channels==2){if(channel != 2 && channel != 3) continue;}

    ntupleVarSet.setRealValue("mass",mass);
    double localW = puW*lepW*lumi;

    ntupleVarSet.setRealValue("myW",localW);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, localW);
    
  }
  //---------

  //cout << "dataset n entries: " << dataset.sumEntries() << endl;
  //cout << "reduced data n entries: " << reducedData->sumEntries() << endl;


  //--- simple CrystalBall
  RooRealVar mean("mean","mean of gaussian",0,-5.0,5.0) ;
  RooRealVar sigma("sigma","width of gaussian",1.5,0.,5.); 
  RooRealVar a("a","a",1.46,0.,4.);
  RooRealVar n("n","n",1.92,0.,25.);   
  RooCBShape CBall("CBall","Crystal ball",x, mean,sigma, a,n);
  

  //--- simple Gaussian
  RooRealVar mean2("mean2","mean2 of gaussian",xInit,xInit*0.80,xInit*1.2) ;
  RooRealVar sigma2("sigma2","width2 of gaussian",10.,5.,500.); 
  RooGaussian tailCatcher("tailCatcher","tailCatcher",x,mean2,sigma2);
  RooRealVar fsig("fsig","signal fraction",0.95,0.7,1.);


  //--- Breit-Wigner
  RooRealVar mean3("mean3","mean3",xInit) ;
  RooRealVar sigma3("sigma3","width3",bwSigma); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooRelBWUFParam bw("bw","bw",x,mean3,scale3);

  //RooAddPdf model("model","model",RooArgList(CBall,tailCatcher),fsig);
  //RooCBShape model("model","model",x, mean,sigma, a,n);
  x.setBins(10000,"fft");
  RooFFTConvPdf model("model","model",x,bw,CBall);

  

  model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channels==0){frameTitle << "4mu, mH = ";}
  if(channels==1){frameTitle << "4el, mH = ";}
  if(channels==2){frameTitle << "2mu2el, mH = ";}
  frameTitle << massBin << " GeV";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  model.plotOn(xframe);
  //model.paramOn(xframe);


  stringstream nameFile;
  nameFile << "fitM" << massBin << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = a.getVal();
    fitValues[1] = mean.getVal();
    fitValues[2] = mean3.getVal();
    fitValues[3] = n.getVal();
    fitValues[4] = sigma.getVal();
  }

  if(fitErrors!=0){
    fitErrors[0] = a.getError();
    fitErrors[1] = mean.getError();
    fitErrors[2] = mean3.getError();
    fitErrors[3] = n.getError();
    fitErrors[4] = sigma.getError();
  }

}



