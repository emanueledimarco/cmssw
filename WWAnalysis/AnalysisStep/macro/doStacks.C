#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
//#include "Math/VectorUtil.h"


#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

int doStacks(int cut=0,char* inputHistoName="projMet",
	     std::string folderName="./",int format=1) {
  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  using namespace std;
  //using namespace ROOT::Math::VectorUtil;

  string sampleString;
  vector<string> backgroundNames;
  vector<string> backgroundIDs;
  vector<string> signalNames;
  vector<string> signalIDs;
  
  // ----- here I set the histogram to be read ------
  stringstream histoNameSS;
  histoNameSS << inputHistoName << cut;
  string histoToRead(histoNameSS.str());

  // ----- here I read the scalers from ascii files -----
  std::ifstream inFileSig("scalers.sig.txt");
  if(!inFileSig) {cout << endl << "Failed to open signal scalers file " << endl; return 1;}
  std::ifstream inFileBkg("scalers.bkg.txt");
  if(!inFileBkg) {cout << endl << "Failed to open bkg scalers file " << endl; return 1;}


  double scaler;
  int counter(1);

  // ----- here I set the signal dataset names and scalers ------
  vector<double> sigScalers;
  while(!inFileSig.eof()){
    char line[256];
    stringstream ss;
    string sampleString;
    string sampleID;
    inFileSig.getline(line,256,'\n');
    if(line[0]=='#' || line[0]=='\0') continue;
    ss << line;
    ss >> sampleID >> sampleString >> scaler; 
    //cout << "counter, scaler: " << counter << " , " << scaler << endl;
    signalNames.push_back(sampleString);
    signalIDs.push_back(sampleID);
    sigScalers.push_back(scaler);
    counter++;
  }


  // ----- here I set the background dataset names and scalers------
  counter=1;
  vector<double> bkgScalers;
  while(!inFileBkg.eof()){
    char line[256];
    stringstream ss;
    string sampleString;
    string sampleID;
    inFileBkg.getline(line,256,'\n');
    if(line[0]=='#' || line[0]=='\0') continue;
    ss << line;
    ss >> sampleID >>sampleString >> scaler; 
    //cout << "counter, scaler: " << counter << " , " << scaler << endl;
    backgroundNames.push_back(sampleString);
    backgroundIDs.push_back(sampleID);
    bkgScalers.push_back(scaler);
    counter++;
  }



  /* ------------------------------------------
     THE REAL ACTION STARTS HERE     
     -------------------------------------------- */
  //cout.precision(3);
  TCanvas* myC = new TCanvas("myCanvas","myCanvas",0,0,600,600); myC->cd(1);

  double entries;
  TFile* f;
  string fileName;

  //fileName=folderName+"histo."+signalNames[0]+".root";
  fileName=folderName+"histo."+signalIDs[0]+".root";
  f = TFile::Open(fileName.c_str());
  TH1F* tmpHisto = (TH1F*) f->Get(histoToRead.c_str());
  int nbins  = tmpHisto->GetXaxis()->GetNbins();
  double min = tmpHisto->GetXaxis()->GetXmin();
  double max = tmpHisto->GetXaxis()->GetXmax();
  string xlabel = tmpHisto->GetTitle();
  f->Close();

  TH1F* sigSumH = new TH1F("sigSumH","sum of signals",nbins,min,max);   
  TH1F* bkgSumH = new TH1F("bkgSumH","sum of backgrounds",nbins,min,max);   
  double sumErrorSig(0.);
  double sumErrorBkg(0.);

  cout << endl << "==== Counting of events for the signal =====" << endl;
  for(unsigned int i=0; i<signalIDs.size();++i){
    fileName=folderName+"histo."+signalIDs[i]+".root";
    cout << "filename to be opened: " << fileName << endl;
    f = TFile::Open(fileName.c_str());
    TH1F* histo = (TH1F*) f->Get(histoToRead.c_str());

    sigSumH->Add(histo,sigScalers[i]);
    entries = histo->GetEntries(); 
    double uncertainty = sqrt(1.*entries);
    sumErrorSig+=uncertainty*uncertainty*sigScalers[i]*sigScalers[i];

    if(entries){
      if(format==0)
	cout << signalNames[i] << ": " << setprecision(0) << fixed << entries << endl;
      else{
	cout << setw(15) << left << signalNames[i] << ": " 
	     << setprecision(3) << setw(8) << right << fixed << entries*sigScalers[i] 
	     << "+- " 
	     << setprecision(3) << setw(8) << fixed << uncertainty*sigScalers[i] << endl;
      }
      
    }
    f->Close();
  }

  if(bkgScalers.size() == backgroundNames.size()){
    //cout << "size of scalers and list of bkg samples are consistent" << endl;
  }else{ cout << "ERROR: bkgScalers.size() != backgroundNames" << endl; return 1;}

  cout << endl << "==== Counting events for the bkg =====" << endl;
  for(unsigned int i=0; i<backgroundIDs.size();++i){
    fileName=folderName+"histo."+backgroundIDs[i]+".root";
    f = TFile::Open(fileName.c_str());
    TH1F* histo = (TH1F*) f->Get(histoToRead.c_str());   

    //histo->Rebin(2);
    bkgSumH->Add(histo,bkgScalers[i]);	
    entries = histo->GetEntries(); 

    double uncertainty = sqrt(1.*entries);
    sumErrorBkg+=uncertainty*uncertainty*bkgScalers[i]*bkgScalers[i];

    if(entries){
      if(format==0)
	cout << backgroundNames[i] << ": " << setprecision(0) << fixed << entries << endl;
      else      
	cout << setw(15) << left << backgroundNames[i] << ": " 
	     << setprecision(3) << setw(8) << right << fixed << entries*bkgScalers[i] 
	     << " +- " 
	     << setprecision(3) << setw(8) << fixed << uncertainty*bkgScalers[i] << endl;
      
    }
    f->Close();
  }
  
  //bkgSumH->SetLineColor(1);
  //sigSumH->SetLineColor(2);

  int sigColor(2);
  int bkgColor(31);

  sigSumH->SetFillColor(sigColor);
  bkgSumH->SetFillColor(bkgColor);  



  //bkgSumH->Draw();
  //sigSumH->Draw("same");

  stringstream hsTitle;
  hsTitle << "signal+background" << " (step " << cut << ")";
  //THStack* hstack = new THStack("hstack",hsTitle.str().c_str());
  THStack* hstack = new THStack("hstack","");

  double integralSig = sigSumH->Integral(0,nbins+1);
  double integralBkg = bkgSumH->Integral(0,nbins+1);


  cout << endl << "============ final sum: ==============" << endl;
  cout <<  "tot sig: " << setprecision(3) << setw(15) << left << integralSig 
       << " +- " << setprecision(3) << setw(8) << right << fixed << sqrt(sumErrorSig) << endl;
  cout << "tot bkg: " << setprecision(3) << setw(15) << left << integralBkg 
       << " +- " << setprecision(3) << setw(8) << right << fixed << sqrt(sumErrorBkg) << endl;
  hstack->Add(bkgSumH);
  hstack->Add(sigSumH);
  hstack->Draw();
  hstack->GetXaxis()->SetTitle(xlabel.c_str());



  stringstream signalLabel; signalLabel << "signal yield: " << integralSig;
  stringstream bkgLabel; bkgLabel << "bkg yield: " << integralBkg;

  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(sigSumH,signalLabel.str().c_str(),"F");
  leg->AddEntry(bkgSumH,bkgLabel.str().c_str(),"F");
  /// leg->Draw();

  /// myC->Print("stack.png");
  cout << "======================================" << endl;
  return 0;
}




