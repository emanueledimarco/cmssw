#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
// #include <RooRealVar.h>
// #include <RooArgSet.h>
// #include <RooDataSet.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TTree.h>
#include <TFile.h>

// #include "FitMaker.h"
// #include "../FakeRateCalculator.h"
#include "../scales2.h"

// using namespace RooFit;


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


void MakePseudoscalarPlots(string filenameSMHiggs, string filenamePSHiggs, vector<pair<string,double> > ZZBkgfilenames,
                           vector<string> filenamesData
  ) {

  //*******************************************************************
  // Setup
  //*******************************************************************


  //*******************************************************************
  // Histograms
  //*******************************************************************
  TH1F *histMass_SMHiggs = new TH1F( "histMass_SMHiggs", " ; m_{4l} [GeV/c^{2} ; Number of Events", 20, 110, 150);
  TH1F *histMass_ZZBkg = new TH1F( "histMass_ZZBkg", " ; m_{4l} [GeV/c^{2} ; Number of Events", 20, 110, 150);
  TH1F *histMass_FakeBkg = new TH1F( "histMass_FakeBkg", " ; m_{4l} [GeV/c^{2} ; Number of Events", 20, 110, 150);
  TH1F *histMass_Data = new TH1F( "histMass_Data", " ; m_{4l} [GeV/c^{2} ; Number of Events", 20, 110, 150);


  TH1F *histPSMELA_SMHiggs = new TH1F( "histPSMELA_SMHiggs", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);
  TH1F *histPSMELA_AfterMELACut_SMHiggs = new TH1F( "histPSMELA_AfterMELACut_SMHiggs", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);
  TH1F *histPSMELA_PSHiggs = new TH1F( "histPSMELA_PSHiggs", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);
  TH1F *histPSMELA_AfterMELACut_PSHiggs = new TH1F( "histPSMELA_AfterMELACut_PSHiggs", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);
  TH1F *histPSMELA_ZZBkg = new TH1F( "histPSMELA_ZZBkg", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);
  TH1F *histPSMELA_AfterMELACut_ZZBkg = new TH1F( "histPSMELA_AfterMELACut_ZZBkg", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);
  TH1F *histPSMELA_Data = new TH1F( "histPSMELA_Data", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);
  TH1F *histPSMELA_AfterMELACut_Data = new TH1F( "histPSMELA_AfterMELACut_Data", " ; Pseudoscalar MELA ; Number of Events", 40, 0, 1);


  TH1F *histCosThetaStar_SMHiggs = new TH1F( "histCosThetaStar_SMHiggs", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_SMHiggs = new TH1F( "histCosTheta1_SMHiggs", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_SMHiggs = new TH1F( "histCosTheta2_SMHiggs", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_SMHiggs = new TH1F( "histPhi_SMHiggs", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_SMHiggs = new TH1F( "histPhiStar1_SMHiggs", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_SMHiggs = new TH1F( "histZ1Mass_SMHiggs", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_SMHiggs = new TH1F( "histZ2Mass_SMHiggs", " ; m_{Z2} ; Number of Events", 35, 0, 70);
  TH1F *histCosThetaStar_AfterMELACut_SMHiggs = new TH1F( "histCosThetaStar_AfterMELACut_SMHiggs", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_AfterMELACut_SMHiggs = new TH1F( "histCosTheta1_AfterMELACut_SMHiggs", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_AfterMELACut_SMHiggs = new TH1F( "histCosTheta2_AfterMELACut_SMHiggs", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_AfterMELACut_SMHiggs = new TH1F( "histPhi_AfterMELACut_SMHiggs", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_AfterMELACut_SMHiggs = new TH1F( "histPhiStar1_AfterMELACut_SMHiggs", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_AfterMELACut_SMHiggs = new TH1F( "histZ1Mass_AfterMELACut_SMHiggs", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_AfterMELACut_SMHiggs = new TH1F( "histZ2Mass_AfterMELACut_SMHiggs", " ; m_{Z2} ; Number of Events", 35, 0, 70);

  TH1F *histCosThetaStar_PSHiggs = new TH1F( "histCosThetaStar_PSHiggs", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_PSHiggs = new TH1F( "histCosTheta1_PSHiggs", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_PSHiggs = new TH1F( "histCosTheta2_PSHiggs", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_PSHiggs = new TH1F( "histPhi_PSHiggs", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_PSHiggs = new TH1F( "histPhiStar1_PSHiggs", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_PSHiggs = new TH1F( "histZ1Mass_PSHiggs", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_PSHiggs = new TH1F( "histZ2Mass_PSHiggs", " ; m_{Z2} ; Number of Events", 35, 0, 70);
  TH1F *histCosThetaStar_AfterMELACut_PSHiggs = new TH1F( "histCosThetaStar_AfterMELACut_PSHiggs", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_AfterMELACut_PSHiggs = new TH1F( "histCosTheta1_AfterMELACut_PSHiggs", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_AfterMELACut_PSHiggs = new TH1F( "histCosTheta2_AfterMELACut_PSHiggs", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_AfterMELACut_PSHiggs = new TH1F( "histPhi_AfterMELACut_PSHiggs", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_AfterMELACut_PSHiggs = new TH1F( "histPhiStar1_AfterMELACut_PSHiggs", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_AfterMELACut_PSHiggs = new TH1F( "histZ1Mass_AfterMELACut_PSHiggs", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_AfterMELACut_PSHiggs = new TH1F( "histZ2Mass_AfterMELACut_PSHiggs", " ; m_{Z2} ; Number of Events", 35, 0, 70);

  TH1F *histCosThetaStar_ZZBkg = new TH1F( "histCosThetaStar_ZZBkg", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_ZZBkg = new TH1F( "histCosTheta1_ZZBkg", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_ZZBkg = new TH1F( "histCosTheta2_ZZBkg", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_ZZBkg = new TH1F( "histPhi_ZZBkg", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_ZZBkg = new TH1F( "histPhiStar1_ZZBkg", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_ZZBkg = new TH1F( "histZ1Mass_ZZBkg", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_ZZBkg = new TH1F( "histZ2Mass_ZZBkg", " ; m_{Z2} ; Number of Events", 35, 0, 70);
  TH1F *histCosThetaStar_AfterMELACut_ZZBkg = new TH1F( "histCosThetaStar_AfterMELACut_ZZBkg", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_AfterMELACut_ZZBkg = new TH1F( "histCosTheta1_AfterMELACut_ZZBkg", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_AfterMELACut_ZZBkg = new TH1F( "histCosTheta2_AfterMELACut_ZZBkg", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_AfterMELACut_ZZBkg = new TH1F( "histPhi_AfterMELACut_ZZBkg", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_AfterMELACut_ZZBkg = new TH1F( "histPhiStar1_AfterMELACut_ZZBkg", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_AfterMELACut_ZZBkg = new TH1F( "histZ1Mass_AfterMELACut_ZZBkg", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_AfterMELACut_ZZBkg = new TH1F( "histZ2Mass_AfterMELACut_ZZBkg", " ; m_{Z2} ; Number of Events", 35, 0, 70);

  TH1F *histCosThetaStar_Data = new TH1F( "histCosThetaStar_Data", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_Data = new TH1F( "histCosTheta1_Data", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_Data = new TH1F( "histCosTheta2_Data", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_Data = new TH1F( "histPhi_Data", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_Data = new TH1F( "histPhiStar1_Data", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_Data = new TH1F( "histZ1Mass_Data", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_Data = new TH1F( "histZ2Mass_Data", " ; m_{Z2} ; Number of Events", 35, 0, 70);
  TH1F *histCosThetaStar_AfterMELACut_Data = new TH1F( "histCosThetaStar_AfterMELACut_Data", " ; CosThetaStar ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta1_AfterMELACut_Data = new TH1F( "histCosTheta1_AfterMELACut_Data", " ; CosTheta1 ; Number of Events", 40, -1, 1);
  TH1F *histCosTheta2_AfterMELACut_Data = new TH1F( "histCosTheta2_AfterMELACut_Data", " ; CosTheta2 ; Number of Events", 40, -1, 1);
  TH1F *histPhi_AfterMELACut_Data = new TH1F( "histPhi_AfterMELACut_Data", " ; Phi ; Number of Events", 40, -3.142, 3.142);
  TH1F *histPhiStar1_AfterMELACut_Data = new TH1F( "histPhiStar1_AfterMELACut_Data", " ; PhiStar1 ; Number of Events", 40, -3.142, 3.142);
  TH1F *histZ1Mass_AfterMELACut_Data = new TH1F( "histZ1Mass_AfterMELACut_Data", " ; m_{Z1} ; Number of Events", 40, 40, 120);
  TH1F *histZ2Mass_AfterMELACut_Data = new TH1F( "histZ2Mass_AfterMELACut_Data", " ; m_{Z2} ; Number of Events", 35, 0, 70);



  //*******************************************************************
  // Compute Signal Yields
  //*******************************************************************
  

  //*******************************************************************
  // Open Input Files
  //*******************************************************************
  TFile* fileSMHiggs = new TFile(filenameSMHiggs.c_str());
  TTree* treeSMHiggs = (TTree*)fileSMHiggs->Get("zz4lTree/probe_tree");
  TFile* filePSHiggs = new TFile(filenamePSHiggs.c_str());
  TTree* treePSHiggs = (TTree*)filePSHiggs->Get("zz4lTree/probe_tree");
            

  float channel   = 0.0;
  float z1mass    = 0.0;
  float z2mass    = 0.0;
  float mass      = 0.0;
  float numsim    = 0.0;
  float l3pt      = 0.0;
  float l3eta     = 0.0;
  float l3id      = 0.0;
  float l3iso     = 0.0;
  float l3pdgId   = 0.0;
  float l4pt      = 0.0;
  float l4eta     = 0.0;
  float l4id      = 0.0;
  float l4iso     = 0.0;
  float l4pdgId   = 0.0;
  float l1pt      = 0.0;
  float l1eta     = 0.0;
  float l1pdgId   = 0.0;
  float l2pt      = 0.0;
  float l2eta     = 0.0;
  float l2pdgId   = 0.0;
  float mela               = 0.0;
  float melaPS             = 0.0;
  float melaSpinTwoMinimal = 0.0;
  float melaPtY            = 0.0;
  float melaCosThetaStar   = 0.0; 
  float melaCosTheta1      = 0.0;
  float melaCosTheta2      = 0.0;
  float melaPhi            = 0.0;
  float melaPhiStar1       = 0.0;
  float higgswgt  = 0.0;
  unsigned event  = 0;
  unsigned run    = 0;
  unsigned lumi   = 0;


  //*******************************************************************
  // Loop Over SM Higgs
  //*******************************************************************
  treeSMHiggs->SetBranchAddress( "channel", &channel);
  treeSMHiggs->SetBranchAddress( "z1mass",&z1mass);
  treeSMHiggs->SetBranchAddress( "z2mass",&z2mass);
  treeSMHiggs->SetBranchAddress("mass" ,&mass);
  treeSMHiggs->SetBranchAddress( "numTrueInteractions",&numsim);
  treeSMHiggs->SetBranchAddress("l3pt" ,&l3pt);
  treeSMHiggs->SetBranchAddress("l3eta" ,&l3eta);
  treeSMHiggs->SetBranchAddress("l3idNew" ,&l3id);
  treeSMHiggs->SetBranchAddress("l3pfIsoComb04EACorr" ,&l3iso);
  treeSMHiggs->SetBranchAddress("l3pdgId" ,&l3pdgId);
  treeSMHiggs->SetBranchAddress( "l4pt",&l4pt);
  treeSMHiggs->SetBranchAddress("l4eta" ,&l4eta);
  treeSMHiggs->SetBranchAddress("l4idNew" ,&l4id);
  treeSMHiggs->SetBranchAddress( "l4pfIsoComb04EACorr",&l4iso);
  treeSMHiggs->SetBranchAddress("l4pdgId" ,&l4pdgId);
  treeSMHiggs->SetBranchAddress( "l1pt",&l1pt);
  treeSMHiggs->SetBranchAddress( "l1eta",&l1eta);
  treeSMHiggs->SetBranchAddress( "l1pdgId",&l1pdgId);
  treeSMHiggs->SetBranchAddress( "l2pt",&l2pt);
  treeSMHiggs->SetBranchAddress( "l2eta",&l2eta);
  treeSMHiggs->SetBranchAddress( "l2pdgId",&l2pdgId);
  treeSMHiggs->SetBranchAddress( "melaLD",&mela);
  treeSMHiggs->SetBranchAddress( "melaPSLD",&melaPS);
  treeSMHiggs->SetBranchAddress( "melaSpinTwoMinimal",&melaSpinTwoMinimal);
  treeSMHiggs->SetBranchAddress( "melaPtY",&melaPtY);
  treeSMHiggs->SetBranchAddress( "melaCosThetaStar",&melaCosThetaStar);
  treeSMHiggs->SetBranchAddress( "melaCosTheta1",&melaCosTheta1);
  treeSMHiggs->SetBranchAddress( "melaCosTheta2",&melaCosTheta2);
  treeSMHiggs->SetBranchAddress( "melaPhi",&melaPhi);
  treeSMHiggs->SetBranchAddress( "melaPhiStar1",&melaPhiStar1);
  treeSMHiggs->SetBranchAddress( "event" ,&event);
  treeSMHiggs->SetBranchAddress( "run",&run);
  treeSMHiggs->SetBranchAddress( "lumi",&lumi);
  treeSMHiggs->SetBranchAddress( "genhiggsmassweight",&higgswgt);
        
  for (int i = 0; i < treeSMHiggs->GetEntries(); i++) {
    treeSMHiggs->GetEntry(i);
    bool existsAlready = false;
       
    if (existsAlready) {
      std::cout << "Run : " << run << " Lumi : " << lumi << " Event : " << event << std::endl;
    } 
 
    if (!existsAlready) {


      float weight = getBkgXsec(1125)/evt_8TeV(1125);
      float weighterr = 1.0;

      weight  *= getPUWeight((int)numsim, 1);
      weight *= getSF(l1pt, l1eta, l1pdgId);
      weight *= getSF(l2pt, l2eta, l2pdgId);
      weight *= getSF(l3pt, l3eta, l3pdgId);
      weight *= getSF(l4pt, l4eta, l4pdgId);

      weighterr *= getSF(l1pt, l1eta, l1pdgId);
      weighterr *= getSF(l2pt, l2eta, l2pdgId);
      weighterr *= getSF(l3pt, l3eta, l3pdgId);
      weighterr *= getSF(l4pt, l4eta, l4pdgId);

      weight    *= higgswgt;
      weighterr *= higgswgt;

      //make plot within 10GeV mass window around the higgs
      if (mass > 120 && mass < 130) {
        histPSMELA_SMHiggs->Fill(melaPS, weight);
        histCosThetaStar_SMHiggs->Fill(melaCosThetaStar, weight);
        histCosTheta1_SMHiggs->Fill(melaCosTheta1, weight);
        histCosTheta2_SMHiggs->Fill(melaCosTheta2, weight);
        histPhi_SMHiggs->Fill(melaPhi, weight);
        histPhiStar1_SMHiggs->Fill(melaPhiStar1, weight);
        histZ1Mass_SMHiggs->Fill(z1mass, weight);
        histZ2Mass_SMHiggs->Fill(z2mass, weight);
        if (mela > 0.5) {
          histPSMELA_AfterMELACut_SMHiggs->Fill(melaPS,weight);
          histCosThetaStar_AfterMELACut_SMHiggs->Fill(melaCosThetaStar, weight);
          histCosTheta1_AfterMELACut_SMHiggs->Fill(melaCosTheta1, weight);
          histCosTheta2_AfterMELACut_SMHiggs->Fill(melaCosTheta2, weight);
          histPhi_AfterMELACut_SMHiggs->Fill(melaPhi, weight);
          histPhiStar1_AfterMELACut_SMHiggs->Fill(melaPhiStar1, weight);
          histZ1Mass_AfterMELACut_SMHiggs->Fill(z1mass, weight);
          histZ2Mass_AfterMELACut_SMHiggs->Fill(z2mass, weight);
        }
      }

    }
  }


  //*******************************************************************
  // Loop Over SM Higgs
  //*******************************************************************
  treePSHiggs->SetBranchAddress( "channel", &channel);
  treePSHiggs->SetBranchAddress( "z1mass",&z1mass);
  treePSHiggs->SetBranchAddress( "z2mass",&z2mass);
  treePSHiggs->SetBranchAddress("mass" ,&mass);
  treePSHiggs->SetBranchAddress( "numTrueInteractions",&numsim);
  treePSHiggs->SetBranchAddress("l3pt" ,&l3pt);
  treePSHiggs->SetBranchAddress("l3eta" ,&l3eta);
  treePSHiggs->SetBranchAddress("l3idNew" ,&l3id);
  treePSHiggs->SetBranchAddress("l3pfIsoComb04EACorr" ,&l3iso);
  treePSHiggs->SetBranchAddress("l3pdgId" ,&l3pdgId);
  treePSHiggs->SetBranchAddress( "l4pt",&l4pt);
  treePSHiggs->SetBranchAddress("l4eta" ,&l4eta);
  treePSHiggs->SetBranchAddress("l4idNew" ,&l4id);
  treePSHiggs->SetBranchAddress( "l4pfIsoComb04EACorr",&l4iso);
  treePSHiggs->SetBranchAddress("l4pdgId" ,&l4pdgId);
  treePSHiggs->SetBranchAddress( "l1pt",&l1pt);
  treePSHiggs->SetBranchAddress( "l1eta",&l1eta);
  treePSHiggs->SetBranchAddress( "l1pdgId",&l1pdgId);
  treePSHiggs->SetBranchAddress( "l2pt",&l2pt);
  treePSHiggs->SetBranchAddress( "l2eta",&l2eta);
  treePSHiggs->SetBranchAddress( "l2pdgId",&l2pdgId);
  treePSHiggs->SetBranchAddress( "melaLD",&mela);
  treePSHiggs->SetBranchAddress( "melaPSLD",&melaPS);
  treePSHiggs->SetBranchAddress( "melaSpinTwoMinimal",&melaSpinTwoMinimal);
  treePSHiggs->SetBranchAddress( "melaPtY",&melaPtY);
  treePSHiggs->SetBranchAddress( "melaCosThetaStar",&melaCosThetaStar);
  treePSHiggs->SetBranchAddress( "melaCosTheta1",&melaCosTheta1);
  treePSHiggs->SetBranchAddress( "melaCosTheta2",&melaCosTheta2);
  treePSHiggs->SetBranchAddress( "melaPhi",&melaPhi);
  treePSHiggs->SetBranchAddress( "melaPhiStar1",&melaPhiStar1);
  treePSHiggs->SetBranchAddress( "event" ,&event);
  treePSHiggs->SetBranchAddress( "run",&run);
  treePSHiggs->SetBranchAddress( "lumi",&lumi);
  treePSHiggs->SetBranchAddress( "genhiggsmassweight",&higgswgt);
        
  for (int i = 0; i < treePSHiggs->GetEntries(); i++) {
    treePSHiggs->GetEntry(i);
    bool existsAlready = false;
       
    if (existsAlready) {
      std::cout << "Run : " << run << " Lumi : " << lumi << " Event : " << event << std::endl;
    } 
 
    if (!existsAlready) {


      float weight = 1;
      float weighterr = 1;

      weight  *= getBkgXsec(1125)/evt_8TeV(8125);
      weighterr  *= getPUWeight((int)numsim, 1);
      weight *= getSF(l1pt, l1eta, l1pdgId);
      weight *= getSF(l2pt, l2eta, l2pdgId);
      weight *= getSF(l3pt, l3eta, l3pdgId);
      weight *= getSF(l4pt, l4eta, l4pdgId);

      weighterr *= getSF(l1pt, l1eta, l1pdgId);
      weighterr *= getSF(l2pt, l2eta, l2pdgId);
      weighterr *= getSF(l3pt, l3eta, l3pdgId);
      weighterr *= getSF(l4pt, l4eta, l4pdgId);

      //make plot within 10GeV mass window around the higgs
      if (mass > 120 && mass < 130) {
        histPSMELA_PSHiggs->Fill(melaPS, weight);
        histCosThetaStar_PSHiggs->Fill(melaCosThetaStar, weight);
        histCosTheta1_PSHiggs->Fill(melaCosTheta1, weight);
        histCosTheta2_PSHiggs->Fill(melaCosTheta2, weight);
        histPhi_PSHiggs->Fill(melaPhi, weight);
        histPhiStar1_PSHiggs->Fill(melaPhiStar1, weight);
        histZ1Mass_PSHiggs->Fill(z1mass, weight);
        histZ2Mass_PSHiggs->Fill(z2mass, weight);
        if (mela > 0.5) {
          histPSMELA_AfterMELACut_PSHiggs->Fill(melaPS,weight);
          histCosThetaStar_AfterMELACut_PSHiggs->Fill(melaCosThetaStar, weight);
          histCosTheta1_AfterMELACut_PSHiggs->Fill(melaCosTheta1, weight);
          histCosTheta2_AfterMELACut_PSHiggs->Fill(melaCosTheta2, weight);
          histPhi_AfterMELACut_PSHiggs->Fill(melaPhi, weight);
          histPhiStar1_AfterMELACut_PSHiggs->Fill(melaPhiStar1, weight);
          histZ1Mass_AfterMELACut_PSHiggs->Fill(z1mass, weight);
          histZ2Mass_AfterMELACut_PSHiggs->Fill(z2mass, weight);
        }
      }

    }
  }


  //*******************************************************************
  // Loop Over ZZBkg
  //*******************************************************************


  for (int f=0; f < ZZBkgfilenames.size(); ++f) {

    TFile* tmpfileZZBkg = new TFile((ZZBkgfilenames[f].first).c_str());
    TTree* treeZZBkg = (TTree*)tmpfileZZBkg->Get("zz4lTree/probe_tree");
    cout << "load bkg file: " << ZZBkgfilenames[f].first << " , " << ZZBkgfilenames[f].second << " : " << treeZZBkg->GetEntries() << endl;

    treeZZBkg->SetBranchAddress( "channel", &channel);
    treeZZBkg->SetBranchAddress( "z1mass",&z1mass);
    treeZZBkg->SetBranchAddress( "z2mass",&z2mass);
    treeZZBkg->SetBranchAddress("mass" ,&mass);
    treeZZBkg->SetBranchAddress( "numTrueInteractions",&numsim);
    treeZZBkg->SetBranchAddress("l3pt" ,&l3pt);
    treeZZBkg->SetBranchAddress("l3eta" ,&l3eta);
    treeZZBkg->SetBranchAddress("l3idNew" ,&l3id);
    treeZZBkg->SetBranchAddress("l3pfIsoComb04EACorr" ,&l3iso);
    treeZZBkg->SetBranchAddress("l3pdgId" ,&l3pdgId);
    treeZZBkg->SetBranchAddress( "l4pt",&l4pt);
    treeZZBkg->SetBranchAddress("l4eta" ,&l4eta);
    treeZZBkg->SetBranchAddress("l4idNew" ,&l4id);
    treeZZBkg->SetBranchAddress( "l4pfIsoComb04EACorr",&l4iso);
    treeZZBkg->SetBranchAddress("l4pdgId" ,&l4pdgId);
    treeZZBkg->SetBranchAddress( "l1pt",&l1pt);
    treeZZBkg->SetBranchAddress( "l1eta",&l1eta);
    treeZZBkg->SetBranchAddress( "l1pdgId",&l1pdgId);
    treeZZBkg->SetBranchAddress( "l2pt",&l2pt);
    treeZZBkg->SetBranchAddress( "l2eta",&l2eta);
    treeZZBkg->SetBranchAddress( "l2pdgId",&l2pdgId);
    treeZZBkg->SetBranchAddress( "melaLD",&mela);
    treeZZBkg->SetBranchAddress( "melaPSLD",&melaPS);
    treeZZBkg->SetBranchAddress( "melaSpinTwoMinimal",&melaSpinTwoMinimal);
    treeZZBkg->SetBranchAddress( "melaPtY",&melaPtY);
    treeZZBkg->SetBranchAddress( "melaCosThetaStar",&melaCosThetaStar);
    treeZZBkg->SetBranchAddress( "melaCosTheta1",&melaCosTheta1);
    treeZZBkg->SetBranchAddress( "melaCosTheta2",&melaCosTheta2);
    treeZZBkg->SetBranchAddress( "melaPhi",&melaPhi);
    treeZZBkg->SetBranchAddress( "melaPhiStar1",&melaPhiStar1);
    treeZZBkg->SetBranchAddress( "event" ,&event);
    treeZZBkg->SetBranchAddress( "run",&run);
    treeZZBkg->SetBranchAddress( "lumi",&lumi);
    treeZZBkg->SetBranchAddress( "genhiggsmassweight",&higgswgt);
        
    for (int i = 0; i < treeZZBkg->GetEntries(); i++) {
      treeZZBkg->GetEntry(i);
      bool existsAlready = false;
       
      if (existsAlready) {
        std::cout << "Run : " << run << " Lumi : " << lumi << " Event : " << event << std::endl;
      } 
 
      if (!existsAlready) {


        float weight = ZZBkgfilenames[f].second;
        float weighterr = 1.0;

        weight  *= getPUWeight((int)numsim, 1);
        weight *= getSF(l1pt, l1eta, l1pdgId);
        weight *= getSF(l2pt, l2eta, l2pdgId);
        weight *= getSF(l3pt, l3eta, l3pdgId);
        weight *= getSF(l4pt, l4eta, l4pdgId);

        weighterr *= getSF(l1pt, l1eta, l1pdgId);
        weighterr *= getSF(l2pt, l2eta, l2pdgId);
        weighterr *= getSF(l3pt, l3eta, l3pdgId);
        weighterr *= getSF(l4pt, l4eta, l4pdgId);

        //make plot within 10GeV mass window around the higgs
        if (mass > 120 && mass < 130) {
          histPSMELA_ZZBkg->Fill(melaPS, weight);
          histCosThetaStar_ZZBkg->Fill(melaCosThetaStar, weight);
          histCosTheta1_ZZBkg->Fill(melaCosTheta1, weight);
          histCosTheta2_ZZBkg->Fill(melaCosTheta2, weight);
          histPhi_ZZBkg->Fill(melaPhi, weight);
          histPhiStar1_ZZBkg->Fill(melaPhiStar1, weight);
          histZ1Mass_ZZBkg->Fill(z1mass, weight);
          histZ2Mass_ZZBkg->Fill(z2mass, weight);
          if (mela > 0.5) {
            histPSMELA_AfterMELACut_ZZBkg->Fill(melaPS,weight);
            histCosThetaStar_AfterMELACut_ZZBkg->Fill(melaCosThetaStar, weight);
            histCosTheta1_AfterMELACut_ZZBkg->Fill(melaCosTheta1, weight);
            histCosTheta2_AfterMELACut_ZZBkg->Fill(melaCosTheta2, weight);
            histPhi_AfterMELACut_ZZBkg->Fill(melaPhi, weight);
            histPhiStar1_AfterMELACut_ZZBkg->Fill(melaPhiStar1, weight);
            histZ1Mass_AfterMELACut_ZZBkg->Fill(z1mass, weight);
            histZ2Mass_AfterMELACut_ZZBkg->Fill(z2mass, weight);
            //cout << "fill : " << melaPS << " " << weight << endl;
          }
        }

      }
    }
    //tmpfileZZBkg->Close();
  }


  //*******************************************************************
  // Loop Over Data
  //*******************************************************************
  vector<UInt_t> runList;
  vector<UInt_t> lumiList;
  vector<UInt_t> eventList;
  runList.clear();
  lumiList.clear();
  eventList.clear();

  for (int f=0; f < filenamesData.size(); ++f) {

    TFile* tmpfileData = new TFile((filenamesData[f]).c_str());
    TTree* treeData = (TTree*)tmpfileData->Get("zz4lTree/probe_tree");
    cout << "load data file: " << filenamesData[f] << " : " << treeData->GetEntries() << endl;

    treeData->SetBranchAddress( "channel", &channel);
    treeData->SetBranchAddress( "z1mass",&z1mass);
    treeData->SetBranchAddress( "z2mass",&z2mass);
    treeData->SetBranchAddress("mass" ,&mass);
    treeData->SetBranchAddress( "numTrueInteractions",&numsim);
    treeData->SetBranchAddress("l3pt" ,&l3pt);
    treeData->SetBranchAddress("l3eta" ,&l3eta);
    treeData->SetBranchAddress("l3idNew" ,&l3id);
    treeData->SetBranchAddress("l3pfIsoComb04EACorr" ,&l3iso);
    treeData->SetBranchAddress("l3pdgId" ,&l3pdgId);
    treeData->SetBranchAddress( "l4pt",&l4pt);
    treeData->SetBranchAddress("l4eta" ,&l4eta);
    treeData->SetBranchAddress("l4idNew" ,&l4id);
    treeData->SetBranchAddress( "l4pfIsoComb04EACorr",&l4iso);
    treeData->SetBranchAddress("l4pdgId" ,&l4pdgId);
    treeData->SetBranchAddress( "l1pt",&l1pt);
    treeData->SetBranchAddress( "l1eta",&l1eta);
    treeData->SetBranchAddress( "l1pdgId",&l1pdgId);
    treeData->SetBranchAddress( "l2pt",&l2pt);
    treeData->SetBranchAddress( "l2eta",&l2eta);
    treeData->SetBranchAddress( "l2pdgId",&l2pdgId);
    treeData->SetBranchAddress( "melaLD",&mela);
    treeData->SetBranchAddress( "melaPSLD",&melaPS);
    treeData->SetBranchAddress( "melaSpinTwoMinimal",&melaSpinTwoMinimal);
    treeData->SetBranchAddress( "melaPtY",&melaPtY);
    treeData->SetBranchAddress( "melaCosThetaStar",&melaCosThetaStar);
    treeData->SetBranchAddress( "melaCosTheta1",&melaCosTheta1);
    treeData->SetBranchAddress( "melaCosTheta2",&melaCosTheta2);
    treeData->SetBranchAddress( "melaPhi",&melaPhi);
    treeData->SetBranchAddress( "melaPhiStar1",&melaPhiStar1);
    treeData->SetBranchAddress( "event" ,&event);
    treeData->SetBranchAddress( "run",&run);
    treeData->SetBranchAddress( "lumi",&lumi);
    treeData->SetBranchAddress( "genhiggsmassweight",&higgswgt);
        
    for (int i = 0; i < treeData->GetEntries(); i++) {
      treeData->GetEntry(i);
      bool existsAlready = false;

      for (UInt_t k=0; k < runList.size(); ++k) {
        if (runList[k] == run && lumiList[k] == lumi && eventList[k] == event) {
          existsAlready = true;
          break;
        }
      }

      if (existsAlready) {
        std::cout << "Exists: Run : " << run << " Lumi : " << lumi << " Event : " << event << std::endl;
        continue;
      } 
 
      //make plot within 10GeV mass window around the higgs
      if (mass > 120 && mass < 130) {
        histPSMELA_Data->Fill(melaPS);
        histCosThetaStar_Data->Fill(melaCosThetaStar);
        histCosTheta1_Data->Fill(melaCosTheta1);
        histCosTheta2_Data->Fill(melaCosTheta2);
        histPhi_Data->Fill(melaPhi);
        histPhiStar1_Data->Fill(melaPhiStar1);
        histZ1Mass_Data->Fill(z1mass);
        histZ2Mass_Data->Fill(z2mass);
        if (mela > 0.5) {
          histPSMELA_AfterMELACut_Data->Fill(melaPS);
          histCosThetaStar_AfterMELACut_Data->Fill(melaCosThetaStar);
          histCosTheta1_AfterMELACut_Data->Fill(melaCosTheta1);
          histCosTheta2_AfterMELACut_Data->Fill(melaCosTheta2);
          histPhi_AfterMELACut_Data->Fill(melaPhi);
          histPhiStar1_AfterMELACut_Data->Fill(melaPhiStar1);
          histZ1Mass_AfterMELACut_Data->Fill(z1mass);
          histZ2Mass_AfterMELACut_Data->Fill(z2mass);
         }
      }

    }
    //tmpfileData->Close();
  }



  //*******************************************************************
  // Make Plots
  //*******************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
//   gStyle->SetStatOpt(0);


  //*******************************************************************
  // Compare With and without cut on MELA
  //*******************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPSMELA_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histPSMELA_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histPSMELA_SMHiggs->SetLineColor(kRed);
  histPSMELA_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histPSMELA_SMHiggs);
  NormalizeHist(histPSMELA_AfterMELACut_SMHiggs);
  histPSMELA_SMHiggs->SetMaximum( 1.2*max(histPSMELA_SMHiggs->GetMaximum(),histPSMELA_AfterMELACut_SMHiggs->GetMaximum()));
  histPSMELA_SMHiggs->Draw("hist");
  histPSMELA_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("PSMELA_SMHiggs_MELACutImpact.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPSMELA_ZZBkg, "No Cuts", "L");
  legend->AddEntry(histPSMELA_AfterMELACut_ZZBkg, "with MELA > 0.5 Cut", "L");

  histPSMELA_ZZBkg->SetLineColor(kRed);
  histPSMELA_AfterMELACut_ZZBkg->SetLineColor(kBlue);

  NormalizeHist(histPSMELA_ZZBkg);
  NormalizeHist(histPSMELA_AfterMELACut_ZZBkg);
  histPSMELA_ZZBkg->SetMaximum( 1.2*max(histPSMELA_ZZBkg->GetMaximum(),histPSMELA_AfterMELACut_ZZBkg->GetMaximum()));
  histPSMELA_ZZBkg->Draw("hist");
  histPSMELA_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("PSMELA_ZZBkg_MELACutImpact.gif");



  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histCosThetaStar_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histCosThetaStar_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histCosThetaStar_SMHiggs->SetLineColor(kRed);
  histCosThetaStar_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histCosThetaStar_SMHiggs);
  NormalizeHist(histCosThetaStar_AfterMELACut_SMHiggs);
  histCosThetaStar_SMHiggs->SetMaximum( 1.2*max(histCosThetaStar_SMHiggs->GetMaximum(),histCosThetaStar_AfterMELACut_SMHiggs->GetMaximum()));
  histCosThetaStar_SMHiggs->Draw("hist");
  histCosThetaStar_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("CosThetaStar_SMHiggs_MELACutImpact.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histCosTheta1_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histCosTheta1_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histCosTheta1_SMHiggs->SetLineColor(kRed);
  histCosTheta1_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histCosTheta1_SMHiggs);
  NormalizeHist(histCosTheta1_AfterMELACut_SMHiggs);
  histCosTheta1_SMHiggs->SetMaximum( 1.2*max(histCosTheta1_SMHiggs->GetMaximum(),histCosTheta1_AfterMELACut_SMHiggs->GetMaximum()));
  histCosTheta1_SMHiggs->Draw("hist");
  histCosTheta1_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("CosTheta1_SMHiggs_MELACutImpact.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histCosTheta2_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histCosTheta2_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histCosTheta2_SMHiggs->SetLineColor(kRed);
  histCosTheta2_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histCosTheta2_SMHiggs);
  NormalizeHist(histCosTheta2_AfterMELACut_SMHiggs);
  histCosTheta2_SMHiggs->SetMaximum( 1.2*max(histCosTheta2_SMHiggs->GetMaximum(),histCosTheta2_AfterMELACut_SMHiggs->GetMaximum()));
  histCosTheta2_SMHiggs->Draw("hist");
  histCosTheta2_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("CosTheta2_SMHiggs_MELACutImpact.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhi_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histPhi_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histPhi_SMHiggs->SetLineColor(kRed);
  histPhi_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histPhi_SMHiggs);
  NormalizeHist(histPhi_AfterMELACut_SMHiggs);
  histPhi_SMHiggs->SetMaximum( 1.2*max(histPhi_SMHiggs->GetMaximum(),histPhi_AfterMELACut_SMHiggs->GetMaximum()));
  histPhi_SMHiggs->Draw("hist");
  histPhi_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Phi_SMHiggs_MELACutImpact.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhiStar1_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histPhiStar1_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histPhiStar1_SMHiggs->SetLineColor(kRed);
  histPhiStar1_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histPhiStar1_SMHiggs);
  NormalizeHist(histPhiStar1_AfterMELACut_SMHiggs);
  histPhiStar1_SMHiggs->SetMaximum( 1.2*max(histPhiStar1_SMHiggs->GetMaximum(),histPhiStar1_AfterMELACut_SMHiggs->GetMaximum()));
  histPhiStar1_SMHiggs->Draw("hist");
  histPhiStar1_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("PhiStar1_SMHiggs_MELACutImpact.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histZ1Mass_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histZ1Mass_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histZ1Mass_SMHiggs->SetLineColor(kRed);
  histZ1Mass_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histZ1Mass_SMHiggs);
  NormalizeHist(histZ1Mass_AfterMELACut_SMHiggs);
  histZ1Mass_SMHiggs->SetMaximum( 1.2*max(histZ1Mass_SMHiggs->GetMaximum(),histZ1Mass_AfterMELACut_SMHiggs->GetMaximum()));
  histZ1Mass_SMHiggs->Draw("hist");
  histZ1Mass_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Z1Mass_SMHiggs_MELACutImpact.gif");

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histZ2Mass_SMHiggs, "No Cuts", "L");
  legend->AddEntry(histZ2Mass_AfterMELACut_SMHiggs, "with MELA > 0.5 Cut", "L");

  histZ2Mass_SMHiggs->SetLineColor(kRed);
  histZ2Mass_AfterMELACut_SMHiggs->SetLineColor(kBlue);

  NormalizeHist(histZ2Mass_SMHiggs);
  NormalizeHist(histZ2Mass_AfterMELACut_SMHiggs);
  histZ2Mass_SMHiggs->SetMaximum( 1.2*max(histZ2Mass_SMHiggs->GetMaximum(),histZ2Mass_AfterMELACut_SMHiggs->GetMaximum()));
  histZ2Mass_SMHiggs->Draw("hist");
  histZ2Mass_AfterMELACut_SMHiggs->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Z2Mass_SMHiggs_MELACutImpact.gif");


  //*******************************************************************
  // Compare PSMELA distribution for signal and bkg
  //*******************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPSMELA_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histPSMELA_PSHiggs, "0- Higgs", "L");
  legend->AddEntry(histPSMELA_ZZBkg, "ZZ Bkg", "L");

  histPSMELA_SMHiggs->SetLineColor(kRed);
  histPSMELA_PSHiggs->SetLineColor(kBlue);
  histPSMELA_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histPSMELA_SMHiggs);
  NormalizeHist(histPSMELA_PSHiggs);
  NormalizeHist(histPSMELA_ZZBkg);
  histPSMELA_SMHiggs->SetMaximum( 1.2*max(max(histPSMELA_SMHiggs->GetMaximum(),histPSMELA_PSHiggs->GetMaximum()),histPSMELA_ZZBkg->GetMaximum()));
  histPSMELA_SMHiggs->Draw("hist");
  histPSMELA_PSHiggs->Draw("hist,same");
  histPSMELA_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("PSMELA_SignalVsBkg_NoMELACut.gif");

  //*******************************************************************
  // Compare PSMELA distribution for signal and bkg
  //*******************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPSMELA_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histPSMELA_AfterMELACut_PSHiggs, "0- Higgs", "L");
  legend->AddEntry(histPSMELA_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histPSMELA_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histPSMELA_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histPSMELA_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histPSMELA_AfterMELACut_SMHiggs);
  NormalizeHist(histPSMELA_AfterMELACut_PSHiggs);
  NormalizeHist(histPSMELA_AfterMELACut_ZZBkg);
  histPSMELA_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histPSMELA_AfterMELACut_SMHiggs->GetMaximum(),histPSMELA_AfterMELACut_PSHiggs->GetMaximum()),histPSMELA_AfterMELACut_ZZBkg->GetMaximum()));
  histPSMELA_AfterMELACut_SMHiggs->Draw("hist");
  histPSMELA_AfterMELACut_PSHiggs->Draw("hist,same");
  histPSMELA_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("PSMELA_SignalVsBkg_AfterMELACut.gif");

  //*******************************************************************
  // Compare Angles distribution for signal and bkg
  //*******************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histCosThetaStar_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histCosThetaStar_AfterMELACut_PSHiggs, "0- Higgs", "L");
  //legend->AddEntry(histCosThetaStar_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histCosThetaStar_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histCosThetaStar_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histCosThetaStar_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histCosThetaStar_AfterMELACut_SMHiggs);
  NormalizeHist(histCosThetaStar_AfterMELACut_PSHiggs);
  NormalizeHist(histCosThetaStar_AfterMELACut_ZZBkg);
  histCosThetaStar_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histCosThetaStar_AfterMELACut_SMHiggs->GetMaximum(),histCosThetaStar_AfterMELACut_PSHiggs->GetMaximum()),histCosThetaStar_AfterMELACut_ZZBkg->GetMaximum()));
  histCosThetaStar_AfterMELACut_SMHiggs->Draw("hist");
  histCosThetaStar_AfterMELACut_PSHiggs->Draw("hist,same");
  //histCosThetaStar_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("CosThetaStar_SignalVsBkg_AfterMELACut.gif");



  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histCosTheta1_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histCosTheta1_AfterMELACut_PSHiggs, "0- Higgs", "L");
  //legend->AddEntry(histCosTheta1_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histCosTheta1_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histCosTheta1_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histCosTheta1_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histCosTheta1_AfterMELACut_SMHiggs);
  NormalizeHist(histCosTheta1_AfterMELACut_PSHiggs);
  NormalizeHist(histCosTheta1_AfterMELACut_ZZBkg);
  histCosTheta1_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histCosTheta1_AfterMELACut_SMHiggs->GetMaximum(),histCosTheta1_AfterMELACut_PSHiggs->GetMaximum()),histCosTheta1_AfterMELACut_ZZBkg->GetMaximum()));
  histCosTheta1_AfterMELACut_SMHiggs->Draw("hist");
  histCosTheta1_AfterMELACut_PSHiggs->Draw("hist,same");
  //histCosTheta1_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("CosTheta1_SignalVsBkg_AfterMELACut.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histCosTheta2_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histCosTheta2_AfterMELACut_PSHiggs, "0- Higgs", "L");
  //legend->AddEntry(histCosTheta2_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histCosTheta2_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histCosTheta2_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histCosTheta2_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histCosTheta2_AfterMELACut_SMHiggs);
  NormalizeHist(histCosTheta2_AfterMELACut_PSHiggs);
  NormalizeHist(histCosTheta2_AfterMELACut_ZZBkg);
  histCosTheta2_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histCosTheta2_AfterMELACut_SMHiggs->GetMaximum(),histCosTheta2_AfterMELACut_PSHiggs->GetMaximum()),histCosTheta2_AfterMELACut_ZZBkg->GetMaximum()));
  histCosTheta2_AfterMELACut_SMHiggs->Draw("hist");
  histCosTheta2_AfterMELACut_PSHiggs->Draw("hist,same");
  //histCosTheta2_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("CosTheta2_SignalVsBkg_AfterMELACut.gif");



  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhi_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histPhi_AfterMELACut_PSHiggs, "0- Higgs", "L");
  //legend->AddEntry(histPhi_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histPhi_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histPhi_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histPhi_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histPhi_AfterMELACut_SMHiggs);
  NormalizeHist(histPhi_AfterMELACut_PSHiggs);
  NormalizeHist(histPhi_AfterMELACut_ZZBkg);
  histPhi_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histPhi_AfterMELACut_SMHiggs->GetMaximum(),histPhi_AfterMELACut_PSHiggs->GetMaximum()),histPhi_AfterMELACut_ZZBkg->GetMaximum()));
  histPhi_AfterMELACut_SMHiggs->Draw("hist");
  histPhi_AfterMELACut_PSHiggs->Draw("hist,same");
  //histPhi_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Phi_SignalVsBkg_AfterMELACut.gif");



  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histPhiStar1_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histPhiStar1_AfterMELACut_PSHiggs, "0- Higgs", "L");
  //legend->AddEntry(histPhiStar1_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histPhiStar1_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histPhiStar1_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histPhiStar1_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histPhiStar1_AfterMELACut_SMHiggs);
  NormalizeHist(histPhiStar1_AfterMELACut_PSHiggs);
  NormalizeHist(histPhiStar1_AfterMELACut_ZZBkg);
  histPhiStar1_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histPhiStar1_AfterMELACut_SMHiggs->GetMaximum(),histPhiStar1_AfterMELACut_PSHiggs->GetMaximum()),histPhiStar1_AfterMELACut_ZZBkg->GetMaximum()));
  histPhiStar1_AfterMELACut_SMHiggs->Draw("hist");
  histPhiStar1_AfterMELACut_PSHiggs->Draw("hist,same");
  //histPhiStar1_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("PhiStar1_SignalVsBkg_AfterMELACut.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histZ1Mass_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histZ1Mass_AfterMELACut_PSHiggs, "0- Higgs", "L");
  //legend->AddEntry(histZ1Mass_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histZ1Mass_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histZ1Mass_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histZ1Mass_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histZ1Mass_AfterMELACut_SMHiggs);
  NormalizeHist(histZ1Mass_AfterMELACut_PSHiggs);
  NormalizeHist(histZ1Mass_AfterMELACut_ZZBkg);
  histZ1Mass_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histZ1Mass_AfterMELACut_SMHiggs->GetMaximum(),histZ1Mass_AfterMELACut_PSHiggs->GetMaximum()),histZ1Mass_AfterMELACut_ZZBkg->GetMaximum()));
  histZ1Mass_AfterMELACut_SMHiggs->Draw("hist");
  histZ1Mass_AfterMELACut_PSHiggs->Draw("hist,same");
  //histZ1Mass_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Z1Mass_SignalVsBkg_AfterMELACut.gif");


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histZ2Mass_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histZ2Mass_AfterMELACut_PSHiggs, "0- Higgs", "L");
  //legend->AddEntry(histZ2Mass_AfterMELACut_ZZBkg, "ZZ Bkg", "L");

  histZ2Mass_AfterMELACut_SMHiggs->SetLineColor(kRed);
  histZ2Mass_AfterMELACut_PSHiggs->SetLineColor(kBlue);
  histZ2Mass_AfterMELACut_ZZBkg->SetLineColor(kGreen+2);

  NormalizeHist(histZ2Mass_AfterMELACut_SMHiggs);
  NormalizeHist(histZ2Mass_AfterMELACut_PSHiggs);
  NormalizeHist(histZ2Mass_AfterMELACut_ZZBkg);
  histZ2Mass_AfterMELACut_SMHiggs->SetMaximum( 1.2*max(max(histZ2Mass_AfterMELACut_SMHiggs->GetMaximum(),histZ2Mass_AfterMELACut_PSHiggs->GetMaximum()),histZ2Mass_AfterMELACut_ZZBkg->GetMaximum()));
  histZ2Mass_AfterMELACut_SMHiggs->Draw("hist");
  histZ2Mass_AfterMELACut_PSHiggs->Draw("hist,same");
  //histZ2Mass_AfterMELACut_ZZBkg->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Z2Mass_SignalVsBkg_AfterMELACut.gif");



  //*******************************************************************
  // Make Plots for Data
  //*******************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histZ2Mass_AfterMELACut_SMHiggs, "0+ Higgs", "L");
  legend->AddEntry(histZ2Mass_AfterMELACut_PSHiggs, "0- Higgs", "L");



}

void SetTemplateDrawOptions( TH2F* h) {
  
  h->SetMaximum(0.07);
  h->GetXaxis()->SetRangeUser(120,130);
  h->GetXaxis()->SetTitle("m_{4l} [GeV/c^{2}]");
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitle("PseudoMELA");
  h->GetYaxis()->SetTitleOffset(1.02);

}

void PlotTemplates(string filename) {

  TFile *file = new TFile(filename.c_str(),"READ");

  TH2F* hist2D_sig_em = (TH2F*)file->Get("hist2D_sig_em");
  TH2F* hist2D_sig_ee = (TH2F*)file->Get("hist2D_sig_ee");
  TH2F* hist2D_sig_mm = (TH2F*)file->Get("hist2D_sig_mm");
  TH2F* hist2D_pssig_em = (TH2F*)file->Get("hist2D_pssig_em");
  TH2F* hist2D_pssig_ee = (TH2F*)file->Get("hist2D_pssig_ee");
  TH2F* hist2D_pssig_mm = (TH2F*)file->Get("hist2D_pssig_mm");
  TH2F* hist2D_bkg_em = (TH2F*)file->Get("hist2D_bkg_em");
  TH2F* hist2D_bkg_ee = (TH2F*)file->Get("hist2D_bkg_ee");
  TH2F* hist2D_bkg_mm = (TH2F*)file->Get("hist2D_bkg_mm");
  TH2F* hist2D_bkg_em_up = (TH2F*)file->Get("hist2D_bkg_em_up");
  TH2F* hist2D_bkg_ee_up = (TH2F*)file->Get("hist2D_bkg_ee_up");
  TH2F* hist2D_bkg_mm_up = (TH2F*)file->Get("hist2D_bkg_mm_up");


  SetTemplateDrawOptions(hist2D_sig_em);
  SetTemplateDrawOptions(hist2D_sig_ee);
  SetTemplateDrawOptions(hist2D_sig_mm);
  SetTemplateDrawOptions(hist2D_pssig_em);
  SetTemplateDrawOptions(hist2D_pssig_ee);
  SetTemplateDrawOptions(hist2D_pssig_mm);
  SetTemplateDrawOptions(hist2D_bkg_em);
  SetTemplateDrawOptions(hist2D_bkg_ee);
  SetTemplateDrawOptions(hist2D_bkg_mm);
  SetTemplateDrawOptions(hist2D_bkg_em_up);
  SetTemplateDrawOptions(hist2D_bkg_ee_up);
  SetTemplateDrawOptions(hist2D_bkg_mm_up);

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  cv->SetRightMargin(0.10);
  cv->SetRightMargin(0.12);
  hist2D_sig_em->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_sig_em.png");
  hist2D_sig_ee->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_sig_ee.png");
  hist2D_sig_mm->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_sig_mm.png");

  hist2D_pssig_em->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_pssig_em.png");
  hist2D_pssig_ee->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_pssig_ee.png");
  hist2D_pssig_mm->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_pssig_mm.png");

  hist2D_bkg_em->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_bkg_em.png");
  hist2D_bkg_ee->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_bkg_ee.png");
  hist2D_bkg_mm->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_bkg_mm.png");

  hist2D_bkg_em_up->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_bkg_em_up.png");
  hist2D_bkg_ee_up->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_bkg_ee_up.png");
  hist2D_bkg_mm_up->Draw("colz");
  cv->SaveAs("PSMELAMassTemplate_bkg_mm_up.png");


}

void MakePseudoMELAPlots() {

  //set up 
  double lumi_7TeV = 5.05;
  double lumi_8TeV = 12.2;

  bool do7TeV=false;
  init(do7TeV);

  //define bkg samples
  vector<pair<string,double> > zzbkgfilenames;
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id102.root",getBkgXsec(102)/evt_8TeV(102)) );
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id103.root",getBkgXsec(103)/evt_8TeV(103)) );
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id104.root",getBkgXsec(104)/evt_8TeV(104)) );
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id105.root",getBkgXsec(105)/evt_8TeV(105)) );
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id106.root",getBkgXsec(106)/evt_8TeV(106)) );
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id107.root",getBkgXsec(107)/evt_8TeV(107)) );
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id101.root",getBkgXsec(101)/evt_8TeV(101)) );
  zzbkgfilenames.push_back( pair<string,double>("/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/MC/8TeV/yesRegrYesCalib/hzzTree_id100.root",getBkgXsec(100)/evt_8TeV(100)) );

  cout << getBkgXsec(100) << endl;
  cout << evt_8TeV(100) << endl;

  //define data samples
  vector<string> filenamesData;
  filenamesData.push_back( "/data1/sixie/ntuples/step2/hcp/2011/hcp.root");
  filenamesData.push_back( "/data1/sixie/ntuples/step2/hcp/2012/hcp.root");


  MakePseudoscalarPlots( "/data1/sixie/ntuples/step2/hcp/2012/hzzTree_id1125.root", 
                         "/data1/sixie/ntuples/step2/hcp/2012/hzzTree_id8125.root", 
                         zzbkgfilenames, filenamesData );

  PlotTemplates("PSMelaVsMass2DShapes.root");

}

