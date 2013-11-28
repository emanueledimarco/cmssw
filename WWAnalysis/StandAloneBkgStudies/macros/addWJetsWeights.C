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

using namespace std;



int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
   }


void setPrMuClaraLP(TH2F* promptH2Mu);
void setPrMuClaraKink(TH2F* promptH2Mu);
void setPrMuClaraKinkLow(TH2F* promptH2Mu);
void setPrMuClaraKinkHigh(TH2F* promptH2Mu);
void setPrMuClaraKink2011B(TH2F* promptH2Mu);

void setPrElClaraLP(TH2F* promptH2El);
void setPrElClaraLH(TH2F* promptH2El);

void setPrElClaraBDT(TH2F* promptH2El);
void setPrElClaraBDTLow(TH2F* promptH2El);
void setPrElClaraBDTHigh(TH2F* promptH2El);
void setPrElClaraBDT2011B(TH2F* promptH2El);

void setPrElEmaLP(TH2F* promptH2El);
void setPrElEmaLH(TH2F* promptH2El);

void setFrMuXieLP(TH2F* fakeH2Mu);
void setFrMuAliciaLP(TH2F* fakeH2Mu);
void setFrMuAliciaKink(TH2F* fakeH2Mu);

void setFrElChiaraLP(TH2F* fakeH2El);
void setFrElChiaraLH(TH2F* fakeH2El);
void setFrElXieBDT(TH2F* fakeH2El);
void setFrElChiaraBDT2011B(TH2F* fakeH2El);
void setFrElChiaraBDT2011(TH2F* fakeH2El);
void setFrElXieLP(TH2F* fakeH2El);


double fake(double pt,double eta,TH2F* fakeH2);
double prompt(double pt,double eta,TH2F* promptH2);

double eps(double fake){
  return fake/(1-fake);
}

double eta(double prompt){
  return (1-prompt)/prompt;
}

double wFailFail(double eps1,double eps2,
		 double eta1,double eta2){
  double result = -2*eps1*eps2;
  return result;
}

double wPassPass(double eps1,double eps2,
		 double eta1,double eta2){

  double result = -(eps1*eta1+eps2*eta2);
  return result;
}

double wPassFail(double epsPass, double etaPass,
		 double epsFail, double etaFail){
  double result = epsFail+epsPass*etaPass*epsFail;
  return result;		 
}


void setHisto2D(TH2F* fakeHistoMu,TH2F* promptHistoMuLow,TH2F* promptHistoMuHigh,
		TH2F* fakeHistoEl,TH2F* promptHistoElLow,TH2F* promptHistoElHigh);

int addWJetsWeights(string input="latino_ll_LP_test.root",bool addWeight=false) {

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
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleOffset(0.80,"X");
  gStyle->SetTitleOffset(0.80,"Y");
  // ------------------------------ 


  string openMode;
  if(addWeight)
    openMode="update";//xAlicia
  else{
    openMode="read";
  }

  TH1F* h_mll = new TH1F("h_mll","h_mll",25,0,200);
  TH1F* h_ptMax = new TH1F("h_ptMax","h_ptMax",85,15,100);
  TH1F* h_ptMin = new TH1F("h_ptMin","h_ptMin",95,5,70);


  //
  TFile* file = new TFile(input.c_str(),openMode.c_str());
  TTree* tree = (TTree*)file->Get("latino");

  /* ---- OUR BINNING ---- 
  const int ptBinsSizeMu(5);
  float ptBinsMu[ptBinsSizeMu+1];
  ptBinsMu[0] = 10;
  ptBinsMu[1] = 15;
  ptBinsMu[2] = 20;
  ptBinsMu[3] = 25;
  ptBinsMu[4] = 50;
  ptBinsMu[5] = 9999999;

  const int etaBinsSizeMu(2);
  float etaBinsMu[etaBinsSizeMu+1];
  etaBinsMu[0] = 0;
  etaBinsMu[1] = 1.5;
  etaBinsMu[2] = 2.5;

  const int ptBinsSizeEl(5);
  float ptBinsEl[ptBinsSizeEl+1];
  ptBinsEl[0] = 10;
  ptBinsEl[1] = 15;
  ptBinsEl[2] = 20;
  ptBinsEl[3] = 25;
  ptBinsEl[4] = 50;
  ptBinsEl[5] = 9999999;

  const int etaBinsSizeElFake(2);
  float etaBinsElFake[etaBinsSizeElFake+1];
  etaBinsElFake[0] = 0;
  etaBinsElFake[1] = 1.476;
  etaBinsElFake[2] = 2.5;

  const int etaBinsSizeElPrompt(3);
  float etaBinsElPrompt[etaBinsSizeElPrompt+1];
  etaBinsElPrompt[0] = 0;
  etaBinsElPrompt[1] = 1.4442;
  etaBinsElPrompt[2] = 1.566;
  etaBinsElPrompt[3] = 2.5;

  TH2F* fakeMuonH2 = new TH2F("fakeMuonH2","fakeMuonH2",ptBinsSizeMu,ptBinsMu,etaBinsSizeMu,etaBinsMu);
  TH2F* promptMuonH2 = new TH2F("promptMuonH2","promptMuonH2",ptBinsSizeMu,ptBinsMu,etaBinsSizeMu,etaBinsMu);

  TH2F* fakeElH2 = new TH2F("fakeElH2","fakeElH2",ptBinsSizeEl,ptBinsEl,etaBinsSizeElFake,etaBinsElFake);
  TH2F* promptElH2 = new TH2F("promptElH2","promptElH2",ptBinsSizeEl,ptBinsEl,etaBinsSizeElPrompt,etaBinsElPrompt);
  --------- */


  // ---- XIE denominator ---- 
  //FakeRate
  const int ptBinsSizeMuFake(5);
  float ptBinsMuFake[ptBinsSizeMuFake+1];
  ptBinsMuFake[0] = 10;
  ptBinsMuFake[1] = 15;
  ptBinsMuFake[2] = 20;
  ptBinsMuFake[3] = 25;
  ptBinsMuFake[4] = 30;
  ptBinsMuFake[5] = 35;

  
  const int etaBinsSizeMuFake(4);
  float etaBinsMuFake[etaBinsSizeMuFake+1];
  etaBinsMuFake[0] = 0;
  etaBinsMuFake[1] = 1.0;
  etaBinsMuFake[2] = 1.479;
  etaBinsMuFake[3] = 2.0;
  etaBinsMuFake[4] = 2.5;



  const int ptBinsSizeElFake(5);
  float ptBinsElFake[ptBinsSizeElFake+1];
  ptBinsElFake[0] = 10;
  ptBinsElFake[1] = 15;
  ptBinsElFake[2] = 20;
  ptBinsElFake[3] = 25;
  ptBinsElFake[4] = 30;
  ptBinsElFake[5] = 35;


  const int etaBinsSizeElFake(4);
  float etaBinsElFake[etaBinsSizeElFake+1];
  etaBinsElFake[0] = 0;
  etaBinsElFake[1] = 1.0;
  etaBinsElFake[2] = 1.479;
  etaBinsElFake[3] = 2.0;
  etaBinsElFake[4] = 2.5;


  //PromptRate
  const int ptBinsSizeMuPrompt(5);
  float ptBinsMuPrompt[ptBinsSizeMuPrompt+1];
  ptBinsMuPrompt[0] = 10;
  ptBinsMuPrompt[1] = 15;
  ptBinsMuPrompt[2] = 20;
  ptBinsMuPrompt[3] = 25;
  ptBinsMuPrompt[4] = 50;
  ptBinsMuPrompt[5] = 100;


  const int etaBinsSizeMuPrompt(2);
  float etaBinsMuPrompt[etaBinsSizeMuPrompt+1];
  etaBinsMuPrompt[0] = 0;
  etaBinsMuPrompt[1] = 1.5;
  etaBinsMuPrompt[2] = 2.4;




  const int ptBinsSizeElPrompt(5);
  float ptBinsElPrompt[ptBinsSizeElPrompt+1];
  ptBinsElPrompt[0] = 10;
  ptBinsElPrompt[1] = 15;
  ptBinsElPrompt[2] = 20;
  ptBinsElPrompt[3] = 25;
  ptBinsElPrompt[4] = 50;
  ptBinsElPrompt[5] = 100;

  const int etaBinsSizeElPrompt(3);
  float etaBinsElPrompt[etaBinsSizeElPrompt+1];
  etaBinsElPrompt[0] = 0;
  etaBinsElPrompt[1] = 1.4442;
  etaBinsElPrompt[2] = 1.566;
  etaBinsElPrompt[3] = 2.5;




  TH2F* fakeMuonH2 = new TH2F("fakeMuonH2","fakeMuonH2",ptBinsSizeMuFake,ptBinsMuFake,
			      etaBinsSizeMuFake,etaBinsMuFake);
  TH2F* promptMuonH2Low = new TH2F("promptMuonH2Low","promptMuonH2Low",ptBinsSizeMuPrompt,ptBinsMuPrompt,
				   etaBinsSizeMuPrompt,etaBinsMuPrompt);
  TH2F* promptMuonH2High = new TH2F("promptMuonH2High","promptMuonH2High",ptBinsSizeMuPrompt,ptBinsMuPrompt,
				    etaBinsSizeMuPrompt,etaBinsMuPrompt);
  TH2F* promptMuonH2;

  TH2F* fakeElH2 = new TH2F("fakeElH2","fakeElH2",ptBinsSizeElFake,ptBinsElFake,
			    etaBinsSizeElFake,etaBinsElFake);
  TH2F* promptElH2Low = new TH2F("promptElH2Low","promptElH2Low",ptBinsSizeElPrompt,ptBinsElPrompt,
				 etaBinsSizeElPrompt,etaBinsElPrompt);
  TH2F* promptElH2High = new TH2F("promptElH2High","promptElH2High",ptBinsSizeElPrompt,ptBinsElPrompt,
				  etaBinsSizeElPrompt,etaBinsElPrompt);
  TH2F* promptElH2;

  setHisto2D(fakeMuonH2,promptMuonH2Low,promptMuonH2High,
	     fakeElH2,promptElH2Low,promptElH2High);

  enum ch {mm=0,ee=1,em=2,me=3};
  
  float weight;
  float channel,pt1,pt2,eta1,eta2,ptll,nvtx;
  int passTight1A,passTight2A,passTight1B,passTight2B,passTight1,passTight2;
  
  float trigger,pfmet,ppfmet,mll,mpmet,njet,nextra,nbjet,dphilljetjet;
  int zveto,dphiveto,bveto_mu,bveto_ip,sameflav;
  unsigned int run,lumi,event;

  TBranch *newBranch;
  if(addWeight)
    newBranch = tree->Branch("fake2W",&weight,"fakeW/F");//xAlicia

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("lumi",&lumi);
  tree->SetBranchAddress("event",&event);

  tree->SetBranchAddress("channel",&channel);
  tree->SetBranchAddress("pt1",&pt1);
  tree->SetBranchAddress("pt2",&pt2);
  tree->SetBranchAddress("eta1",&eta1);
  tree->SetBranchAddress("eta2",&eta2);
  //tree->SetBranchAddress("passTight1",&passTight1); //for scenario1
  //tree->SetBranchAddress("passTight2",&passTight2);
  //tree->SetBranchAddress("passCB1",&passTight1); //for scenario2
  //tree->SetBranchAddress("passCB2",&passTight2);
  //tree->SetBranchAddress("passLH1",&passTight1); //for scenario3
  //tree->SetBranchAddress("passLH2",&passTight2);
  tree->SetBranchAddress("passBDT1",&passTight1);
  tree->SetBranchAddress("passBDT2",&passTight2);



  tree->SetBranchAddress("trigger",&trigger);
  tree->SetBranchAddress("pfmet",&pfmet);
  tree->SetBranchAddress("ppfmet",&ppfmet);
  tree->SetBranchAddress("mll",&mll);
  tree->SetBranchAddress("nvtx",&nvtx);
  tree->SetBranchAddress("zveto",&zveto);
  tree->SetBranchAddress("mpmet",&mpmet);
  tree->SetBranchAddress("sameflav",&sameflav);
  tree->SetBranchAddress("njet",&njet);
  tree->SetBranchAddress("dphiveto",&dphiveto);
  tree->SetBranchAddress("bveto_mu",&bveto_mu);
  tree->SetBranchAddress("nextra",&nextra);
  tree->SetBranchAddress("bveto_ip",&bveto_ip);
  tree->SetBranchAddress("nbjet",&nbjet);
  tree->SetBranchAddress("ptll",&ptll);
  tree->SetBranchAddress("dphilljetjet",&dphilljetjet);
  //tree->SetBranchAddress("",&);
  
  int mmPP(0),mmPF(0),mmFF(0);
  int mePP(0),mePF(0),meFF(0);
  int emPP(0),emPF(0),emFF(0);
  int eePP(0),eePF(0),eeFF(0);

  double sumOfWeightsMuMu(0.);
  double sumOfWeightsMuEl(0.);
  double sumOfWeightsElMu(0.);
  double sumOfWeightsElEl(0.);
  double sumOfWeightsTrueElFakeMu(0.);
  double sumOfWeightsTrueMuFakeEl(0.);
  double sum00(0),sum10(0.),sum11(0.);
  double checksum(0.);

  unsigned int nentries = tree->GetEntries();
  cout << "nentries: " << nentries << endl;
  unsigned int countEntriesByHand(0);
  unsigned int countEntriesByHand2(0);

  for(unsigned int i = 0; i<nentries; i++){
  //for(unsigned int i = 0; i<1000; i++){
    weight = 1.;//xAlicia
    tree->GetEntry(i);

    // -- new ww selection
    bool ww0JetSel = trigger && pfmet>20 && mll > (12 + 8*sameflav) && zveto && mpmet > (20+(17+nvtx/2.)*sameflav);
    ww0JetSel = ww0JetSel && njet==0 && (dphiveto || !sameflav) && bveto_mu && nextra == 0 && bveto_ip;
    ww0JetSel = ww0JetSel && (pt2 > 15||!sameflav) && ptll > 45;
    //
    // -- old ww selection
    //bool ww0JetSel = trigger && pfmet>20 && mll > 12 && zveto && mpmet > (20+(17+nvtx/2.)*sameflav);
    //ww0JetSel = ww0JetSel && njet==0 && (dphiveto || !sameflav) && bveto_mu && nextra == 0 && bveto_ip;



    // -- new ww selection
    bool ww1JetSel = trigger && pfmet>20 && mll > (12 + 8*sameflav) && zveto && mpmet > (20+(17+nvtx/2.)*sameflav);
    ww1JetSel = ww1JetSel && njet==1 && (dphiveto || !sameflav) && bveto_mu && nextra == 0 && bveto_ip && nbjet==0;
    ww1JetSel = ww1JetSel && (pt2>15 || !sameflav) && ptll > 45;
    //
    // -- old ww selection
    //bool ww1JetSel = trigger && pfmet>20 && mll > 12 && zveto && mpmet > (20+(17+nvtx/2.)*sameflav);
    //ww1JetSel = ww1JetSel && njet==1 && (dphiveto || !sameflav) && bveto_mu && nextra == 0 && bveto_ip && nbjet==0;


    
    bool ww2JetSel = trigger && pfmet>20 && mll > (12 + 8*sameflav) && zveto && mpmet > (20+(17+nvtx/2.)*sameflav);
    ww2JetSel = ww2JetSel && njet==2 && (dphilljetjet< M_PI/180.*165. || !sameflav) && bveto_mu && nextra == 0 && bveto_ip && nbjet==0;
    ww2JetSel = ww2JetSel && (pt2>15 || !sameflav) && ptll > 45;


    //if(!ww0JetSel) continue;
    //if(lumi!=56 || event!=38975488) continue;

    double wSumUp(999999999.);
    double wTrueElFakeMu(999999999.);
    double wTrueMuFakeEl(999999999.);
    //======== HERE IS THE PART COMMON TO ALL CHANNELS
    if(nvtx<=9){
      promptMuonH2 = promptMuonH2Low;
      promptElH2 = promptElH2Low;
    }else{
      promptMuonH2 = promptMuonH2High;
      promptElH2 = promptElH2High;
    }


    double fake1 = fake(pt1,eta1,(channel == mm || channel == me) ? fakeMuonH2 : fakeElH2);
    double fake2 = fake(pt2,eta2,(channel == ee || channel == me) ? fakeElH2 : fakeMuonH2);
    double prompt1 = prompt(pt1,eta1,(channel == mm || channel == me) ? promptMuonH2 : promptElH2);
    double prompt2 = prompt(pt2,eta2,(channel == ee || channel == me) ? promptElH2 : promptMuonH2);

    double Eps1 = eps(fake1);
    double Eps2 = eps(fake2);

    double Eta1 = eta(prompt1);
    double Eta2 = eta(prompt2);
    
    if(passTight1 && passTight2){
      wSumUp = -(Eps1*Eta1+Eps2*Eta2);
    }else if(!passTight1 && !passTight2){
      wSumUp = -2*Eps1*Eps2;
    }else if(passTight1 && !passTight2){
      wSumUp = Eps2+Eps1*Eta1*Eps2;
    }else if(!passTight1 && passTight2){ 
      wSumUp = Eps1+Eps2*Eta2*Eps1;
    }else{
      cout << "ERROR: logical problem" << endl;
      return 1;
    }
    double norm = 1./(1-Eps1*Eta1)/(1-Eps2*Eta2);
    wSumUp = wSumUp*norm;
    weight = wSumUp; //converting to float for saving into tree //xAlicia

    //------ summing the weights only for those events passing the full selection
    if(ww1JetSel) {
      countEntriesByHand2++;
      if(channel == mm) sumOfWeightsMuMu += wSumUp; 
      if(channel == me) {sumOfWeightsMuEl += wSumUp; 
	char* label;
	if(passTight1 && passTight2) label="PP";
	if(!passTight1 && !passTight2) label="FF";
	if((!passTight1 && passTight2) || (passTight1 && !passTight2)) label="PF";
	/*
	  cout << "run,lumi,event,w,type: " 
	  << run << "  "
	  << lumi << "  "
	  << event << "  "
	  << wSumUp << " "
	  << label << endl;
	*/
      }
      if(channel == em) sumOfWeightsElMu += wSumUp; 
      if(channel == ee) sumOfWeightsElEl += wSumUp; 
      
      if(passTight1 && passTight2){
	if(channel == mm) mmPP++;
	if(channel == me) mePP++;
	if(channel == em) emPP++;
	if(channel == ee) eePP++;
      }
      if((passTight1 && !passTight2) || (!passTight1 && passTight2)){
	if(channel == mm) mmPF++;
	if(channel == me) mePF++;
	if(channel == em) emPF++;
	if(channel == ee) eePF++;
      }
      if(!passTight1 && !passTight2){
	if(channel == mm) mmFF++;
	if(channel == me) meFF++;
	if(channel == em) emFF++;
	if(channel == ee) eeFF++;
      }
    }

    /*
    //======== HERE IS THE PART FOR MIXED CHANNEL
    if( channel == em || channel == me){
      double ptEl = (channel == em) ? pt1 : pt2;
      double ptMu = (channel == em) ? pt2 : pt1;

      double psrapEl = (channel == em) ? eta1 : eta2;
      double psrapMu = (channel == em) ? eta2 : eta1;

      double fakeEl = fake(ptEl,psrapEl, fakeElH2);
      double fakeMu = fake(ptMu,psrapMu, fakeMuonH2);
      double promptEl = prompt(ptEl,psrapEl,promptElH2);
      double promptMu = prompt(ptMu,psrapMu, promptMuonH2);

      double epsEl = eps(fakeEl);
      double epsMu = eps(fakeMu);

      double etaEl = eta(promptEl);
      double etaMu = eta(promptMu);

      bool passEl,passMu;
      if(channel == em)
	{passEl = passTight1;passMu = passTight2;}
      else
	{passEl = passTight2; passMu = passTight1;}
      

      if(passEl && passMu){
	wTrueElFakeMu = -epsEl*epsMu;
	wTrueMuFakeEl = -epsMu*epsEl;
      }else if(!passEl && !passMu){
	wTrueElFakeMu = -epsMu*etaMu;
	wTrueMuFakeEl = -epsEl*etaEl;
      }else if(passEl && !passMu){
	wTrueElFakeMu = epsMu;
	wTrueMuFakeEl = epsEl*etaEl*epsMu;
      }else if(!passEl && passMu){ 
	wTrueElFakeMu = epsMu*epsEl*etaMu;
	wTrueMuFakeEl = epsEl;
      }else{
	cout << "ERROR: logical problem" << endl;
	return 1;
      }
      double norm = 1./(1-epsEl*etaEl)/(1-epsMu*etaMu);
      wTrueElFakeMu = wTrueElFakeMu*norm;
      wTrueMuFakeEl = wTrueMuFakeEl*norm;
      //weight = wSF; //TO BE FIXED

      //------ summing the weights only for those events passing the full selection
      if(ww0JetSel) {
	sumOfWeightsTrueElFakeMu += wTrueElFakeMu; 
	sumOfWeightsTrueMuFakeEl += wTrueMuFakeEl; 
      }
    }
    */
    
    countEntriesByHand++;
    //cout << "entry #: " << countEntriesByHand << ", weight " << weight << endl;
    if(addWeight) newBranch->Fill();//xAlicia
  }

  cout << "# total: " << countEntriesByHand << endl;
  cout << "# passing WW 0-jet: " << countEntriesByHand2 << endl;
  

  cout << "sumOfWeights mm: " << sumOfWeightsMuMu << endl;
  cout << "sumOfWeights me: " << sumOfWeightsMuEl << endl;
  cout << "sumOfWeights em: " << sumOfWeightsElMu << endl;
  cout << "sumOfWeights ee: " << sumOfWeightsElEl << endl;

  /*
  cout << endl << "sumOfWeights em+me: " << sumOfWeightsElMu+sumOfWeightsMuEl << endl;
  cout << "sumOfWeights m-FakeEl: " << sumOfWeightsTrueMuFakeEl << endl;
  cout << "sumOfWeights e-FakeMu: " << sumOfWeightsTrueElFakeMu << endl;
  cout << "sumOfWeights e-FakeMu+m-FakeEl: " << (sumOfWeightsTrueElFakeMu+sumOfWeightsTrueMuFakeEl) << endl;
  */

  /*
  float scaleToLP(7.38);
  cout << "=== Extrapolation to LP ===" << endl;
  cout << "mm to LP: " << sumOfWeightsMuMu*scaleToLP << endl;
  cout << "me to LP: " << sumOfWeightsMuEl*scaleToLP << endl;
  cout << "em to LP: " << sumOfWeightsElMu*scaleToLP << endl;
  cout << "ee to LP: " << sumOfWeightsElEl*scaleToLP << endl;

  cout << "m-FakeEl to LP: " << sumOfWeightsTrueMuFakeEl*scaleToLP << endl;
  cout << "e-FakeMu to LP: " << sumOfWeightsTrueElFakeMu*scaleToLP << endl;  

  cout << "checksum mm: " << checksum << endl;
  */

  cout << "mm PP,PF,FF: " << mmPP << " , " << mmPF << " , " << mmFF << endl;
  cout << "me PP,PF,FF: " << mePP << " , " << mePF << " , " << meFF << endl;
  cout << "em PP,PF,FF: " << emPP << " , " << emPF << " , " << emFF << endl;
  cout << "ee PP,PF,FF: " << eePP << " , " << eePF << " , " << eeFF << endl;

  //====================================
  /*
  cout << "mll integral: " << h_mll->Integral() << endl;
  h_mll->Draw();gPad->Update(); Wait();
  cout << "ptMax integral: " << h_ptMax->Integral() << endl;
  h_ptMax->Draw();gPad->Update(); Wait();
  cout << "ptMin integral: " << h_ptMin->Integral() << endl;
  h_ptMin->Draw();gPad->Update(); Wait();
  */

  if(addWeight) tree->Write("",TObject::kOverwrite);//xAlicia
  delete file;//xAlicia
  return 0;
}	      



double fake(double pt,double eta,TH2F* h2){
  double etaMax = h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins());
  double ptMax = h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins());
  double eps = 0.00000001;
  //cout << "etaMax, ptmax: " << etaMax << " , " << ptMax << endl;

  int bin = h2->FindFixBin(std::min(pt,ptMax-eps),std::min(fabs(eta),etaMax-eps));
  double result = h2->GetBinContent(bin);
  //cout << "fake: " << result << endl;
  return result;
}

double prompt(double pt,double eta,TH2F* h2){
  double etaMax = h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins());
  double ptMax = h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins());
  double eps = 0.00000001;

  int bin = h2->FindFixBin(std::min(pt,ptMax-eps),std::min(fabs(eta),etaMax-eps));
  double result = h2->GetBinContent(bin);
  //cout << "prompt: " << result << endl;
  return result;
}

void setHisto2D(TH2F* fakeH2Mu,TH2F* promptH2MuLow,TH2F* promptH2MuHigh,
		TH2F* fakeH2El,TH2F* promptH2ElLow,TH2F* promptH2ElHigh){
  setPrMuClaraKinkLow(promptH2MuLow);    
  setPrMuClaraKinkHigh(promptH2MuHigh);    
  setPrElClaraBDTLow(promptH2ElLow);
  setPrElClaraBDTHigh(promptH2ElHigh);

  setFrMuAliciaKink(fakeH2Mu);
  setFrElChiaraBDT2011(fakeH2El);

}


void setHisto2D(TH2F* fakeH2Mu,TH2F* promptH2Mu,
		TH2F* fakeH2El,TH2F* promptH2El){

  // --- for scenario 1 (using LL0)
  // --- for scenario 1' (using LL1)
  //be careful! these numbers can only be used with LL1 trees
  //setPrMuClaraLP(promptH2Mu);      
  //setFrMuAliciaLP(fakeH2Mu);       
  //setFrElChiaraLP(fakeH2El);
  //setPrElEmaLP(promptH2El);
  ///setPrElClaraLP(promptH2El);

  // --- for scenario 2 (using LL2)
  //be careful! these numbers can only be used with LL2 trees
  //setPrMuClaraKink(promptH2Mu);    
  //setFrMuAliciaKink(fakeH2Mu);
  //setFrElChiaraLP(fakeH2El);
  //setPrElEmaLP(promptH2El);

  // --- for scenario 3 (using LL2)
  //be careful! these numbers can only be used with LL2 trees
  //setPrMuClaraKink(promptH2Mu);    
  //setFrMuAliciaKink(fakeH2Mu);
  //setFrElChiaraLH(fakeH2El);
  //setPrElEmaLH(promptH2El);

  // --- for scenario 4 (using LL2)
  //be careful! these numbers can only be used with LL2 trees
  //setPrMuClaraKink(promptH2Mu);    
  //setFrMuAliciaKink(fakeH2Mu);
  //setFrElXieBDT(fakeH2El);
  //setPrElClaraBDT(promptH2El);
  //
  setPrMuClaraKink2011B(promptH2Mu);    
  setFrMuAliciaKink(fakeH2Mu);
  setFrElChiaraBDT2011B(fakeH2El);
  setPrElClaraBDT2011B(promptH2El);


  // --- others
  //setPrElClaraLH(promptH2El);
  //setFrMuXieLP(fakeH2Mu);


}


void setPrMuClaraLP(TH2F* promptH2Mu){
  promptH2Mu->SetBinContent(1,1,0.740);
  promptH2Mu->SetBinContent(2,1,0.790);
  promptH2Mu->SetBinContent(3,1,0.953);
  promptH2Mu->SetBinContent(4,1,0.985);
  promptH2Mu->SetBinContent(5,1,0.995);
  //
  promptH2Mu->SetBinContent(1,2,0.749);
  promptH2Mu->SetBinContent(2,2,0.773);
  promptH2Mu->SetBinContent(3,2,0.919);
  promptH2Mu->SetBinContent(4,2,0.968);
  promptH2Mu->SetBinContent(5,2,0.988);
}


void setPrMuClaraKink(TH2F* promptH2Mu){
  //be careful! these numbers can only be used with LL2 trees
  promptH2Mu->SetBinContent(1,1,0.741);
  promptH2Mu->SetBinContent(2,1,0.789);
  promptH2Mu->SetBinContent(3,1,0.953);
  promptH2Mu->SetBinContent(4,1,0.985);
  promptH2Mu->SetBinContent(5,1,0.995);
  //
  promptH2Mu->SetBinContent(1,2,0.749);
  promptH2Mu->SetBinContent(2,2,0.773);
  promptH2Mu->SetBinContent(3,2,0.919);
  promptH2Mu->SetBinContent(4,2,0.968);
  promptH2Mu->SetBinContent(5,2,0.988);
}

void setPrMuClaraKinkLow(TH2F* promptH2Mu){
  //be careful! these numbers can only be used with LL2 trees
  promptH2Mu->SetBinContent(1,1,0.731);
  promptH2Mu->SetBinContent(2,1,0.780);
  promptH2Mu->SetBinContent(3,1,0.949);
  promptH2Mu->SetBinContent(4,1,0.984);
  promptH2Mu->SetBinContent(5,1,0.995);
  //
  promptH2Mu->SetBinContent(1,2,0.740);
  promptH2Mu->SetBinContent(2,2,0.767);
  promptH2Mu->SetBinContent(3,2,0.913);
  promptH2Mu->SetBinContent(4,2,0.965);
  promptH2Mu->SetBinContent(5,2,0.989);
}

void setPrMuClaraKinkHigh(TH2F* promptH2Mu){
  //be careful! these numbers can only be used with LL2 trees
  promptH2Mu->SetBinContent(1,1,0.684);
  promptH2Mu->SetBinContent(2,1,0.718);
  promptH2Mu->SetBinContent(3,1,0.921);
  promptH2Mu->SetBinContent(4,1,0.976);
  promptH2Mu->SetBinContent(5,1,0.992);
  //
  promptH2Mu->SetBinContent(1,2,0.658);
  promptH2Mu->SetBinContent(2,2,0.688);
  promptH2Mu->SetBinContent(3,2,0.856);
  promptH2Mu->SetBinContent(4,2,0.938);
  promptH2Mu->SetBinContent(5,2,0.982);
}



void setPrMuClaraKink2011B(TH2F* promptH2Mu){
  //be careful! these numbers can only be used with LL2 trees
  promptH2Mu->SetBinContent(1,1,0.690);
  promptH2Mu->SetBinContent(2,1,0.731);
  promptH2Mu->SetBinContent(3,1,0.934);
  promptH2Mu->SetBinContent(4,1,0.980);
  promptH2Mu->SetBinContent(5,1,0.993);
  //
  promptH2Mu->SetBinContent(1,2,0.686);
  promptH2Mu->SetBinContent(2,2,0.718);
  promptH2Mu->SetBinContent(3,2,0.886);
  promptH2Mu->SetBinContent(4,2,0.949);
  promptH2Mu->SetBinContent(5,2,0.985);
}


void setPrElClaraLP(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.58);
  promptH2El->SetBinContent(2,1,0.65);
  promptH2El->SetBinContent(3,1,0.78);
  promptH2El->SetBinContent(4,1,0.88);
  promptH2El->SetBinContent(5,1,0.91);
  //
  promptH2El->SetBinContent(1,2,0.41);
  promptH2El->SetBinContent(2,2,0.50);
  promptH2El->SetBinContent(3,2,0.70);
  promptH2El->SetBinContent(4,2,0.80);
  promptH2El->SetBinContent(5,2,0.88);
  //
  promptH2El->SetBinContent(1,3,0.40);
  promptH2El->SetBinContent(2,3,0.48);
  promptH2El->SetBinContent(3,3,0.71);
  promptH2El->SetBinContent(4,3,0.81);
  promptH2El->SetBinContent(5,3,0.87);
}

void setPrElClaraLH(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.58);
  promptH2El->SetBinContent(2,1,0.65);
  promptH2El->SetBinContent(3,1,0.78);
  promptH2El->SetBinContent(4,1,0.88);
  promptH2El->SetBinContent(5,1,0.91);
  //
  promptH2El->SetBinContent(1,2,0.41);
  promptH2El->SetBinContent(2,2,0.50);
  promptH2El->SetBinContent(3,2,0.70);
  promptH2El->SetBinContent(4,2,0.80);
  promptH2El->SetBinContent(5,2,0.88);
  //
  promptH2El->SetBinContent(1,3,0.40);
  promptH2El->SetBinContent(2,3,0.48);
  promptH2El->SetBinContent(3,3,0.71);
  promptH2El->SetBinContent(4,3,0.81);
  promptH2El->SetBinContent(5,3,0.87);
}

void setPrElClaraBDT(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.61);
  promptH2El->SetBinContent(2,1,0.68);
  promptH2El->SetBinContent(3,1,0.78);
  promptH2El->SetBinContent(4,1,0.88);
  promptH2El->SetBinContent(5,1,0.90);
  //
  promptH2El->SetBinContent(1,2,0.43);
  promptH2El->SetBinContent(2,2,0.54);
  promptH2El->SetBinContent(3,2,0.72);
  promptH2El->SetBinContent(4,2,0.83);
  promptH2El->SetBinContent(5,2,0.88);
  //
  promptH2El->SetBinContent(1,3,0.38);
  promptH2El->SetBinContent(2,3,0.46);
  promptH2El->SetBinContent(3,3,0.68);
  promptH2El->SetBinContent(4,3,0.79);
  promptH2El->SetBinContent(5,3,0.84);
}

void setPrElClaraBDTLow(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.59);
  promptH2El->SetBinContent(2,1,0.67);
  promptH2El->SetBinContent(3,1,0.78);
  promptH2El->SetBinContent(4,1,0.87);
  promptH2El->SetBinContent(5,1,0.91);
  //
  promptH2El->SetBinContent(1,2,0.47);
  promptH2El->SetBinContent(2,2,0.54);
  promptH2El->SetBinContent(3,2,0.71);
  promptH2El->SetBinContent(4,2,0.82);
  promptH2El->SetBinContent(5,2,0.87);
  //
  promptH2El->SetBinContent(1,3,0.37);
  promptH2El->SetBinContent(2,3,0.50);
  promptH2El->SetBinContent(3,3,0.68);
  promptH2El->SetBinContent(4,3,0.78);
  promptH2El->SetBinContent(5,3,0.84);
}

void setPrElClaraBDTHigh(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.53);
  promptH2El->SetBinContent(2,1,0.62);
  promptH2El->SetBinContent(3,1,0.74);
  promptH2El->SetBinContent(4,1,0.86);
  promptH2El->SetBinContent(5,1,0.90);
  //
  promptH2El->SetBinContent(1,2,0.27);
  promptH2El->SetBinContent(2,2,0.48);
  promptH2El->SetBinContent(3,2,0.66);
  promptH2El->SetBinContent(4,2,0.79);
  promptH2El->SetBinContent(5,2,0.87);
  //
  promptH2El->SetBinContent(1,3,0.30);
  promptH2El->SetBinContent(2,3,0.35);
  promptH2El->SetBinContent(3,3,0.59);
  promptH2El->SetBinContent(4,3,0.72);
  promptH2El->SetBinContent(5,3,0.79);
}


void setPrElClaraBDT2011B(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.54);
  promptH2El->SetBinContent(2,1,0.64);
  promptH2El->SetBinContent(3,1,0.75);
  promptH2El->SetBinContent(4,1,0.87);
  promptH2El->SetBinContent(5,1,0.90);
  //
  promptH2El->SetBinContent(1,2,0.41);
  promptH2El->SetBinContent(2,2,0.58);
  promptH2El->SetBinContent(3,2,0.69);
  promptH2El->SetBinContent(4,2,0.81);
  promptH2El->SetBinContent(5,2,0.85);
  //
  promptH2El->SetBinContent(1,3,0.31);
  promptH2El->SetBinContent(2,3,0.38);
  promptH2El->SetBinContent(3,3,0.62);
  promptH2El->SetBinContent(4,3,0.75);
  promptH2El->SetBinContent(5,3,0.82);
}


void setPrElEmaLP(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.61);
  promptH2El->SetBinContent(2,1,0.71);
  promptH2El->SetBinContent(3,1,0.82);
  promptH2El->SetBinContent(4,1,0.92);
  promptH2El->SetBinContent(5,1,0.96);
  //
  promptH2El->SetBinContent(1,2,0.43);
  promptH2El->SetBinContent(2,2,0.47);
  promptH2El->SetBinContent(3,2,0.67);
  promptH2El->SetBinContent(4,2,0.80);
  promptH2El->SetBinContent(5,2,0.90);
  //
  promptH2El->SetBinContent(1,3,0.42);
  promptH2El->SetBinContent(2,3,0.59);
  promptH2El->SetBinContent(3,3,0.75);
  promptH2El->SetBinContent(4,3,0.85);
  promptH2El->SetBinContent(5,3,0.92);
}

void setPrElEmaLH(TH2F* promptH2El){
  promptH2El->SetBinContent(1,1,0.28);
  promptH2El->SetBinContent(2,1,0.53);
  promptH2El->SetBinContent(3,1,0.66);
  promptH2El->SetBinContent(4,1,0.79);
  promptH2El->SetBinContent(5,1,0.85);
  //
  promptH2El->SetBinContent(1,2,0.33);
  promptH2El->SetBinContent(2,2,0.41);
  promptH2El->SetBinContent(3,2,0.54);
  promptH2El->SetBinContent(4,2,0.72);
  promptH2El->SetBinContent(5,2,0.82);
  //
  promptH2El->SetBinContent(1,3,0.26);
  promptH2El->SetBinContent(2,3,0.37);
  promptH2El->SetBinContent(3,3,0.46);
  promptH2El->SetBinContent(4,3,0.70);
  promptH2El->SetBinContent(5,3,0.80);

}


void setFrMuXieLP(TH2F* fakeH2Mu){
  fakeH2Mu->SetBinContent(1,1,0.192);
  fakeH2Mu->SetBinContent(2,1,0.136);
  fakeH2Mu->SetBinContent(3,1,0.254);
  fakeH2Mu->SetBinContent(4,1,0.242);
  fakeH2Mu->SetBinContent(5,1,0.237);
  //
  fakeH2Mu->SetBinContent(1,2,0.209);
  fakeH2Mu->SetBinContent(2,2,0.156);
  fakeH2Mu->SetBinContent(3,2,0.289);
  fakeH2Mu->SetBinContent(4,2,0.276);
  fakeH2Mu->SetBinContent(5,2,0.267);
  //
  fakeH2Mu->SetBinContent(1,3,0.229);
  fakeH2Mu->SetBinContent(2,3,0.173);
  fakeH2Mu->SetBinContent(3,3,0.256);
  fakeH2Mu->SetBinContent(4,3,0.253);
  fakeH2Mu->SetBinContent(5,3,0.253);
  //
  fakeH2Mu->SetBinContent(1,4,0.259);
  fakeH2Mu->SetBinContent(2,4,0.208);
  fakeH2Mu->SetBinContent(3,4,0.306);
  fakeH2Mu->SetBinContent(4,4,0.295);
  fakeH2Mu->SetBinContent(5,4,0.321);
}

void setFrMuAliciaLP(TH2F* fakeH2Mu){
  fakeH2Mu->SetBinContent(1,1,0.195);
  fakeH2Mu->SetBinContent(2,1,0.125);
  fakeH2Mu->SetBinContent(3,1,0.238);
  fakeH2Mu->SetBinContent(4,1,0.202);
  fakeH2Mu->SetBinContent(5,1,0.203);
  //
  fakeH2Mu->SetBinContent(1,2,0.222);
  fakeH2Mu->SetBinContent(2,2,0.150);
  fakeH2Mu->SetBinContent(3,2,0.277);
  fakeH2Mu->SetBinContent(4,2,0.249);
  fakeH2Mu->SetBinContent(5,2,0.224);
  //
  fakeH2Mu->SetBinContent(1,3,0.248);
  fakeH2Mu->SetBinContent(2,3,0.153);
  fakeH2Mu->SetBinContent(3,3,0.236);
  fakeH2Mu->SetBinContent(4,3,0.215);
  fakeH2Mu->SetBinContent(5,3,0.199);
  //
  fakeH2Mu->SetBinContent(1,4,0.291);
  fakeH2Mu->SetBinContent(2,4,0.190);
  fakeH2Mu->SetBinContent(3,4,0.265);
  fakeH2Mu->SetBinContent(4,4,0.241);
  fakeH2Mu->SetBinContent(5,4,0.261);
}

void setFrMuAliciaKink(TH2F* fakeH2Mu){
  /* as used for the note update  */
  fakeH2Mu->SetBinContent(1,1,0.197);
  fakeH2Mu->SetBinContent(2,1,0.128);
  fakeH2Mu->SetBinContent(3,1,0.239);
  fakeH2Mu->SetBinContent(4,1,0.197);
  fakeH2Mu->SetBinContent(5,1,0.179);
  //
  fakeH2Mu->SetBinContent(1,2,0.223);
  fakeH2Mu->SetBinContent(2,2,0.157);
  fakeH2Mu->SetBinContent(3,2,0.284);
  fakeH2Mu->SetBinContent(4,2,0.261);
  fakeH2Mu->SetBinContent(5,2,0.241);
  //
  fakeH2Mu->SetBinContent(1,3,0.252);
  fakeH2Mu->SetBinContent(2,3,0.165);
  fakeH2Mu->SetBinContent(3,3,0.250);
  fakeH2Mu->SetBinContent(4,3,0.250);
  fakeH2Mu->SetBinContent(5,3,0.257);
  //
  fakeH2Mu->SetBinContent(1,4,0.290);
  fakeH2Mu->SetBinContent(2,4,0.203);
  fakeH2Mu->SetBinContent(3,4,0.281);
  fakeH2Mu->SetBinContent(4,4,0.283);
  fakeH2Mu->SetBinContent(5,4,0.328);  

  // TO BE FILLED using 
  // /afs/cern.ch/user/c/calderon/public/WJetsSynch/NovAN/MuFR_all2011_LPcuts_AND_kink_jet30.txt 
  // (Jet 30 cut)
  /*
  fakeH2Mu->SetBinContent(1,1,0.126);
  fakeH2Mu->SetBinContent(2,1,0.083);
  fakeH2Mu->SetBinContent(3,1,0.179);
  fakeH2Mu->SetBinContent(4,1,0.189);
  fakeH2Mu->SetBinContent(5,1,0.175);
  //
  fakeH2Mu->SetBinContent(1,2,0.156);
  fakeH2Mu->SetBinContent(2,2,0.084);
  fakeH2Mu->SetBinContent(3,2,0.234);
  fakeH2Mu->SetBinContent(4,2,0.256);
  fakeH2Mu->SetBinContent(5,2,0.239);
  //
  fakeH2Mu->SetBinContent(1,3,0.169);
  fakeH2Mu->SetBinContent(2,3,0.087);
  fakeH2Mu->SetBinContent(3,3,0.194);
  fakeH2Mu->SetBinContent(4,3,0.242);
  fakeH2Mu->SetBinContent(5,3,0.257);
  //
  fakeH2Mu->SetBinContent(1,4,0.183);
  fakeH2Mu->SetBinContent(2,4,0.142);
  fakeH2Mu->SetBinContent(3,4,0.228);
  fakeH2Mu->SetBinContent(4,4,0.276);
  fakeH2Mu->SetBinContent(5,4,0.327);
  */
}


void setFrElChiaraLP(TH2F* fakeH2El){
  fakeH2El->SetBinContent(1,1,0.087);
  fakeH2El->SetBinContent(2,1,0.083);
  fakeH2El->SetBinContent(3,1,0.122);
  fakeH2El->SetBinContent(4,1,0.126);
  fakeH2El->SetBinContent(5,1,0.109);
  //
  fakeH2El->SetBinContent(1,2,0.057);
  fakeH2El->SetBinContent(2,2,0.067);
  fakeH2El->SetBinContent(3,2,0.116);
  fakeH2El->SetBinContent(4,2,0.099);
  fakeH2El->SetBinContent(5,2,0.114);
  //
  fakeH2El->SetBinContent(1,3,0.053);
  fakeH2El->SetBinContent(2,3,0.037);
  fakeH2El->SetBinContent(3,3,0.126);
  fakeH2El->SetBinContent(4,3,0.082);
  fakeH2El->SetBinContent(5,3,0.105);
  //
  fakeH2El->SetBinContent(1,4,0.064);
  fakeH2El->SetBinContent(2,4,0.050);
  fakeH2El->SetBinContent(3,4,0.120);
  fakeH2El->SetBinContent(4,4,0.115);
  fakeH2El->SetBinContent(5,4,0.130);
}

void setFrElChiaraLH(TH2F* fakeH2El){
  fakeH2El->SetBinContent(1,1,0.032);
  fakeH2El->SetBinContent(2,1,0.053);
  fakeH2El->SetBinContent(3,1,0.065);
  fakeH2El->SetBinContent(4,1,0.057);
  fakeH2El->SetBinContent(5,1,0.057);
  //
  fakeH2El->SetBinContent(1,2,0.017);
  fakeH2El->SetBinContent(2,2,0.035);
  fakeH2El->SetBinContent(3,2,0.039);
  fakeH2El->SetBinContent(4,2,0.044);
  fakeH2El->SetBinContent(5,2,0.029);
  //
  fakeH2El->SetBinContent(1,3,0.053);
  fakeH2El->SetBinContent(2,3,0.030);
  fakeH2El->SetBinContent(3,3,0.041);
  fakeH2El->SetBinContent(4,3,0.023);
  fakeH2El->SetBinContent(5,3,0.023);
  //
  fakeH2El->SetBinContent(1,4,0.040);
  fakeH2El->SetBinContent(2,4,0.031);
  fakeH2El->SetBinContent(3,4,0.035);
  fakeH2El->SetBinContent(4,4,0.034);
  fakeH2El->SetBinContent(5,4,0.026);
}

void setFrElXieBDT(TH2F* fakeH2El){
  fakeH2El->SetBinContent(1,1,0.072);
  fakeH2El->SetBinContent(2,1,0.077);
  fakeH2El->SetBinContent(3,1,0.086);
  fakeH2El->SetBinContent(4,1,0.077);
  fakeH2El->SetBinContent(5,1,0.079);
  //
  fakeH2El->SetBinContent(1,2,0.039);
  fakeH2El->SetBinContent(2,2,0.043);
  fakeH2El->SetBinContent(3,2,0.057);
  fakeH2El->SetBinContent(4,2,0.055);
  fakeH2El->SetBinContent(5,2,0.091);
  //
  fakeH2El->SetBinContent(1,3,0.021);
  fakeH2El->SetBinContent(2,3,0.015);
  fakeH2El->SetBinContent(3,3,0.050);
  fakeH2El->SetBinContent(4,3,0.036);
  fakeH2El->SetBinContent(5,3,0.077);
  //
  fakeH2El->SetBinContent(1,4,0.026);
  fakeH2El->SetBinContent(2,4,0.035);
  fakeH2El->SetBinContent(3,4,0.049);
  fakeH2El->SetBinContent(4,4,0.065);
  fakeH2El->SetBinContent(5,4,0.061);
}

void setFrElChiaraBDT2011B(TH2F* fakeH2El){
  fakeH2El->SetBinContent(1,1,0.065);
  fakeH2El->SetBinContent(2,1,0.051);
  fakeH2El->SetBinContent(3,1,0.078);
  fakeH2El->SetBinContent(4,1,0.069);
  fakeH2El->SetBinContent(5,1,0.049);
  //
  fakeH2El->SetBinContent(1,2,0.039);
  fakeH2El->SetBinContent(2,2,0.052);
  fakeH2El->SetBinContent(3,2,0.060);
  fakeH2El->SetBinContent(4,2,0.057);
  fakeH2El->SetBinContent(5,2,0.037);
  //
  fakeH2El->SetBinContent(1,3,0.016);
  fakeH2El->SetBinContent(2,3,0.016);
  fakeH2El->SetBinContent(3,3,0.054);
  fakeH2El->SetBinContent(4,3,0.046);
  fakeH2El->SetBinContent(5,3,0.060);
  //
  fakeH2El->SetBinContent(1,4,0.023);
  fakeH2El->SetBinContent(2,4,0.040);
  fakeH2El->SetBinContent(3,4,0.040);
  fakeH2El->SetBinContent(4,4,0.066);
  fakeH2El->SetBinContent(5,4,0.060);
  //
}

void setFrElChiaraBDT2011(TH2F* fakeH2El){
  /*
  fakeH2El->SetBinContent(1,1,0.065);
  fakeH2El->SetBinContent(2,1,0.051);
  fakeH2El->SetBinContent(3,1,0.078);
  fakeH2El->SetBinContent(4,1,0.069);
  fakeH2El->SetBinContent(5,1,0.049);
  //
  fakeH2El->SetBinContent(1,2,0.039);
  fakeH2El->SetBinContent(2,2,0.052);
  fakeH2El->SetBinContent(3,2,0.060);
  fakeH2El->SetBinContent(4,2,0.057);
  fakeH2El->SetBinContent(5,2,0.037);
  //
  fakeH2El->SetBinContent(1,3,0.016);
  fakeH2El->SetBinContent(2,3,0.016);
  fakeH2El->SetBinContent(3,3,0.054);
  fakeH2El->SetBinContent(4,3,0.046);
  fakeH2El->SetBinContent(5,3,0.060);
  //
  fakeH2El->SetBinContent(1,4,0.023);
  fakeH2El->SetBinContent(2,4,0.040);
  fakeH2El->SetBinContent(3,4,0.040);
  fakeH2El->SetBinContent(4,4,0.066);
  fakeH2El->SetBinContent(5,4,0.060);
  */

  //-- REFERENCE
  fakeH2El->SetBinContent(1,1,0.066);
  fakeH2El->SetBinContent(2,1,0.057);
  fakeH2El->SetBinContent(3,1,0.085);
  fakeH2El->SetBinContent(4,1,0.085);
  fakeH2El->SetBinContent(5,1,0.085);
  //
  fakeH2El->SetBinContent(1,2,0.040);
  fakeH2El->SetBinContent(2,2,0.054);
  fakeH2El->SetBinContent(3,2,0.064);
  fakeH2El->SetBinContent(4,2,0.067);
  fakeH2El->SetBinContent(5,2,0.059);
  //
  fakeH2El->SetBinContent(1,3,0.016);
  fakeH2El->SetBinContent(2,3,0.018);
  fakeH2El->SetBinContent(3,3,0.057);
  fakeH2El->SetBinContent(4,3,0.052);
  fakeH2El->SetBinContent(5,3,0.072);
  //
  fakeH2El->SetBinContent(1,4,0.023);
  fakeH2El->SetBinContent(2,4,0.041);
  fakeH2El->SetBinContent(3,4,0.046);
  fakeH2El->SetBinContent(4,4,0.070);
  fakeH2El->SetBinContent(5,4,0.068);

  //to be updated using 
  // https://crovelli.web.cern.ch/crovelli/Analysis/HWW/fakeRate/lastFR_2011AandB_bdt_withEWKsub_halfLumi.txt
  // (Et 15)
  /*
  fakeH2El->SetBinContent(1,1,0.109);
  fakeH2El->SetBinContent(2,1,0.103);
  fakeH2El->SetBinContent(3,1,0.123);
  fakeH2El->SetBinContent(4,1,0.099);
  fakeH2El->SetBinContent(5,1,0.102);
  //
  fakeH2El->SetBinContent(1,2,0.056);
  fakeH2El->SetBinContent(2,2,0.074);
  fakeH2El->SetBinContent(3,2,0.093);
  fakeH2El->SetBinContent(4,2,0.070);
  fakeH2El->SetBinContent(5,2,0.058);
  //
  fakeH2El->SetBinContent(1,3,0.028);
  fakeH2El->SetBinContent(2,3,0.028);
  fakeH2El->SetBinContent(3,3,0.085);
  fakeH2El->SetBinContent(4,3,0.059);
  fakeH2El->SetBinContent(5,3,0.071);
  //
  fakeH2El->SetBinContent(1,4,0.035);
  fakeH2El->SetBinContent(2,4,0.048);
  fakeH2El->SetBinContent(3,4,0.069);
  fakeH2El->SetBinContent(4,4,0.080);
  fakeH2El->SetBinContent(5,4,0.073);
  */

  //
  //
  // Chiara FR with original mT cut, no contamination subtraction
  /*
  fakeH2El->SetBinContent(1,1,0.093);
  fakeH2El->SetBinContent(2,1,0.056);
  fakeH2El->SetBinContent(3,1,0.081);
  fakeH2El->SetBinContent(4,1,0.072);
  fakeH2El->SetBinContent(5,1,0.076);
  //
  fakeH2El->SetBinContent(1,2,0.040);
  fakeH2El->SetBinContent(2,2,0.040);
  fakeH2El->SetBinContent(3,2,0.066);
  fakeH2El->SetBinContent(4,2,0.053);
  fakeH2El->SetBinContent(5,2,0.037);
  //
  fakeH2El->SetBinContent(1,3,0.019);
  fakeH2El->SetBinContent(2,3,0.018);
  fakeH2El->SetBinContent(3,3,0.052);
  fakeH2El->SetBinContent(4,3,0.038);
  fakeH2El->SetBinContent(5,3,0.049);
  //
  fakeH2El->SetBinContent(1,4,0.030);
  fakeH2El->SetBinContent(2,4,0.028);
  fakeH2El->SetBinContent(3,4,0.044);
  fakeH2El->SetBinContent(4,4,0.050);
  fakeH2El->SetBinContent(5,4,0.051);
  */

}



void setFrElXieLP(TH2F* fakeH2El){
  fakeH2El->SetBinContent(1,1,0.133);
  fakeH2El->SetBinContent(2,1,0.091);
  fakeH2El->SetBinContent(3,1,0.113);
  fakeH2El->SetBinContent(4,1,0.102);
  fakeH2El->SetBinContent(5,1,0.123);
  //
  fakeH2El->SetBinContent(1,2,0.075);
  fakeH2El->SetBinContent(2,2,0.063);
  fakeH2El->SetBinContent(3,2,0.115);
  fakeH2El->SetBinContent(4,2,0.089);
  fakeH2El->SetBinContent(5,2,0.123);
  //
  fakeH2El->SetBinContent(1,3,0.033);
  fakeH2El->SetBinContent(2,3,0.039);
  fakeH2El->SetBinContent(3,3,0.114);
  fakeH2El->SetBinContent(4,3,0.090);
  fakeH2El->SetBinContent(5,3,0.099);
  //
  fakeH2El->SetBinContent(1,4,0.057);
  fakeH2El->SetBinContent(2,4,0.057);
  fakeH2El->SetBinContent(3,4,0.112);
  fakeH2El->SetBinContent(4,4,0.120);
  fakeH2El->SetBinContent(5,4,0.120);

}

  /*
  //dummy mu PR
  promptH2Mu->SetBinContent(1,1,1.);
  promptH2Mu->SetBinContent(2,1,1.);
  promptH2Mu->SetBinContent(3,1,1.);
  promptH2Mu->SetBinContent(4,1,1.);
  promptH2Mu->SetBinContent(5,1,1.);
  //
  promptH2Mu->SetBinContent(1,2,1.);
  promptH2Mu->SetBinContent(2,2,1.);
  promptH2Mu->SetBinContent(3,2,1.);
  promptH2Mu->SetBinContent(4,2,1.);
  promptH2Mu->SetBinContent(5,2,1.);
  // ----

  //dummy el PR
  promptH2El->SetBinContent(1,1,1.);
  promptH2El->SetBinContent(2,1,1.);
  promptH2El->SetBinContent(3,1,1.);
  promptH2El->SetBinContent(4,1,1.);
  promptH2El->SetBinContent(5,1,1.);
  //
  promptH2El->SetBinContent(1,2,1.);
  promptH2El->SetBinContent(2,2,1.);
  promptH2El->SetBinContent(3,2,1.);
  promptH2El->SetBinContent(4,2,1.);
  promptH2El->SetBinContent(5,2,1.);
  //
  promptH2El->SetBinContent(1,3,1.);
  promptH2El->SetBinContent(2,3,1.);
  promptH2El->SetBinContent(3,3,1.);
  promptH2El->SetBinContent(4,3,1.);
  promptH2El->SetBinContent(5,3,1.);
  */
