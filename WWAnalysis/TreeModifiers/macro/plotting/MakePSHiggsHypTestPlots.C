#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Math/DistFunc.h"


// void PlotHypSeparation(string filename ) {

//   TFile * fitFile = TFile::Open(filename.c_str());
//   TTree * tree = (TTree*) fitFile->Get("q");

//   TH1D * SMHiggsTestStat = new TH1D("SMHiggsTestStat", "; Test Statistic = -2 #times ln(L_{1}/L_{2}); Number of Toy Experiments", 80,-20,20);
//   TH1D * PSHiggsTestStat = new TH1D("PSHiggsTestStat", "; Test Statistic = -2 #times ln(L_{1}/L_{2}); Number of Toy Experiments", 80,-20,20);

//   tree->Draw("q>>SMHiggsTestStat","type == 1");
//   tree->Draw("q>>PSHiggsTestStat","type == -1");

//   cout << PSHiggsTestStat->GetMean() << " " << SMHiggsTestStat->GetMean() << "\n";
//   cout << SMHiggsTestStat->Integral(0,SMHiggsTestStat->GetXaxis()->FindFixBin(PSHiggsTestStat->GetMean())) /  SMHiggsTestStat->Integral() << endl;

//   double avgPValue = SMHiggsTestStat->Integral(0,SMHiggsTestStat->GetXaxis()->FindFixBin(PSHiggsTestStat->GetMean())) /  SMHiggsTestStat->Integral();

//   cout << "Average Separation: " << RooStats::PValueToSignificance(avgPValue) << endl;



//   TCanvas* cv = new TCanvas("cv","cv",610,0,600,600);
//   SMHiggsTestStat->SetLineColor(kMagenta+1);
//   PSHiggsTestStat->SetLineColor(kBlue);

//   SMHiggsTestStat->GetYaxis()->SetTitleOffset(1.3);
//   SMHiggsTestStat->GetXaxis()->SetTitleOffset(1.05);

//   SMHiggsTestStat->Draw("hist");
//   PSHiggsTestStat->Draw("hist,same");

//   cv->SaveAs("PSHiggsHypSeparation.gif");

// }



void PlotHypSeparation(string filename ) {

  const float lumi7TeV=5.05;
  const float lumi8TeV=12.2;

  TFile *fq=new TFile(filename.c_str(),"READ");
  TTree *t=(TTree*)fq->Get("q");

  float q,m,w;
  int type;
  t->SetBranchAddress("q",&q);
  t->SetBranchAddress("mh",&m);
  t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("type",&type);

  TH1F *hSM=new TH1F("hSM;S = -2 #times ln(L_{1}/L_{2});Number of Toys","",8000,-40,40);
  TH1F *hPS=new TH1F("hPS;S = -2 #times ln(L_{1}/L_{2});Number of Toys","",8000,-40,40);
  TH1F *hObs=new TH1F("hObserved","",100,-20,20);
  cout<<"Start to lopp on tree in file "<<filename<<endl;

  std::vector<float> v_SM, v_PS,v_Obs;

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if(i==0)cout<<"MASS in the TREE = "<<m<<endl<<endl;

    if(type<0){ //SM hypothesis 
      hSM->Fill(-q);
      v_SM.push_back(-q);
    }
    else if(type>0){//ALT hypothesis
      hPS->Fill(-q);
      v_PS.push_back(-q);
    }
    else{
      hObs->Fill(q);
      v_Obs.push_back(q);
    }
  }//end loop on tree entries
  cout<<"Finished to loop, sorting vectors "<<v_SM.size()<<" "<<v_PS.size()<<" "<<v_Obs.size()<<endl;
  sort(v_SM.begin(),v_SM.end());//sort in ascending order
  sort(v_PS.begin(),v_PS.end()); 
  sort(v_Obs.begin(),v_Obs.end());
  int ntoysSM= hSM->GetEntries();
  int ntoysPS= hPS->GetEntries();

  //we assume that SM is on the right and PS on the left of zero
  if(v_PS.at(0)>v_SM.at(ntoysSM-1)){
    cout<<"Swapped distributions !!! The alternative model shouldstay on the negative side of the significance."<<endl;
    cout<<"Please edit the code and change the sign of q when filling histos and vectors in the loop on tree entries"<<endl;
    return ;
  }
  float medianSM=v_SM.at(int(ntoysSM/2));
  float medianPS=v_PS.at(int(ntoysPS/2));
  cout<<"Toys generated "<<ntoysSM<<"\t"<<ntoysPS<<endl;
  cout<<"Mean of SM/PS hypothesis: "<<hSM->GetMean()<<"\t"<<hPS->GetMean()<<endl;
  cout<<"RMS  of SM/PS hypothesis: "<<hSM->GetRMS()<<"\t"<<hPS->GetRMS()<<endl;
  cout<<"Median of SM/PS hypothesis: "<<medianSM<<"\t"<<medianPS<<endl;

  const float step=0.05;
  float coverage=0.0;
  float diff=10.0;
  float cut=v_PS.at(0)-step;
  float crosspoint=-99.0;
  int startSM=ntoysSM-1, startPS=0;
  cout<<"Starting to loop with cut at "<<cut<<endl;



  cout<<"Finished loop on vector elements, min is "<<diff<<" cut is at "<<cut<<endl;
  cout<<"q value where SM and ALT distributions have same area on opposite sides: "<<crosspoint<<endl;
  cout<<"Coverage "<<coverage<<endl;
  float separation=2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation: "<<separation<<endl<<endl<<endl;

  float integralSM=hSM->Integral();
  float integralPS=hPS->Integral();
 
  float tailSM=hSM->Integral(1,hSM->FindBin(medianPS))/integralSM;
  float tailPS=hPS->Integral(hPS->FindBin(medianSM),hPS->GetNbinsX())/integralPS;
  cout<<"Tail prob SM: "<<tailSM<<"  ("<<ROOT::Math::normal_quantile_c(tailSM,1.0) <<" sigma)"<<endl;
  cout<<"Tail prob PS: "<<tailPS<<"  ("<<ROOT::Math::normal_quantile_c(tailPS,1.0) <<" sigma)"<<endl;

  diff=10.0;
  coverage=0.0;
  for(int i=1;i<hSM->GetNbinsX();i++){
    
    float fracSM=hSM->Integral(1,i) / integralSM;
    float fracPS=hPS->Integral(i,hPS->GetNbinsX()) / integralPS;
    if(fabs(fracSM-fracPS)<diff){
      diff=fabs(fracSM-fracPS);
      coverage=(fracSM+fracPS)/2.0;
    }

  }

  float sepH= 2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation from histograms = "<<sepH<<" with coverage "<<coverage<<endl;

  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c1","c1",800,800);
  c1->cd();
  hSM->Rebin(50);
  hPS->Rebin(50);
  float maxhSM=hSM->GetBinContent(hSM->GetMaximumBin());
  float maxhPS=hPS->GetBinContent(hPS->GetMaximumBin());
  if(maxhPS>maxhSM){
    hSM->SetMaximum(maxhPS*1.15);
    hPS->SetMaximum(maxhPS*1.15);
  }
  else{
    hSM->SetMaximum(maxhSM*1.15);
    hPS->SetMaximum(maxhSM*1.15);
  }
  hSM->SetXTitle("S = -2 #times ln(L_{1}/L_{2})");
  hSM->SetYTitle("Number of Pseudo-experiments");
  hPS->SetXTitle("S = -2 #times ln(L_{1}/L_{2})");
  hPS->SetYTitle("Number of Pseudo-experiments");
  hSM->SetLineColor(kMagenta-3);
  hSM->SetFillColor(kMagenta-3);
  hSM->SetLineWidth(2);
  hSM->SetFillStyle(3605);
  hPS->SetLineColor(kBlue+1);
  hPS->SetFillColor(kBlue+1);
  hPS->SetLineWidth(2);
  hPS->SetFillStyle(3695);
  hObs->SetLineColor(kGreen+3);
  hObs->SetLineWidth(2);


  hSM->GetXaxis()->SetRangeUser(-15,15);

  hSM->Draw();
  hPS->Draw("sames");
  //  hObs->Draw("sames");

  TLegend *leg = new TLegend(0.2,0.6,0.45,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hSM,"  SM, 0+","f");
  leg->AddEntry(hPS,"  PS, 0-","f");
  leg->Draw();


  TPaveText pt(0.16,0.95,0.45,0.99,"NDC");
  pt.SetFillColor(0);
  pt.AddText("CMS Expected");
  pt.SetBorderSize(0);
  TPaveText pt2(0.55,0.95,0.99,0.99,"NDC");
  pt2.SetFillColor(0);
  pt2.AddText(Form(" #sqrt{s} = 7 TeV, L = %.3f fb^{-1}; #sqrt{s} = 8 TeV, L = %.3f fb^{-1}",lumi7TeV,lumi8TeV));
  pt2.SetBorderSize(0);

  TPaveText pt3(0.65,0.75,0.85,0.90,"NDC");
  pt3.SetFillColor(0);
  pt3.AddText(Form("Average Separation: %.2f #sigma",fabs(sepH)));
  pt3.SetTextSize(0.03);
  pt3.SetBorderSize(0);

  pt.Draw();
  pt2.Draw();
  pt3.Draw();

  c1->SaveAs("PSHiggsHypSeparation.gif");

}




void MakePSHiggsHypTestPlots() {

  PlotHypSeparation("/afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/v5/qmu.FixedMu.root");
//   PlotHypSeparation("/afs/cern.ch/work/s/sixie/public/HZZ4l/cards/HCP/pseudoscalar/CJLST/qmu.FixedMu.root");

}
