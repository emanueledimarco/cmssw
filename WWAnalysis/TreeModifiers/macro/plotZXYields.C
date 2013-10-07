#include <iostream>
#include <sstream>

#include "TMath.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1F.h"

#include "RooPlot.h"
#include "RooLandau.h"
#include "RooRealVar.h"

using namespace std;
using namespace RooFit;

double wavg(double x, double y, double ex, double ey) {
  return (x/ex/ex + y/ey/ey)/(1/ex/ex + 1/ey/ey);
}

pair<double,double> werror(double x, double y, double ex, double ey) {
  double max = TMath::Max(x+ex,y+ey);
  double min = TMath::Min(x-ex,y-ey);

  double val = wavg(x,y,ex,ey);
  double eu=fabs(max-val);
  double ed=fabs(val-min);
  
  return std::make_pair<double,double>(ed,eu);
}

void plotZXYields(bool do7TeV) {

  // SS method
  double ssyield_7TeV[3] = {0.92, 2.52, 3.10};
  double ssyielde_7TeV[3] = {0.20, 0.30, 0.48};
  double ssyield_8TeV[3] = {3.75, 9.52, 12.8};
  double ssyielde_8TeV[3] = {0.52, 0.64, 0.91};

  // OS method
  double osyield_7TeV[3] = {1.18, 0.47, 0.89};
  double osyielde_7TeV[3] = {0.42, 0.21, 0.58};
  double osyield_8TeV[3] = {1.52, 3.98, 5.18};
  double osyielde_8TeV[3] = {0.71, 0.61, 1.02};

  // average
  double  yield_7TeV[3], yieldeup_7TeV[3], yieldedn_7TeV[3], yield_8TeV[3], yieldeup_8TeV[3], yieldedn_8TeV[3];
  
  for(int i=0;i<3;++i) {
    yield_7TeV[i] = wavg(ssyield_7TeV[i],osyield_7TeV[i],ssyielde_7TeV[i],osyielde_7TeV[i]);
    pair<double,double> err=werror(ssyield_7TeV[i],osyield_7TeV[i],ssyielde_7TeV[i],osyielde_7TeV[i]);
    yieldedn_7TeV[i] = err.first;
    yieldeup_7TeV[i] = err.second;

    yield_8TeV[i] = wavg(ssyield_8TeV[i],osyield_8TeV[i],ssyielde_8TeV[i],osyielde_8TeV[i]);
    pair<double,double> error=werror(ssyield_8TeV[i],osyield_8TeV[i],ssyielde_8TeV[i],osyielde_8TeV[i]);
    yieldedn_8TeV[i] = error.first;
    yieldeup_8TeV[i] = error.second;

    cout << "channel = " << i << "\tyield(7 TeV) = " << yield_7TeV[i] << " + " << yieldeup_7TeV[i] << " - " << yieldedn_7TeV[i] << endl;
    cout << "channel = " << i << "\tyield(8 TeV) = " << yield_8TeV[i] << " + " << yieldeup_8TeV[i] << " - " << yieldedn_8TeV[i] << endl;
  }

  gStyle->SetLineWidth(1.);
  gStyle->SetMarkerStyle(8);

  TGraphAsymmErrors *grSS = new TGraphAsymmErrors(3);
  TGraphAsymmErrors *grOS = new TGraphAsymmErrors(3);
  TGraphAsymmErrors *grAvg = new TGraphAsymmErrors(3);

  TGraph *grDummy = new TGraph(2);
  grDummy->SetPoint(0,-1.,6);
  grDummy->SetPoint(1,3,6);
  grDummy->SetMarkerColor(0);
  grDummy->GetXaxis()->SetBinLabel(0,"");
  grDummy->SetTitle("");

  for(int ch=0;ch<3;++ch) {
    if(do7TeV) {
      grSS->SetPoint(ch,ch-0.1,ssyield_7TeV[ch]);
      grSS->SetPointError(ch,0,0,ssyielde_7TeV[ch],ssyielde_7TeV[ch]);

      grOS->SetPoint(ch,ch+0.1,osyield_7TeV[ch]);
      grOS->SetPointError(ch,0,0,osyielde_7TeV[ch],osyielde_7TeV[ch]);

      grAvg->SetPoint(ch,ch,yield_7TeV[ch]);
      grAvg->SetPointError(ch,0,0,yieldedn_7TeV[ch],yieldeup_7TeV[ch]);
    } else {
      grSS->SetPoint(ch,ch-0.1,ssyield_8TeV[ch]);
      grSS->SetPointError(ch,0,0,ssyielde_8TeV[ch],ssyielde_8TeV[ch]);

      grOS->SetPoint(ch,ch+0.1,osyield_8TeV[ch]);
      grOS->SetPointError(ch,0,0,osyielde_8TeV[ch],osyielde_8TeV[ch]);

      grAvg->SetPoint(ch,ch,yield_8TeV[ch]);
      grAvg->SetPointError(ch,0,0,yieldedn_8TeV[ch],yieldeup_8TeV[ch]);
    }
  }

  grSS->SetMarkerColor(kBlue+2);
  grOS->SetMarkerColor(kCyan+2);
  grAvg->SetMarkerColor(kViolet-2);

  grSS->SetLineColor(kBlue+2);
  grOS->SetLineColor(kCyan+2);
  grAvg->SetLineColor(kViolet-2);

  TLegend* leg = new TLegend(0.2, 0.6, 0.5, 0.8);
  leg->SetBorderSize(     0);
  leg->SetFillColor (  4000);
  leg->SetTextAlign (    12);
  leg->SetTextFont  (    42);
  leg->SetTextSize  (0.05);
  leg->AddEntry(grOS, "OS method","pl");
  leg->AddEntry(grSS, "SS method","pl");
  leg->AddEntry(grAvg, "Average","pl");

  TLatex* CP = new TLatex(-1.,do7TeV ? 10.5 : 20.5, do7TeV ? "CMS Preliminary                                                  #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                                              #sqrt{s} = 8 TeV, L = 12.2 fb^{-1}");
  CP->SetTextSize(0.035);

  TLatex* xax0 = new TLatex(0-.1,-1, "4#mu");
  TLatex* xax1 = new TLatex(1-.1,-1, "4e");
  TLatex* xax2 = new TLatex(2-.1,-1, "2e2#mu");

  TCanvas *c1 = new TCanvas("c1","c1");
  grDummy->GetYaxis()->SetRangeUser(0,((do7TeV) ? 10 : 20));
  grDummy->Draw("ap");
  grSS->Draw("PE");
  grOS->Draw("PE");
  grAvg->Draw("PE");
  CP->Draw();
  xax0->Draw();
  xax1->Draw();
  xax2->Draw();
  leg->Draw();

  stringstream fss;
  if(do7TeV) fss << "ZXYields-7TeV.pdf";
  else fss << "ZXYields-8TeV.pdf";
  c1->SaveAs(fss.str().c_str());

}


void plotZXShapes(bool do7TeV, bool zoom=false) {

  // SS method
  double ssme_7TeV[3] = { 140.4, 154.4, 152.5 }; 
  double sssi_7TeV[3] = { 18.0,   22.1,  22.1 };
  double ssme_8TeV[3] = { 147.8, 156.9, 153.2 }; 
  double sssi_8TeV[3] = {  29.1,  22.0,  22.5 };

  // OS method
  double osme_7TeV[3] = { 133.9, 112.1, 130.0 }; 
  double ossi_7TeV[3] = {  20.1,  23.8,  18.8 };
  double osme_8TeV[3] = { 145.7, 150.1, 141.9 }; 
  double ossi_8TeV[3] = {  24.5,  28.1,  27.2 };

  // average
  double me_7TeV[3], si_7TeV[3], me_8TeV[3], si_8TeV[3];
  for(int i=0;i<3;++i) {
    me_7TeV[i] = wavg(ssme_7TeV[i],osme_7TeV[i],1,1);
    si_7TeV[i] = wavg(sssi_7TeV[i],ossi_7TeV[i],1,1);

    me_8TeV[i] = wavg(ssme_8TeV[i],osme_8TeV[i],1,1);
    si_8TeV[i] = wavg(sssi_8TeV[i],ossi_8TeV[i],1,1);

    if(do7TeV) {
      cout << "==>\t7TeV. Channel = " << i << endl;
      cout << "\t zxme = " << me_7TeV[i] << endl
	   << "\t zxsi = " << si_7TeV[i] << endl;
    } else {
      cout << "==>\t8TeV. Channel = " << i << endl;
      cout << "\t zxme = " << me_8TeV[i] << endl
	   << "\t zxsi = " << si_8TeV[i] << endl;
    }
  }

  float maxmass = zoom ? 160 : 600;
  RooRealVar mass("mass","m_{4l}",100,maxmass,"GeV/c^{2}");

  std::vector<RooLandau*> sspdfs,ospdfs,pdfs;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);

  TLegend* legend = new TLegend(0.60, 0.70, 0.80, 0.85);
  
  legend->SetBorderSize(     0);
  legend->SetFillColor (  4000);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);

  TH1F *grSS = new TH1F("grSS","grSS",1,0,1);
  TH1F *grOS = new TH1F("grOS","grOS",1,0,1);
  TH1F *grAvg = new TH1F("grAvg","grAvg",1,0,1);

  grSS->SetLineColor(kBlue+2);
  grOS->SetLineColor(kCyan+2);
  grAvg->SetLineColor(kViolet-2);
  
  legend->AddEntry(grOS, "OS method","pl");
  legend->AddEntry(grSS, "SS method","pl");
  legend->AddEntry(grAvg, "Average","pl");

  for(int i=0;i<3;++i) {

    stringstream sspdf,ospdf,pdf,sme,sssme,ossme,ssi,ssssi,osssi;
    sspdf << "sspdf_cha" << i; ospdf << "ospdf_cha" << i; pdf << "pdf_cha" << i; 
    sme << "me_cha" << i;    sssme << "ssme_cha" << i;    ossme << "osme_cha" << i;
    ssi << "si_cha" << i;    ssssi << "sssi_cha" << i;    osssi << "ossi_cha" << i;
    RooRealVar ssme(sssme.str().c_str(),sssme.str().c_str(),120,100,200);
    RooRealVar sssi(ssssi.str().c_str(),ssssi.str().c_str(),25,10,50);

    RooRealVar osme(ossme.str().c_str(),ossme.str().c_str(),120,100,200);
    RooRealVar ossi(osssi.str().c_str(),osssi.str().c_str(),25,10,50);

    RooRealVar me(sme.str().c_str(),sme.str().c_str(),120,100,200);
    RooRealVar si(ssi.str().c_str(),ssi.str().c_str(),25,10,50);

    if(do7TeV) {
      ssme.setVal(ssme_7TeV[i]);
      sssi.setVal(sssi_7TeV[i]);
      RooLandau *sslnd = new RooLandau(sspdf.str().c_str(),sspdf.str().c_str(),mass,ssme,sssi);      
      sspdfs.push_back(sslnd);      

      osme.setVal(osme_7TeV[i]);
      ossi.setVal(ossi_7TeV[i]);
      RooLandau *oslnd = new RooLandau(ospdf.str().c_str(),ospdf.str().c_str(),mass,osme,ossi);      
      ospdfs.push_back(oslnd);      

      me.setVal(me_7TeV[i]);
      si.setVal(si_7TeV[i]);
      RooLandau *lnd = new RooLandau(pdf.str().c_str(),pdf.str().c_str(),mass,me,si);      
      pdfs.push_back(lnd);      
    } else {
      ssme.setVal(ssme_8TeV[i]);
      sssi.setVal(sssi_8TeV[i]);
      RooLandau *sslnd = new RooLandau(sspdf.str().c_str(),sspdf.str().c_str(),mass,ssme,sssi);      
      sspdfs.push_back(sslnd);      

      osme.setVal(osme_8TeV[i]);
      ossi.setVal(ossi_8TeV[i]);
      RooLandau *oslnd = new RooLandau(ospdf.str().c_str(),ospdf.str().c_str(),mass,osme,ossi);      
      ospdfs.push_back(oslnd);      

      me.setVal(me_8TeV[i]);
      si.setVal(si_8TeV[i]);
      RooLandau *lnd = new RooLandau(pdf.str().c_str(),pdf.str().c_str(),mass,me,si);      
      pdfs.push_back(lnd);      
    }

    RooPlot *frame = mass.frame(50);
    frame->SetTitle("");

    sspdfs.back()->plotOn(frame,LineColor(kBlue+2));
    ospdfs.back()->plotOn(frame,LineColor(kCyan+2));
    pdfs.back()->plotOn(frame,LineColor(kViolet-2));

    stringstream cname;
    if(do7TeV) cname << "zxpdf-cha" << i << "-7TeV" << (zoom ? "-zoom" : "") << ".pdf";
    else cname << "zxpdf-cha" << i << "-8TeV" << (zoom ? "-zoom" : "") << ".pdf";

    c1->cd();
    frame->Draw();
    legend->Draw();
    c1->SaveAs(cname.str().c_str());
  }


}
