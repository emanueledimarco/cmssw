#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>
#include <RooHistFunc.h>

#include "scales2.h"

using namespace RooFit;

std::vector<TH1*> createHists(int ch, bool is7) {

    init(is7);

    std::string chstr;
    if (ch == 0) chstr = "4mu";
    if (ch == 1) chstr = "4e";
    if (ch == 2) chstr = "2e2mu";

    std::string tevstr = is7 ? "_7TeV" : "_8TeV"; 

    float mass = 710.;

    TH1F* histxsecbrggh;
    TH1F* histxsecbrvbf;
    TH1F* histxsecbrwhi;
    TH1F* histxsecbrzhi;
    TH1F* histxsecbrtth;


    if (mass >= 110 && mass < 160) {
        Float_t gghxsecbry[50];
        Float_t vbfxsecbry[50];
        Float_t whixsecbry[50];
        Float_t zhixsecbry[50];
        Float_t tthxsecbry[50];
        
        
        for (float i = 110.; i < 160.; i = i+1.) {
            gghxsecbry[(int)(i-110.)] = getXsecggHByChannel(i, ch); 
            vbfxsecbry[(int)(i-110.)] = getXsecVBFByChannel(i, ch); 
            whixsecbry[(int)(i-110.)] = getXsecWHiByChannel(i, ch); 
            zhixsecbry[(int)(i-110.)] = getXsecZHiByChannel(i, ch); 
            tthxsecbry[(int)(i-110.)] = getXsecttHByChannel(i, ch); 
        }
        
        histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
        histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
        histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
        histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
        histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
        
        for (int i = 1; i <= 50; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
        for (int i = 1; i <= 50; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
        for (int i = 1; i <= 50; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
        for (int i = 1; i <= 50; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
        for (int i = 1; i <= 50; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
    }

    else if (mass >= 160 && mass < 290) {
        Float_t gghxsecbry[65];
        Float_t vbfxsecbry[65];
        Float_t whixsecbry[65];
        Float_t zhixsecbry[65];
        Float_t tthxsecbry[65];


        for (float i = 160.; i < 290.; i = i+2.) {
            gghxsecbry[((int)(i-160.))/2] = getXsecggHByChannel(i, ch);
            vbfxsecbry[((int)(i-160.))/2] = getXsecVBFByChannel(i, ch);
            whixsecbry[((int)(i-160.))/2] = getXsecWHiByChannel(i, ch);
            zhixsecbry[((int)(i-160.))/2] = getXsecZHiByChannel(i, ch);
            tthxsecbry[((int)(i-160.))/2] = getXsecttHByChannel(i, ch);
        }

        histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
        histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
        histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
        histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
        histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 65, 160., 290.);

        for (int i = 1; i <= 65; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
        for (int i = 1; i <= 65; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
        for (int i = 1; i <= 65; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
        for (int i = 1; i <= 65; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
        for (int i = 1; i <= 65; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
    }

    else if (mass >= 290 && mass < 350) {
        Float_t gghxsecbry[12];
        Float_t vbfxsecbry[12];
        Float_t whixsecbry[12];
        Float_t zhixsecbry[12];
        Float_t tthxsecbry[12];


        for (float i = 290.; i < 350.; i = i+5.) {
            gghxsecbry[((int)(i-290.))/5] = getXsecggHByChannel(i, ch);
            vbfxsecbry[((int)(i-290.))/5] = getXsecVBFByChannel(i, ch);
            whixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecWHiByChannel(i, ch) : 0.0;
            zhixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecZHiByChannel(i, ch) : 0.0;
            tthxsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecttHByChannel(i, ch) : 0.0;
        }

        histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 12, 290., 350);
        histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 12, 290., 350);
        histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 12, 290., 350);
        histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 12, 290., 350);
        histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 12, 290., 350);

        for (int i = 1; i <= 12; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
        for (int i = 1; i <= 12; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
        for (int i = 1; i <= 12; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
        for (int i = 1; i <= 12; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
        for (int i = 1; i <= 12; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
    }

    else if (mass >= 350 && mass < 400) {
        Float_t gghxsecbry[5];
        Float_t vbfxsecbry[5];
        Float_t whixsecbry[5];
        Float_t zhixsecbry[5];
        Float_t tthxsecbry[5];


        for (float i = 350.; i < 400.; i = i+10.) {
            gghxsecbry[((int)(i-350.))/10] = getXsecggHByChannel(i, ch);
            vbfxsecbry[((int)(i-350.))/10] = getXsecVBFByChannel(i, ch);
            whixsecbry[((int)(i-350.))/10] = 0.0;
            zhixsecbry[((int)(i-350.))/10] = 0.0;
            tthxsecbry[((int)(i-350.))/10] = 0.0;
        }

        histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
        histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
        histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
        histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
        histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 5, 350., 400.);

        for (int i = 1; i <= 5; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
        for (int i = 1; i <= 5; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
        for (int i = 1; i <= 5; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
        for (int i = 1; i <= 5; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
        for (int i = 1; i <= 5; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
    }

    else  {
        Float_t gghxsecbry[30];
        Float_t vbfxsecbry[30];
        Float_t whixsecbry[30];
        Float_t zhixsecbry[30];
        Float_t tthxsecbry[30];


        for (float i = 400.; i < 1000.; i = i+20.) {
            gghxsecbry[((int)(i-400.))/20] = getXsecggHByChannel(i, ch);
            vbfxsecbry[((int)(i-400.))/20] = getXsecVBFByChannel(i, ch);
            whixsecbry[((int)(i-400.))/20] = 0.;
            zhixsecbry[((int)(i-400.))/20] = 0.;
            tthxsecbry[((int)(i-400.))/20] = 0.;
        }

        histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
        histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
        histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
        histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
        histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);

        for (int i = 1; i <= 30; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
        for (int i = 1; i <= 30; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
        for (int i = 1; i <= 30; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
        for (int i = 1; i <= 30; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
        for (int i = 1; i <= 30; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
    }


    std::vector<TH1*> histvector;
    histvector.push_back(histxsecbrggh);
    histvector.push_back(histxsecbrvbf);
    histvector.push_back(histxsecbrwhi);
    histvector.push_back(histxsecbrzhi);
    histvector.push_back(histxsecbrtth);
    return histvector;

}



void makeXsecHistos() {

    std::vector<TH1*> histos;
    std::vector<TH1*> histos1;
    std::vector<TH1*> histos2;
    std::vector<TH1*> histos3;
    std::vector<TH1*> histos4;
    std::vector<TH1*> histos5;
    std::vector<TH1*> histos6;

    histos1 = createHists(0, true);
    histos2 = createHists(1, true);
    histos3 = createHists(2, true);
    histos4 = createHists(0, false);
    histos5 = createHists(1, false);
    histos6 = createHists(2, false);

    for (unsigned int i = 0; i < histos1.size(); i++) histos.push_back(histos1[i]);
    for (unsigned int i = 0; i < histos2.size(); i++) histos.push_back(histos2[i]);
    for (unsigned int i = 0; i < histos3.size(); i++) histos.push_back(histos3[i]);
    for (unsigned int i = 0; i < histos4.size(); i++) histos.push_back(histos4[i]);
    for (unsigned int i = 0; i < histos5.size(); i++) histos.push_back(histos5[i]);
    for (unsigned int i = 0; i < histos6.size(); i++) histos.push_back(histos6[i]);

    TFile* outfile = new TFile("xsechistos.root", "RECREATE");
    for (unsigned int i = 0; i < histos.size(); i++) histos[i]->Write();
    outfile->Close();

}
