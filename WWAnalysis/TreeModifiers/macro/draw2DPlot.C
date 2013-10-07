#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TPad.h>
#include <TFrame.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooLandau.h>
#include <sstream>
#include <vector>
#include <utility>
#include <cmath>
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "SignalInterpolationStrings.h"

using namespace RooFit;

TH2* drawBkgHist(bool is7, int ch) {

    RooRealVar m("m", "", 100, 1600);

    std::stringstream cardss;
    cardss << "/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/TreeModifiers/cards/card_2D_m126_";

    std::string tevstr = is7 ? "7TeV" : "8TeV";
    std::string chstr;
    if (ch == 0) chstr = "4mu";
    if (ch == 1) chstr = "4e";
    if (ch == 2) chstr = "2e2mu";

    cardss << tevstr << "_" << chstr << "_workspace.root";

    TFile* file = TFile::Open("mela2DShapes.root");
    TFile* card = TFile::Open(cardss.str().c_str());

    RooWorkspace* w = (RooWorkspace*)card->Get("w");

    std::stringstream qqzzvarss;
    std::stringstream ggzzvarss;
    std::stringstream zjetvarss;
    qqzzvarss << "bkg_qqzz_" << chstr << "_" << tevstr;
    ggzzvarss << "bkg_ggzz_" << chstr << "_" << tevstr;
    zjetvarss << "bkg_zjets_"<< chstr << "_" << tevstr;

    std::cout << "Reading qq->ZZ shape parameters for " << chstr << " " << tevstr << std::endl;
    RooRealVar qqzz_a0  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a0").c_str()));
    RooRealVar qqzz_a1  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a1").c_str()));
    RooRealVar qqzz_a2  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a2").c_str()));
    RooRealVar qqzz_a3  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a3").c_str()));
    RooRealVar qqzz_a4  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a4").c_str()));
    RooRealVar qqzz_a5  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a5").c_str()));
    RooRealVar qqzz_a6  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a6").c_str()));
    RooRealVar qqzz_a7  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a7").c_str()));
    RooRealVar qqzz_a8  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a8").c_str()));
    RooRealVar qqzz_a9  = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a9").c_str()));
    RooRealVar qqzz_a10 = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a10").c_str()));
    RooRealVar qqzz_a11 = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a11").c_str()));
    RooRealVar qqzz_a12 = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a12").c_str()));
    RooRealVar qqzz_a13 = (*(RooRealVar*)w->var((qqzzvarss.str()+"_a13").c_str()));

    std::cout << "Reading gg->ZZ shape parameters for " << chstr << " " << tevstr << std::endl;
    RooRealVar ggzz_a0  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a0").c_str()));
    RooRealVar ggzz_a1  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a1").c_str()));
    RooRealVar ggzz_a2  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a2").c_str()));
    RooRealVar ggzz_a3  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a3").c_str()));
    RooRealVar ggzz_a4  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a4").c_str()));
    RooRealVar ggzz_a5  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a5").c_str()));
    RooRealVar ggzz_a6  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a6").c_str()));
    RooRealVar ggzz_a7  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a7").c_str()));
    RooRealVar ggzz_a8  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a8").c_str()));
    RooRealVar ggzz_a9  = (*(RooRealVar*)w->var((ggzzvarss.str()+"_a9").c_str()));

    std::cout << "Reading Z+X shape parameters for " << chstr << " " << tevstr << std::endl;
    RooRealVar zjet_me  = (*(RooRealVar*)w->var((zjetvarss.str()+"_mean_zx").c_str()));
    RooRealVar zjet_si  = (*(RooRealVar*)w->var((zjetvarss.str()+"_sigma_zx").c_str()));

    RooqqZZPdf_v2 pdf_qqzz("pdf_qqzz","",m,qqzz_a0,qqzz_a1,qqzz_a2,qqzz_a3,qqzz_a4,qqzz_a5,qqzz_a6,qqzz_a7,qqzz_a8,qqzz_a9, qqzz_a10, qqzz_a11, qqzz_a12, qqzz_a13);
    RooggZZPdf_v2 pdf_ggzz("pdf_ggzz","",m,ggzz_a0,ggzz_a1,ggzz_a2,ggzz_a3,ggzz_a4,ggzz_a5,ggzz_a6,ggzz_a7,ggzz_a8,ggzz_a9);
    RooLandau     pdf_zjet("pdf_zjet","",m,zjet_me,zjet_si);

    std::cout << "Creating histograms from background shapes for " << chstr << " " << tevstr << std::endl;
    TH1* hist_qqzz = pdf_qqzz.createHistogram(("hist_qqzz_"+tevstr+"_"+chstr).c_str(), m, Binning(750, 100, 1600));
    TH1* hist_ggzz = pdf_ggzz.createHistogram(("hist_ggzz_"+tevstr+"_"+chstr).c_str(), m, Binning(750, 100, 1600));
    TH1* hist_zjet = pdf_zjet.createHistogram(("hist_zjet_"+tevstr+"_"+chstr).c_str(), m, Binning(750, 100, 1600));

    std::stringstream qqzz2Dss;
    std::stringstream ggzz2Dss;
    std::stringstream zjet2Dss;
    if (ch == 0) qqzz2Dss << "hist2D_bkg_mm";
    if (ch == 1) qqzz2Dss << "hist2D_bkg_ee";
    if (ch == 2) qqzz2Dss << "hist2D_bkg_em";
    if (ch == 0) ggzz2Dss << "hist2D_b2g_mm";
    if (ch == 1) ggzz2Dss << "hist2D_b2g_ee";
    if (ch == 2) ggzz2Dss << "hist2D_b2g_em";
    if (ch == 0) zjet2Dss << "hist2D_b3g_ee";
    if (ch == 1) zjet2Dss << "hist2D_b3g_em";
    if (ch == 2) zjet2Dss << "hist2D_b3g_mm";

    

    std::cout << "Reading 2D templates for " << chstr << " " << tevstr << std::endl;
    TH2* hist2D_bkg = (TH2*)((TH2*)file->Get(qqzz2Dss.str().c_str()))->Clone(("hist2D_qqzz_"+tevstr+"_"+chstr).c_str());
    TH2* hist2D_b2g = (TH2*)((TH2*)file->Get(ggzz2Dss.str().c_str()))->Clone(("hist2D_ggzz_"+tevstr+"_"+chstr).c_str());
    TH2* hist2D_b3g = (TH2*)(hist2D_bkg->Clone(("hist2D_zjet_"+tevstr+"_"+chstr).c_str()));

    hist2D_bkg->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_b2g->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_b3g->GetXaxis()->SetTitle("m_{4l} [GeV]");

    hist2D_bkg->GetYaxis()->SetTitle("K_{D}");
    hist2D_b2g->GetYaxis()->SetTitle("K_{D}");
    hist2D_b3g->GetYaxis()->SetTitle("K_{D}");

    std::cout << "Taking the product of the template and m(4l) shape for " << chstr << " " << tevstr << std::endl;
    for (int i = 1; i <= 750; i++) {
        for (int j = 1; j <= 30; j++) {
            hist2D_bkg->SetBinContent(i,j,hist2D_bkg->GetBinContent(i, j) * hist_qqzz->GetBinContent(i));
            hist2D_b2g->SetBinContent(i,j,hist2D_b2g->GetBinContent(i, j) * hist_ggzz->GetBinContent(i));
            hist2D_b3g->SetBinContent(i,j,hist2D_b3g->GetBinContent(i, j) * hist_zjet->GetBinContent(i));
        }
    }
    

    float qqzz_fullyield = 0.0;
    float ggzz_fullyield = 0.0;
    float zjet_fullyield = 0.0;
    
    if (ch == 0) {
        qqzz_fullyield = is7 ? 21.0 : 57.0;
        ggzz_fullyield = is7 ? 1.30 : 3.11;
        zjet_fullyield = is7 ? 0.97 : 2.97;
    }
    
    else if (ch == 1) {
        qqzz_fullyield = is7 ? 13.5 : 37.9;
        ggzz_fullyield = is7 ? 0.92 : 2.28;
        zjet_fullyield = is7 ? 1.14 : 6.62;
    }
    
    else {
        qqzz_fullyield = is7 ? 32.2 : 91.6;
        ggzz_fullyield = is7 ? 2.18 : 6.63;
        zjet_fullyield = is7 ? 2.20 : 9.42;
    }

    std::cout << "Scaling 1D background histograms by yields for " << chstr << " " << tevstr << std::endl;
    hist2D_bkg->Scale(qqzz_fullyield);
    hist2D_b2g->Scale(ggzz_fullyield);
    hist2D_b3g->Scale(zjet_fullyield);

    hist2D_bkg->Add(hist2D_b2g);
    hist2D_bkg->Add(hist2D_b3g);

    return hist2D_bkg;
}

TH2* drawSigHist(bool is7, int ch) {

    RooRealVar m("m", "", 100, 180);
    RooRealVar mH("mH", "", 126);
    RooRealVar e1("e1", "", 0);
    RooRealVar e2("e2", "", 0);

    std::string tevstr = is7 ? "7TeV" : "8TeV";
    std::string chstr;
    if (ch == 0) chstr = "4mu";
    if (ch == 1) chstr = "4e";
    if (ch == 2) chstr = "2e2mu";

    TFile* file = TFile::Open("mela2DShapes.root");

    std::stringstream varss;
    varss << "sig_ggh_" << chstr << "_" << tevstr;

    std::cout << "Reading signal shape parameters for " << chstr << " " << tevstr << std::endl;

    RooFormulaVar sig_mean   ("sig_mean_CB" , "", getSignalCBMeanString(126., ch, is7, false).c_str() , RooArgList(mH,e1,e2));
    RooFormulaVar sig_sigma  ("sig_sigma_CB", "", getSignalCBSigmaString(126., ch, is7).c_str()      , RooArgList(mH,e1,e2));
    RooFormulaVar sig_alphaL ("sig_alphaL"  , "", getSignalCBAlphaLString(126., ch, is7).c_str()     , RooArgList(mH));
    RooFormulaVar sig_alphaR ("sig_alphaR"  , "", getSignalCBAlphaRString(126., ch, is7).c_str()     , RooArgList(mH));
    RooFormulaVar sig_nL     ("sig_nL"      , "", getSignalCBNLString(126., ch, is7).c_str()         , RooArgList(mH));
    RooFormulaVar sig_nR     ("sig_nR"      , "", getSignalCBNRString(126., ch, is7).c_str()         , RooArgList(mH));

    RooDoubleCB pdf_sig("pdf_sig","",m,sig_mean,sig_sigma,sig_alphaL,sig_nL,sig_alphaR,sig_nR);

    std::cout << "Creating histograms from signal shapes for " << chstr << " " << tevstr << std::endl;
    TH1* hist_sig = pdf_sig.createHistogram(("hist_sig_"+tevstr+"_"+chstr).c_str(), m, Binning(40, 100, 180));

    std::stringstream sig2Dss;
    if (ch == 0) sig2Dss << "hist2D_sig_mm";
    if (ch == 1) sig2Dss << "hist2D_sig_ee";
    if (ch == 2) sig2Dss << "hist2D_sig_em";

    std::cout << "Reading 2D templates for " << chstr << " " << tevstr << std::endl;
    TH2* hist2D_sig = (TH2*)((TH2*)file->Get(sig2Dss.str().c_str()))->Clone(("hist2D_sig_"+tevstr+"_"+chstr).c_str());
    hist2D_sig->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_sig->GetYaxis()->SetTitle("K_{D}");

    std::cout << "Taking the product of the template and m(4l) shape for " << chstr << " " << tevstr << std::endl;
    for (int i = 1; i <= 750; i++) {
        for (int j = 1; j <= 30; j++) {
            if (i <= 40) hist2D_sig->SetBinContent(i,j,hist2D_sig->GetBinContent(i, j) * hist_sig->GetBinContent(i));
            else         hist2D_sig->SetBinContent(i,j,0);
        }
    }
   

    return hist2D_sig;

}

TGraphErrors* drawData(bool is7, int ch) {

    std::vector<std::pair<unsigned, unsigned> > runeventinfo;
    std::vector<float> massvector;    
    std::vector<float> melavector;    
    std::vector<float> massErrvector;    

    std::string filename; 
    if (is7) filename = "root://pcmssd12.cern.ch//data/hzz4l/step2/HZZ4L_42X_S1_V13_S2_V06/7TeV/DATA/data2011.root"; 
    else     filename = "root://pcmssd12.cern.ch//data/hzz4l/step2/HZZ4L_53X_S1_V13_S2_V03/DATA/8TeV/data2012.root";

    TFile* file = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)file->Get("zz4lTree/probe_tree");

    TFile* ebefile = TFile::Open("/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/TreeModifiers/data/ebeOverallCorrections.HCP2012.v1.root");
    std::string ebemuhistname = is7 ? "mu_reco42x" : "mu_reco53x";
    std::string ebeelhistname = is7 ? "el_reco42x" : "el_reco53x";
    TH2* ebemuhist = (TH2*)ebefile->Get(ebemuhistname.c_str());
    TH2* ebeelhist = (TH2*)ebefile->Get(ebeelhistname.c_str());

    TBranch *bchannel      = tree->GetBranch("channel");
    TBranch *bmass         = tree->GetBranch("mass");
    TBranch *bmela         = tree->GetBranch("ME_SMH_ZZ");
    TBranch *bmassErr      = tree->GetBranch("massErr");
    TBranch *bl1pt         = tree->GetBranch("l1pt");
    TBranch *bl2pt         = tree->GetBranch("l2pt");
    TBranch *bl3pt         = tree->GetBranch("l3pt");
    TBranch *bl4pt         = tree->GetBranch("l4pt");
    TBranch *bl1eta        = tree->GetBranch("l1eta");
    TBranch *bl2eta        = tree->GetBranch("l2eta");
    TBranch *bl3eta        = tree->GetBranch("l3eta");
    TBranch *bl4eta        = tree->GetBranch("l4eta");
    TBranch *bl1massErr    = tree->GetBranch("l1massErr");
    TBranch *bl2massErr    = tree->GetBranch("l2massErr");
    TBranch *bl3massErr    = tree->GetBranch("l3massErr");
    TBranch *bl4massErr    = tree->GetBranch("l4massErr");
    TBranch *bevent        = tree->GetBranch("event");
    TBranch *brun          = tree->GetBranch("run");

    float channel          = 0.0;
    float mass             = 0.0;
    float massErr          = 0.0;
    float mela             = 0.0;
    float l1pt             = 0.0;
    float l2pt             = 0.0;
    float l3pt             = 0.0;
    float l4pt             = 0.0;
    float l1eta            = 0.0;
    float l2eta            = 0.0;
    float l3eta            = 0.0;
    float l4eta            = 0.0;
    float l1massErr        = 0.0;
    float l2massErr        = 0.0;
    float l3massErr        = 0.0;
    float l4massErr        = 0.0;
    unsigned event         = 0;
    unsigned run           = 0;

    bchannel   ->SetAddress(&channel);
    bmass      ->SetAddress(&mass);
    bmassErr   ->SetAddress(&massErr);
    bmela      ->SetAddress(&mela);
    bl1pt      ->SetAddress(&l1pt);
    bl2pt      ->SetAddress(&l2pt);
    bl3pt      ->SetAddress(&l3pt);
    bl4pt      ->SetAddress(&l4pt);
    bl1eta     ->SetAddress(&l1eta);
    bl2eta     ->SetAddress(&l2eta);
    bl3eta     ->SetAddress(&l3eta);
    bl4eta     ->SetAddress(&l4eta);
    bl1massErr ->SetAddress(&l1massErr);
    bl2massErr ->SetAddress(&l2massErr);
    bl3massErr ->SetAddress(&l3massErr);
    bl4massErr ->SetAddress(&l4massErr);
    bevent     ->SetAddress(&event);
    brun       ->SetAddress(&run);

    for (int i = 0; i < tree->GetEntries(); i++) {
        bchannel            ->GetEvent(i);
        bmass               ->GetEvent(i);
        bmassErr            ->GetEvent(i);
        bl1pt               ->GetEvent(i);
        bl2pt               ->GetEvent(i);
        bl3pt               ->GetEvent(i);
        bl4pt               ->GetEvent(i);
        bl1eta              ->GetEvent(i);
        bl2eta              ->GetEvent(i);
        bl3eta              ->GetEvent(i);
        bl4eta              ->GetEvent(i);
        bl1massErr          ->GetEvent(i);
        bl2massErr          ->GetEvent(i);
        bl3massErr          ->GetEvent(i);
        bl4massErr          ->GetEvent(i);
        bevent              ->GetEvent(i);
        brun                ->GetEvent(i);
        bmela               ->GetEvent(i);
    
        bool existsAlready = false;
        for (std::size_t k = 0; k < runeventinfo.size(); k++) {
            if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
        }
   
        if (!existsAlready) {
            runeventinfo.push_back(std::pair<int, int>(run, event));
            if (mass>100. && (channel == ch || (ch == 2 && channel == 3))) {
                massvector.push_back(mass);
                melavector.push_back(mela);

                l1eta = fabs(l1eta);
                l2eta = fabs(l2eta);
                l3eta = fabs(l3eta);
                l4eta = fabs(l4eta);

                float pt1Bin = (channel == 0 || channel == 3) ? std::min(std::max(1,ebemuhist->GetXaxis()->FindBin(l1pt )), ebemuhist->GetNbinsX()) : std::min(std::max(1,ebeelhist->GetXaxis()->FindBin(l1pt )), ebeelhist->GetNbinsX());
                float pt2Bin = (channel == 0 || channel == 3) ? std::min(std::max(1,ebemuhist->GetXaxis()->FindBin(l2pt )), ebemuhist->GetNbinsX()) : std::min(std::max(1,ebeelhist->GetXaxis()->FindBin(l2pt )), ebeelhist->GetNbinsX());
                float pt3Bin = (channel == 0 || channel == 2) ? std::min(std::max(1,ebemuhist->GetXaxis()->FindBin(l3pt )), ebemuhist->GetNbinsX()) : std::min(std::max(1,ebeelhist->GetXaxis()->FindBin(l3pt )), ebeelhist->GetNbinsX());
                float pt4Bin = (channel == 0 || channel == 2) ? std::min(std::max(1,ebemuhist->GetXaxis()->FindBin(l4pt )), ebemuhist->GetNbinsX()) : std::min(std::max(1,ebeelhist->GetXaxis()->FindBin(l4pt )), ebeelhist->GetNbinsX());
                float et1Bin = (channel == 0 || channel == 3) ? std::min(std::max(1,ebemuhist->GetYaxis()->FindBin(l1eta)), ebemuhist->GetNbinsY()) : std::min(std::max(1,ebeelhist->GetYaxis()->FindBin(l1eta)), ebeelhist->GetNbinsY());
                float et2Bin = (channel == 0 || channel == 3) ? std::min(std::max(1,ebemuhist->GetYaxis()->FindBin(l2eta)), ebemuhist->GetNbinsY()) : std::min(std::max(1,ebeelhist->GetYaxis()->FindBin(l2eta)), ebeelhist->GetNbinsY());
                float et3Bin = (channel == 0 || channel == 2) ? std::min(std::max(1,ebemuhist->GetYaxis()->FindBin(l3eta)), ebemuhist->GetNbinsY()) : std::min(std::max(1,ebeelhist->GetYaxis()->FindBin(l3eta)), ebeelhist->GetNbinsY());
                float et4Bin = (channel == 0 || channel == 2) ? std::min(std::max(1,ebemuhist->GetYaxis()->FindBin(l4eta)), ebemuhist->GetNbinsY()) : std::min(std::max(1,ebeelhist->GetYaxis()->FindBin(l4eta)), ebeelhist->GetNbinsY());

                l1massErr *= (channel == 0 || channel == 3) ? ebemuhist->GetBinContent(pt1Bin, et1Bin) : ebeelhist->GetBinContent(pt1Bin, et1Bin);
                l2massErr *= (channel == 0 || channel == 3) ? ebemuhist->GetBinContent(pt2Bin, et2Bin) : ebeelhist->GetBinContent(pt2Bin, et2Bin);
                l3massErr *= (channel == 0 || channel == 2) ? ebemuhist->GetBinContent(pt3Bin, et3Bin) : ebeelhist->GetBinContent(pt3Bin, et3Bin);
                l4massErr *= (channel == 0 || channel == 2) ? ebemuhist->GetBinContent(pt4Bin, et4Bin) : ebeelhist->GetBinContent(pt4Bin, et4Bin);

                massErr = sqrt(l1massErr*l1massErr + l2massErr*l2massErr + l3massErr*l3massErr + l4massErr*l4massErr);

                //std::cout << "Event : " << event << " Run : " << run << " massErr : " << massErr << std::endl;

                massErrvector.push_back(massErr);
            }
        }
 
    }    

    float* x    = new float[massvector.size()];
    float* y    = new float[massvector.size()];
    float* xerr = new float[massvector.size()];
    float* yerr = new float[massvector.size()];
    
    for (std::size_t i = 0; i < massvector.size(); i++) {
        x[i] = massvector[i];
        y[i] = melavector[i];
        xerr[i] = massErrvector[i];
        yerr[i] = 0.0;
    }

    TGraphErrors* gr = new TGraphErrors(massvector.size(), x, y, xerr, yerr);

    return gr;
}


void drawTemplates() {

    TFile* file = TFile::Open("mela2DShapes.root");

    TH2* hist2D_sig_em = (TH2*)file->Get("hist2D_sig_em");
    TH2* hist2D_sig_ee = (TH2*)file->Get("hist2D_sig_ee");
    TH2* hist2D_sig_mm = (TH2*)file->Get("hist2D_sig_mm");
    TH2* hist2D_bkg_em = (TH2*)file->Get("hist2D_bkg_em");
    TH2* hist2D_bkg_ee = (TH2*)file->Get("hist2D_bkg_ee");
    TH2* hist2D_bkg_mm = (TH2*)file->Get("hist2D_bkg_mm");
    TH2* hist2D_b2g_em = (TH2*)file->Get("hist2D_b2g_em");
    TH2* hist2D_b2g_ee = (TH2*)file->Get("hist2D_b2g_ee");
    TH2* hist2D_b2g_mm = (TH2*)file->Get("hist2D_b2g_mm");

    hist2D_sig_em->GetXaxis()->SetRange(1, 40);
    hist2D_sig_ee->GetXaxis()->SetRange(1, 40);
    hist2D_sig_mm->GetXaxis()->SetRange(1, 40);
    hist2D_bkg_em->GetXaxis()->SetRange(1, 40);
    hist2D_bkg_ee->GetXaxis()->SetRange(1, 40);
    hist2D_bkg_mm->GetXaxis()->SetRange(1, 40);
    hist2D_b2g_em->GetXaxis()->SetRange(1, 40);
    hist2D_b2g_ee->GetXaxis()->SetRange(1, 40);
    hist2D_b2g_mm->GetXaxis()->SetRange(1, 40);

    hist2D_sig_em->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_sig_ee->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_sig_mm->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_bkg_em->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_bkg_ee->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_bkg_mm->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_b2g_em->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_b2g_ee->GetXaxis()->SetTitle("m_{4l} [GeV]");
    hist2D_b2g_mm->GetXaxis()->SetTitle("m_{4l} [GeV]");

    hist2D_sig_em->GetYaxis()->SetTitle("MELA");
    hist2D_sig_ee->GetYaxis()->SetTitle("MELA");
    hist2D_sig_mm->GetYaxis()->SetTitle("MELA");
    hist2D_bkg_em->GetYaxis()->SetTitle("MELA");
    hist2D_bkg_ee->GetYaxis()->SetTitle("MELA");
    hist2D_bkg_mm->GetYaxis()->SetTitle("MELA");
    hist2D_b2g_em->GetYaxis()->SetTitle("MELA");
    hist2D_b2g_ee->GetYaxis()->SetTitle("MELA");
    hist2D_b2g_mm->GetYaxis()->SetTitle("MELA");

    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);

    hist2D_sig_em->Draw("COLZ");
    canvas->Print("hist2D_sig_em_lowmass.pdf");
    canvas->Print("hist2D_sig_em_lowmass.png");

    hist2D_sig_ee->Draw("COLZ");
    canvas->Print("hist2D_sig_ee_lowmass.pdf");
    canvas->Print("hist2D_sig_ee_lowmass.png");

    hist2D_sig_mm->Draw("COLZ");
    canvas->Print("hist2D_sig_mm_lowmass.pdf");
    canvas->Print("hist2D_sig_mm_lowmass.png");

    hist2D_bkg_em->Draw("COLZ");
    canvas->Print("hist2D_qqzz_em_lowmass.pdf");
    canvas->Print("hist2D_qqzz_em_lowmass.png");

    hist2D_bkg_ee->Draw("COLZ");
    canvas->Print("hist2D_qqzz_ee_lowmass.pdf");
    canvas->Print("hist2D_qqzz_ee_lowmass.png");

    hist2D_bkg_mm->Draw("COLZ");
    canvas->Print("hist2D_qqzz_mm_lowmass.pdf");
    canvas->Print("hist2D_qqzz_mm_lowmass.png");

    hist2D_b2g_em->Draw("COLZ");
    canvas->Print("hist2D_ggzz_em_lowmass.pdf");
    canvas->Print("hist2D_ggzz_em_lowmass.png");

    hist2D_b2g_ee->Draw("COLZ");
    canvas->Print("hist2D_ggzz_ee_lowmass.pdf");
    canvas->Print("hist2D_ggzz_ee_lowmass.png");

    hist2D_b2g_mm->Draw("COLZ");
    canvas->Print("hist2D_ggzz_mm_lowmass.pdf");
    canvas->Print("hist2D_ggzz_mm_lowmass.png");

    hist2D_sig_em->GetXaxis()->SetRange(21, 350);
    hist2D_sig_ee->GetXaxis()->SetRange(21, 350);
    hist2D_sig_mm->GetXaxis()->SetRange(21, 350);
    hist2D_bkg_em->GetXaxis()->SetRange(21, 350);
    hist2D_bkg_ee->GetXaxis()->SetRange(21, 350);
    hist2D_bkg_mm->GetXaxis()->SetRange(21, 350);
    hist2D_b2g_em->GetXaxis()->SetRange(21, 350);
    hist2D_b2g_ee->GetXaxis()->SetRange(21, 350);
    hist2D_b2g_mm->GetXaxis()->SetRange(21, 350);

    hist2D_sig_em->Draw("COLZ");
    canvas->Print("hist2D_sig_highmass.pdf");
    canvas->Print("hist2D_sig_highmass.png");

    hist2D_bkg_em->Draw("COLZ");
    canvas->Print("hist2D_qqzz_highmass.pdf");
    canvas->Print("hist2D_qqzz_highmass.png");

    hist2D_b2g_em->Draw("COLZ");
    canvas->Print("hist2D_ggzz_highmass.pdf");
    canvas->Print("hist2D_ggzz_highmass.png");

    file->Close();
}

void draw2DPlot() {

    bool plotsig = false;

    Int_t colors[38];
    for (std::size_t i = 0; i < 38; i++) colors[i] = gStyle->GetColorPalette(i)+12;
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);

    TH2* bkghist = drawBkgHist(true, 0);
    bkghist->Add(drawBkgHist(true, 1));
    bkghist->Add(drawBkgHist(true, 2));
    bkghist->Add(drawBkgHist(false, 0));
    bkghist->Add(drawBkgHist(false, 1));
    bkghist->Add(drawBkgHist(false, 2));

    bkghist->GetXaxis()->SetRange(1, 40);

    TH2* sighist = drawSigHist(true, 0);
    sighist->Add(drawSigHist(true, 1));
    sighist->Add(drawSigHist(true, 2));
    sighist->Add(drawSigHist(false, 0));
    sighist->Add(drawSigHist(false, 1));
    sighist->Add(drawSigHist(false, 2));

    sighist->GetXaxis()->SetRange(1, 40);

    TGraphErrors* gr70 = drawData(true , 0);    
    gr70->SetName("gr70");
    TGraphErrors* gr71 = drawData(true , 1);    
    gr71->SetName("gr71");
    TGraphErrors* gr72 = drawData(true , 2);    
    gr72->SetName("gr72");
    TGraphErrors* gr80 = drawData(false, 0);    
    gr80->SetName("gr80");
    TGraphErrors* gr81 = drawData(false, 1);    
    gr81->SetName("gr81");
    TGraphErrors* gr82 = drawData(false, 2); 
    gr82->SetName("gr82");

    TLegend* leg = new TLegend(0.75, 0.75, 0.95, 0.95);
    leg->AddEntry(gr71, "4e"    , "lp");
    leg->AddEntry(gr70, "4#mu"  , "lp");
    leg->AddEntry(gr72, "2e2#mu", "lp");

    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);


    gr72->SetMarkerStyle(21);   
    gr82->SetMarkerStyle(21);   
    gr71->SetMarkerStyle(22);   
    gr81->SetMarkerStyle(22);   

    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    sighist->Draw("COLZ");
    gr70->Draw("P SAME");
    gr80->Draw("P SAME");
    gr71->Draw("P SAME");
    gr81->Draw("P SAME");
    gr72->Draw("P SAME");
    gr82->Draw("P SAME");
    leg->Draw("SAME");
    canvas->Print("plot2DWithData_sig.pdf");
    canvas->Print("plot2DWithData_sig.png");

    bkghist->Draw("COLZ");
    gr70->Draw("P SAME");
    gr80->Draw("P SAME");
    gr71->Draw("P SAME");
    gr81->Draw("P SAME");
    gr72->Draw("P SAME");
    gr82->Draw("P SAME");
    leg->Draw("SAME");
    canvas->Print("plot2DWithData_bkg.pdf");
    canvas->Print("plot2DWithData_bkg.png");

    bkghist->GetXaxis()->SetRange(21, 350);
    bkghist->Draw("COLZ");
    gr70->Draw("P SAME");
    gr80->Draw("P SAME");
    gr71->Draw("P SAME");
    gr81->Draw("P SAME");
    gr72->Draw("P SAME");
    gr82->Draw("P SAME");
    leg->Draw("SAME");
    canvas->Print("plot2DWithData_bkg_highmass.pdf");
    canvas->Print("plot2DWithData_bkg_highmass.png");

    bkghist->GetXaxis()->SetRange(21, 350);
    bkghist->Draw("COLZ");
    gr70->Draw("P SAME");
    gr80->Draw("P SAME");
    gr71->Draw("P SAME");
    gr81->Draw("P SAME");
    gr72->Draw("P SAME");
    gr82->Draw("P SAME");
    leg->Draw("SAME");
    canvas->Print("plot2DWithData_bkg_fullmass.pdf");
    canvas->Print("plot2DWithData_bkg_fullmass.png");

    drawTemplates();
}
