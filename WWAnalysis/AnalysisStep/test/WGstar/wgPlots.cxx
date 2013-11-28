#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TSystem.h>
#include <math.h>
#include <iostream>

const int nprocesses = 5;
enum Proc { data = -1, wgstar = 0, background = 1 };
const char *processes[nprocesses] = { "DoubleMu_2011", "WGstarMM", "WZ",       "ZZ",       "Top" };
int         proctype [nprocesses] = { data,            wgstar,     background, background, background };
int         proccol  [nprocesses] = { 1,               93,         210,        215,        28 };
TTree *trees[nprocesses];
TFile *fOut = 0;
TString pOut = ""; TCanvas *c1;
double wgstarNorm = 1.6*1.5; //PREP cross section
bool   lumiNorm   = true;

void makeOnePlot(TString expr, TString name, TString bins, TString cut) {
    TH1F *histos[nprocesses];
    double normD = 0, normM = 0;
    for (int i = 0; i < nprocesses; ++i) {
        TString what = TString::Format("%s >> %s_%s(%s)", expr.Data(), name.Data(), processes[i], bins.Data());
        TString mycut = "("+cut+")*scalef";
        if (proctype[i] == wgstar) { mycut += "*(z_m < 12)"; }
        else if (TString(processes[i]) == "WZ") { mycut += "*(z_m > 12)"; }
        else if (TString(processes[i]) == "ZZ") { mycut += "*(z_m > 12)"; }
        trees[i]->Draw(what, mycut);
        histos[i] = (TH1F*) gROOT->FindObject(Form("%s_%s",name.Data(), processes[i]));
        int n = histos[i]->GetNbinsX();
        if (!name.Contains("z_m_low")) {
            histos[i]->SetBinContent(1, histos[i]->GetBinContent(1) + histos[i]->GetBinContent(0));
            histos[i]->SetBinContent(n, histos[i]->GetBinContent(n) + histos[i]->GetBinContent(n+1));
        }
        if (proctype[i] == data) {
            normD += histos[i]->Integral();
            histos[i]->SetMarkerStyle(20);
            histos[i]->SetMarkerSize(1.0);
            histos[i]->SetLineWidth(2);
        } else {
            if (proctype[i] == wgstar) histos[i]->Scale(wgstarNorm);
            histos[i]->SetFillColor(proccol[i]);
            normM += histos[i]->Integral();
        }
    }
    THStack *stack = new THStack(name,expr);
    // should we normalize wgstar?
    for (int i = nprocesses -1; i >= 0; --i) {
        if (!lumiNorm && normM > 0 && normD > 0 && i > 0) histos[i]->Scale(normD/normM);
        gDirectory->WriteTObject(histos[i]);
        if (proctype[i] != data) {
            stack->Add(histos[i]);
        }
    }
    if (!lumiNorm) std::cout << "SF: " << (normM ? normD/normM : -1)<< std::endl;
    gDirectory->WriteTObject(stack);
    c1->Clear();
    histos[0]->Draw("E");
    histos[0]->GetXaxis()->SetTitle(expr);
    if (expr.Contains("rel_iso")) {
        c1->SetLogy(1);
        histos[0]->GetYaxis()->SetRangeUser(0.8, 3*histos[0]->GetMaximum());
    } else {
        c1->SetLogy(0);
        histos[0]->GetYaxis()->SetRangeUser(0, 1.4*TMath::Max(histos[0]->GetMaximum(), histos[1]->GetMaximum()));
    }
    histos[0]->GetYaxis()->SetTitle(lumiNorm ? "Events" : "Events (norm. to data)");
    stack->Draw("SAME");
    histos[0]->Draw("E SAME");
    c1->Print(pOut+"/"+name+".png");
}

void makePlots(TString name, TString selection) {
    TDirectory *dOut = fOut->mkdir(name, selection);
    dOut->cd();
    pOut = "plots/"+name; gSystem->mkdir(pOut); gSystem->Exec("cp ~/php/index.php "+pOut+"/");
    makeOnePlot("z_m",      "z_m_low", "12,0,12",   selection);  
    makeOnePlot("z_m",      "z_m_all", "40,0,120", selection); 
    makeOnePlot("w_mt",     "w_mt",    "25,0,120", selection);
    makeOnePlot("w_pt",     "w_pt",    "20,0,150", selection);
    makeOnePlot("z_pt",     "z_pt",    "20,0,150", selection);
    makeOnePlot("z_m < 15", "norm",    "2,-0.5,1.5", selection);
    makeOnePlot("met",      "met",     "20,0,150", selection);
    makeOnePlot("chmet",      "metch",     "20,0,150", selection);
    makeOnePlot("min(met,chmet)", "metmin",    "20,0,150", selection);
    makeOnePlot("l1_pt",    "l1_pt",     "20,0,150", selection);
    makeOnePlot("l1_eta",   "l1_eta",    "12,-2.4,2.4", selection);
    makeOnePlot("l2_pt",    "l2_pt",     "20,0,150", selection);
    makeOnePlot("l2_eta",   "l2_eta",    "12,-2.4,2.4", selection);
    makeOnePlot("l3_pt",    "l3_pt",     "60,0,60", selection);
    makeOnePlot("l3_eta",   "l3_eta",    "12,-2.4,2.4", selection);
    //makeOnePlot("rel_iso_wM",      "iso_w",    "20,0,0.2", selection);
    makeOnePlot("worse_rel_isoA_z", "iso_z",    "20,0,0.2", selection);
    makeOnePlot("worse_rel_iso_3", "iso_3",    "20,0,0.2", selection);
}

void init(int hadd=0) {
    if (hadd) gSystem->Exec("hadd -f wgStarPlots_Top.root  wgStarPlots_TT.root   wgStarPlots_TW.root   wgStarPlots_TbarW.root");
    for (int i = 0; i < nprocesses; ++i) {
        TFile *fi = TFile::Open(Form("wgStarPlots_%s.root", processes[i]));
        trees[i] = (TTree *) fi->Get("wgTree/probe_tree");
    }
    fOut = TFile::Open("wgStarPlots_Final.root","RECREATE");
}

void wgPlots(int hadd=0) {
    init(hadd);
    c1 = new TCanvas("c1","c1");
    gStyle->SetOptStat(111110);
    TString relax    = "MMM && OS && nleps == 3 ";
    TString baseline = "MMM && OS && nleps == 3 && ncands == 1 && worse_rel_isoA_z < 0.2 && min(met,chmet) > 25 && w_mt > 30 && nbjets == 0";
    TString topcheck = "MMM && OS && nleps == 3 && ncands == 1 && worse_rel_isoA_z < 0.2 && min(met,chmet) > 25 && w_mt > 30 && nbjets >= 1";
    TString sscheck  = "MMM && SS && nleps == 3 && ncands == 1 && worse_rel_isoA_z < 0.2 && min(met,chmet) > 25 && w_mt > 30 && nbjets == 0";
    TString relaxSS  = "MMM && SS && nleps == 3";
    lumiNorm = true;
    //makePlots("MMM",        "MMM");
    makePlots("MMM_3l",     baseline);
    makePlots("MMM_relax",  relax);
    makePlots("MMM_relaxSS",  relaxSS);
    makePlots("MMM_sscheck",  sscheck);
    makePlots("MMM_topcheck",  topcheck);
    //makePlots("MMM_4l",     "MMM && l2_pt > 10 && nleps == 4");
    //makePlots("MMM_3l_0b",  baseline+" && nbjets == 0");
    //makePlots("MMM_3l_0j",  baseline+" && njets == 0");
    //makePlots("MMM_3l_0bj", baseline+" && nbjets == 0 && njets == 0");
    lumiNorm = false;
    makePlots("MMM_3l_m12",  baseline+" && z_m < 12");
    makePlots("MMM_3l_m60",  baseline+" && z_m > 60");
}
