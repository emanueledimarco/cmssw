#include "/afs/cern.ch/user/g/gpetrucc/cpp/ii.c"
#include <math.h>

TTree *bkg = TFile::Open("wgStarPlots_DoubleMu_2011.root")->Get("wgTree/probe_tree");
TTree *sig = TFile::Open("wgStarPlots_WGstarMM.root")->Get("wgTree/probe_tree");

TGraph* wgROC(TString expr, TString basecut) {
    TGraph *roc = rocCurve2(expr, sig, basecut+" && OS && z_m < 20", bkg, basecut+" && SS && z_m < 20");
    roc->Draw("APL");
    return roc;
}

TH1F* wgFOM(TString expr, TString basecut, bool highIsGood, double deltaB=0.3) {
    double nS = sig->Draw(expr, "("+basecut+" && OS && z_m < 20)*scalef");
    TH1F *hS = (TH1*) gROOT->FindObject("htemp")->Clone("hS");
    TH1F *hSclone = (TH1*) hS->Clone("hSclone");
    double S0 = hS->Integral();
    iih1(hS,hS->Integral(),highIsGood);
    if (expr.Contains(">>htemp")) bkg->Draw(expr, "("+basecut+" && SS && z_m < 20)" );
    else bkg->Draw(expr+TString::Format(">>htemp(%d,%g,%g)", hS->GetNbinsX(), hS->GetXaxis()->GetXmin(), hS->GetXaxis()->GetXmax()),  "("+basecut+" && SS && z_m < 20)");
    TH1F *hB = (TH1*) gROOT->FindObject("htemp")->Clone("hB");
    TH1F *hBclone = (TH1*) hB->Clone("hBclone");
    double B0 = hB->Integral();
    iih1(hB,hB->Integral(),highIsGood); 
    for (int i = 1; i <= hS->GetNbinsX(); ++i) {
        double S = hS->GetBinContent(i), B = hB->GetBinContent(i);
        std::cout << "For a cut at " << hS->GetBinLowEdge(i) << ": S = " << S << " (eff = " << S/(0.01*S0) << "%); B = " << B << " (eff = " << B/(0.01*B0) << "%)" << std::endl;
        if (B == 0) {
            hS->SetBinContent(i, 0);
        } else {
            hS->SetBinContent(i, S/sqrt(S+B+deltaB*deltaB*B*B));
        }
        hS->SetBinError(i, 0);
    }
    hS->SetLineWidth(3);
    hS->Draw();
    hSclone->SetLineColor(kBlue); hSclone->Scale(hS->GetMaximum()/hSclone->GetMaximum()); hSclone->Draw("SAME");
    hBclone->SetLineColor(kRed);  hBclone->Scale(hS->GetMaximum()/hBclone->GetMaximum()); hBclone->Draw("SAME");
    hS->Draw("SAME");
    return hS;
}


void wgROCS() {
    
}
