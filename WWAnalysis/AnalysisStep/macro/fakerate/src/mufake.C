#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveText.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TF1.h>
#include <string>
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"

void prep(TH1* hist, TPaveText* t, TCanvas* canvas, int pad, const char* info) {
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.035);
    hist->GetXaxis()->SetTitleSize(0.035);
    hist->GetYaxis()->SetTitleSize(0.035);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    hist->GetYaxis()->SetTitle("Fake Rate (Reco#rightarrow ID+Isolation)");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(0.7);
    hist->SetLineColor(4);
    hist->SetLineWidth(1);
    hist->SetStats(0);
    hist->GetYaxis()->SetRangeUser(0.0,0.4);

    t->AddText("#sqrt{s} = 7 TeV 2011, L = 5 fb^{-1}");
    t->AddText(info);
    t->SetBorderSize(0);
    t->SetFillStyle(0);
    t->SetTextAlign(12);
    t->SetTextFont(42);
    t->SetTextSize(0.03);

    canvas->cd(pad);
    hist->Draw();
    t->Draw();


}

void dofakes(const char* path, const char* outfilename) {
    TFile* file = new TFile(path);
    TTree* tree = (TTree*)file->Get("zllmtree/probe_tree");

    Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
    Int_t binnum = 10;
    Float_t ybins[] = {0.0,1.0,2.0};
    Int_t ybinnum = 2;


    TH2D* hfake = new TH2D("hfake", "", binnum, bins, ybinnum, ybins);
    TH2D* hall  = new TH2D("hall",  "", binnum, bins, ybinnum, ybins);
    TH2D* hpass = new TH2D("hpass", "", binnum, bins, ybinnum, ybins);

    TH1D* hall1  = new TH1D("hden1",  "", binnum, bins);
    TH1D* hall2  = new TH1D("hden2",  "", binnum, bins);
    TH1D* hall3  = new TH1D("hden3",  "", binnum, bins);
    TH1D* hpass1 = new TH1D("hnum1",  "", binnum, bins);
    TH1D* hpass2 = new TH1D("hnum2",  "", binnum, bins);
    TH1D* hpass3 = new TH1D("hnum3",  "", binnum, bins);
    TH1D* hfake1 = new TH1D("hfake1", "", binnum, bins);
    TH1D* hfake2 = new TH1D("hfake2", "", binnum, bins);
    TH1D* hfake3 = new TH1D("hfake3", "", binnum, bins);

    hall  ->Sumw2();
    hall1 ->Sumw2();
    hall2 ->Sumw2();
    hpass ->Sumw2();
    hpass1->Sumw2();
    hpass2->Sumw2();
    hpass3->Sumw2();
    hfake ->Sumw2();
    hfake1->Sumw2();
    hfake2->Sumw2();
    hfake3->Sumw2();

    TBranch *bmet     = tree->GetBranch("met");
    TBranch *bzmass   = tree->GetBranch("zmass");
    TBranch *bpt      = tree->GetBranch("pt");
    TBranch *beta     = tree->GetBranch("eta");
    TBranch *bid      = tree->GetBranch("pfid");
    TBranch *bbdtiso  = tree->GetBranch("bdtIso");
    TBranch *bchiso   = tree->GetBranch("pfIsoChHad04");
    TBranch *bphiso   = tree->GetBranch("pfIsoPhoton04_NoEA");
    TBranch *bnhiso   = tree->GetBranch("pfIsoNHad04_NoEA");
    TBranch *brho     = tree->GetBranch("rhoAA");
    TBranch *brun     = tree->GetBranch("run");
    TBranch *blum     = tree->GetBranch("lumi");
    TBranch *bevt     = tree->GetBranch("event");

    float met     = 0.0;
    float zmass   = 0.0;
    float pt      = 0.0;
    float eta     = 0.0;
    int   id      = 0;
    float bdtiso  = 0.0;
    float chiso   = 0.0;
    float phiso   = 0.0;
    float nhiso   = 0.0;
    float rho     = 0.0;
    int   run     = 0;
    int   lum     = 0;
    int   evt     = 0;

    bmet      ->SetAddress(&met);
    bzmass    ->SetAddress(&zmass);
    bpt       ->SetAddress(&pt);
    beta      ->SetAddress(&eta);
    bid       ->SetAddress(&id);
    bbdtiso   ->SetAddress(&bdtiso);
    bchiso    ->SetAddress(&chiso);
    bphiso    ->SetAddress(&phiso);
    bnhiso    ->SetAddress(&nhiso);
    brho      ->SetAddress(&rho);
    brun      ->SetAddress(&run);
    blum      ->SetAddress(&lum);
    bevt      ->SetAddress(&evt);

    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = MuonEffectiveArea::kMuEAData2011;
    MuonEffectiveArea::MuonEffectiveAreaType   effArea_ = MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;

    for (int i = 0; i < tree->GetEntries(); i++) {
        bmet     ->GetEvent(i);
        bzmass   ->GetEvent(i);
        bpt      ->GetEvent(i);
        beta     ->GetEvent(i);
        bid      ->GetEvent(i);
        bbdtiso  ->GetEvent(i);
        bchiso   ->GetEvent(i);
        bphiso   ->GetEvent(i);
        bnhiso   ->GetEvent(i);
        brho     ->GetEvent(i);
        brun     ->GetEvent(i);
        blum     ->GetEvent(i);
        bevt     ->GetEvent(i);

        float eff_area = MuonEffectiveArea::GetMuonEffectiveArea(effArea_, fabs(eta), effAreaTarget_);

        float iso = (nhiso + phiso - (eff_area)*rho);
        if (iso<0.0) iso = 0.0;
        iso += chiso;
        iso /= pt;

        if (zmass>60 && zmass<120) {
            if (fabs(eta)<1.479) {
                hall1->Fill(pt);
                hall3->Fill(pt);
                hall->Fill(pt,0.5);
                if (id>0 && iso<0.4) {hpass1->Fill(pt); hpass3->Fill(pt);hpass->Fill(pt, 0.5);} 
            }                  
            else {
                hall2->Fill(pt);
                hall3->Fill(pt);
                hall->Fill(pt,1.5);
                if (id>0 && iso<0.4) {hpass2->Fill(pt); hpass3->Fill(pt);hpass->Fill(pt, 1.5);}
            }          
        }        

    }

    hfake1->Divide(hpass1, hall1);
    hfake2->Divide(hpass2, hall2);
    hfake3->Divide(hpass3, hall3);
    hfake ->Divide(hpass,  hall);

    std::string rootfilename = outfilename;
    std::string pdffilename  = outfilename;
    rootfilename += ".root";
    pdffilename  += ".pdf";
        

    TFile* out = new TFile(rootfilename.c_str(), "RECREATE");
    hfake  ->Write();
    hfake1 ->Write();
    hfake2 ->Write();
    hfake3 ->Write();
    hpass  ->Write();
    hpass1 ->Write();
    hpass2 ->Write();
    hpass3 ->Write();
    hall   ->Write();
    hall1  ->Write();
    hall2  ->Write();
    hall3  ->Write();
    out->Close();

    TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");
    plain->SetPadGridX(kTRUE);
    plain->SetPadGridY(kTRUE);
    plain->SetPadRightMargin(0.07);
    plain->SetPadLeftMargin(0.13);
    plain->cd();

    TCanvas* canvas = new TCanvas("fakerate","Fake Rates",500,1500);

    canvas->Divide(1, 3);

    TPaveText* text1  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
    TPaveText* text2  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
    TPaveText* text3  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");

    prep(hfake1, text1,  canvas, 1, "Muons In Barrel"); 
    prep(hfake2, text2,  canvas, 2, "Muons In Endcap"); 
    prep(hfake3, text3,  canvas, 3, "Muons In Barrel + Endcap"); 

    canvas->Print(pdffilename.c_str());

    delete out;
    delete hall1 ;
    delete hall2 ;
    delete hall3 ;
    delete hfake1;
    delete hfake2;
    delete hfake3;
    delete file  ;
    delete text1 ;
    delete text2 ;
    delete text3 ;
    delete canvas;

}

void mufake() {

    dofakes("/home/avartak/CMS/Higgs/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/conf2/hzzTree.root", "mufakes");

}
