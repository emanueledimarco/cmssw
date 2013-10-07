#ifndef __CINT__
#include "LatinoPlot.C"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

int main(int argc, const char * argv[])
#endif
{

    gInterpreter->ExecuteMacro("scripts/LatinoStyle.C");
#ifdef __CINT__
    gROOT->LoadMacro("scripts/LatinoPlot.C");
#endif 
    TFile* file = new TFile("hists.160.root", "read");

    LatinoPlot myPlot;
    myPlot.setLumi(1300);
    myPlot.setLabel("#Delta#phi");
    myPlot.addLabel("dummy");
    myPlot.setUnits("#circ");

    TH1F* hWW    = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/022.GluGluToWWTo4L");
    TH1F* hZJets = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/004.DYtoMuMu");
    TH1F* hTop   = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/011.tTtoBLNu");
    TH1F* hWZ    = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/001.WZtoAny");
    TH1F* hWJets = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/026.WJetsToLNuMad");

    myPlot.setMCHist(iWW,    (TH1F*)hWW   ->Clone("hWW"));
    myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));
    myPlot.setMCHist(iTop,   (TH1F*)hTop  ->Clone("hTop"));
    myPlot.setMCHist(iWZ,    (TH1F*)hWZ   ->Clone("hWZ")); 
    myPlot.setMCHist(iWJets, (TH1F*)hWJets->Clone("hWJets"));

    // Just as a test
    TH1F *t = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/022.GluGluToWWTo4L");
    myPlot.setDataHist((TH1F*)t->Clone("data"));

    // and just a test
    t = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/101160.ggToH160toWWto2L2Nu");
    t->Scale(0.05);
    myPlot.setMCHist(iHWW, (TH1F*)t->Clone("hww"));

    TCanvas* c1 = new TCanvas("c1", "c1");

    myPlot.Draw(10);  // Can pass a rebin 

    c1->GetFrame()->DrawClone();
    c1->SaveAs("test1.png");


    // Data - MC
    //--------------------------------------------------------------------------
    TCanvas* c2 = new TCanvas("c2", "c2", 575, 10, 550, 600);

    TH1F* dt = (TH1F*)hWW->Clone();
    TH1F* mc = (TH1F*)hWW->Clone();

    dt->SetNameTitle("dt", "dt");
    mc->SetNameTitle("mc", "mc");

    mc->Add(hZJets);
    mc->Add(hTop);
    mc->Add(hWZ);
    mc->Add(hWJets);

    dt->Rebin(10);
    mc->Rebin(10);
    
    TH1F* diff = (TH1F*)dt->Clone();

    diff->SetNameTitle("diff", "diff");

    diff->Reset();

    for (Int_t i=1; i<=diff->GetNbinsX(); i++) {

      Float_t a = dt->GetBinContent(i);
      Float_t b = mc->GetBinContent(i);

      Float_t sa = dt->GetBinError(i);

      Float_t d  = 1e2 * (a - b) / (a + b);
      Float_t sd = 1e2 * 2 * b * sa / (a + b) / (a + b);

      diff->SetBinContent(i, d);
      diff->SetBinError  (i, sd);
    }

    diff->Draw();

    diff->SetMinimum(-50);
    diff->SetMaximum( 50);

    diff->SetMarkerStyle(kFullCircle);

    diff->SetTitle();

    gStyle->SetOptStat(0);

    TBox* tenPercentBox = new TBox(diff->GetXaxis()->GetXmin(), -10,
				   diff->GetXaxis()->GetXmax(), +10);

    tenPercentBox->SetFillColor(kAzure-9);

    tenPercentBox->Draw("same");

    AxisFonts((TAxis*)diff->GetXaxis(), "x", "#Delta#phi [#circ]");
    AxisFonts((TAxis*)diff->GetYaxis(), "y", "100 #upoint (data - MC) / (data + MC)");

    TLine* zeroLine = new TLine(diff->GetXaxis()->GetXmin(), 0,
				diff->GetXaxis()->GetXmax(), 0);

    zeroLine->SetLineStyle(3);
    zeroLine->SetLineWidth(3);

    zeroLine->Draw("same");

    diff->Draw("sameaxis");
    diff->Draw("same");

    c2->GetFrame()->DrawClone();

    c2->SaveAs("test2.png");
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis*  axis,
        TString coordinate,
        TString title)
{
    axis->SetLabelFont  (   42);
    axis->SetLabelOffset(0.015);
    axis->SetLabelSize  (0.050);
    axis->SetNdivisions (  505);
    axis->SetTitleFont  (   42);
    axis->SetTitleOffset(  1.5);
    axis->SetTitleSize  (0.050);

    if (coordinate == "y") axis->SetTitleOffset(1.6);

    axis->SetTitle(title);
}
