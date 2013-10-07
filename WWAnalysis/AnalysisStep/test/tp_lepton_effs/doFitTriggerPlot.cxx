#include <iostream>
#include <cstdio>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TStyle.h>

TF1 *userf1 = 0;
TGraph *theGraph = 0;
FILE *logfile;

TObject *getFromPrefix(TDirectory *dir, TString prefixName, bool allowInTheMiddle=false) {
    if (dir == 0) return 0;
    TIter next(dir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *) next())) {   
        if (key->GetName() == 0) continue;
        TString keyName(key->GetName());
        int index = keyName.Index(prefixName);
        if (allowInTheMiddle ? (index != -1) : (index == 0)) {
            return dir->Get(key->GetName());
        }
    }
    return 0;
}

void sanityCheck(TGraphAsymmErrors *fit) {
    for (int i = 0; i < fit->GetN(); ++i) {
        double x = fit->GetX()[i], y = fit->GetY()[i], eyl = fit->GetErrorYlow(i), eyh = fit->GetErrorYhigh(i);
        if (eyl == 0 && y > 0) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g\n", fit->GetName(), x,y,eyl,eyh);
            if (y < 0.5 && (y-3*eyl < 0)) { eyl = y; } else { eyl = eyh; }
        }
        if (eyl == 0 && y == 1) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g\n", fit->GetName(), x,y,eyl,eyh);
            if (y > 0.5 && (y + 3*eyh > 0)) { eyh = 1-y; } else { eyh = eyl; }
        }
        if ((eyl > 3*eyh) && (y+eyh < 0.9997 || eyh > 0.01)) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g: big eyl\n", fit->GetName(), x,y,eyl,eyh);
            eyl = eyh;
            
        }
        if ((eyh > 3*eyl) && (y-eyl > 0.0003 || eyl > 0.01)) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g: big eyh\n", fit->GetName(), x,y,eyl,eyh);
            eyh = eyl;
        }
        if (eyh > 0.2) eyh = TMath::Min(1-y, (i > 0 ? fit->GetErrorYhigh(i-1) : 0.01));
        if (eyl > 0.2) eyl = TMath::Min(  y, (i > 0 ? fit->GetErrorYlow(i-1)  : 0.01));
        if (i > 0 && i < fit->GetN()-1) {
            double y1 = fit->GetY()[i-1];
            double y2 = fit->GetY()[i+1];
            if (fabs(y - y1) > 0.2 && fabs(y - y2) > 0.2 && ((y - y1)*(y - y2) > 0)) {
                fit->SetPoint(i, x, 0.5*(y1+y2));
            }
        }
        fit->SetPointEYhigh(i, eyh);
        fit->SetPointEYlow(i, eyl);
    }
    fflush(stdout);
}

double cmsprel_xoffs = 0;
double cmsprel_yoffs = 0;
void cmsprelim(double xoffs=-99,double yoffs=-99) {
    bool doSquare = true;
    if (xoffs == -99) xoffs = cmsprel_xoffs;
    if (yoffs == -99) yoffs = cmsprel_yoffs;
    TPaveText *cmsprel = new TPaveText(0.17+xoffs, 0.86-yoffs, 0.84+xoffs, 0.93-yoffs,"NDC");
    cmsprel->SetTextSize(doSquare ? 0.042 : 0.05);
    cmsprel->SetFillColor(0);
    cmsprel->SetFillStyle(0);
    cmsprel->SetLineStyle(2);
    cmsprel->SetLineColor(0);
    cmsprel->SetTextAlign(12);
    cmsprel->SetTextFont(42);
    cmsprel->AddText("CMS Preliminary,  #sqrt{s} = 7 TeV");
    cmsprel->Draw("same");
}
double doLegend_xoffs = 0;
double doLegend_yoffs = 0;
void doLegend(TGraphAsymmErrors *g1, TF1 *g2, TString lab1, TString lab2, TString extraSpam = "") {
    bool doSquare = true;
    double legend_y_offset = 0;
    double legend_y_size   = (extraSpam == "" ? .12 : .18);
    TLegend *leg = new TLegend(doLegend_xoffs+(doSquare ? .62 : .68),.15 + legend_y_offset,doLegend_xoffs+.92,.15 + legend_y_size + legend_y_offset);
    if (extraSpam != "") leg->SetHeader(extraSpam);
    leg->AddEntry(g1, lab1, "LP");
    leg->AddEntry(g2, lab2, "L");
    leg->SetTextSize(doSquare ? 0.04 : 0.05);
    leg->SetTextFont(42);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    if (preliminary != "") cmsprelim();
    leg->Draw();
}



void setFunc(TGraph *graph, TF1 *func, double cut) {
    func->SetLineWidth(2);
    func->SetLineColor(kRed);

    double ymax = 0, ymax30 = 0;
    int N = graph->GetN(); double *x = graph->GetX(), *y = graph->GetY();
    for (int i = 0; i < N; ++i) { 
        if (x[i] < cut) {
            if (y[i] > ymax) ymax = y[i]; 
        } else {
            if (y[i] > ymax30) ymax30 = y[i]; 
        }
    }
    func->SetParameter(0, ymax);
    func->SetParLimits(0, ymax*0.8, ymax*1.2);
    func->SetParameter(3, ymax30);
    func->SetParLimits(3, ymax30*0.8, ymax30*1.2);
    func->FixParameter(4, cut);

    if (cut < 10) {
        // fit with a constant
        func->FixParameter(0, 0);
        func->FixParameter(1, 0);
        func->FixParameter(2, 0);
    } else {
        if (cut > 100) func->FixParameter(3, 0);

        double xmin = 0, xmax = 0;
        for (int i = 0; i < N; ++i) {
            if (y[i] < 0.3*ymax) { xmin = x[i]; } else { break; } 
        }
        for (int i = N-1; i > 0; --i) {
            if (y[i] > 0.3*ymax) { xmax = x[i]; } else { break; } 
        }
        func->SetParameter(1, 0.5*(xmax+xmin));
        func->SetParLimits(1, xmin-1, xmax+1);

        func->SetParameter(2, 0.5*(xmax-xmin)+1);
        func->SetParLimits(2, 0.1, 5);

        std::cout << "Guessed plateau    " << ymax              << std::endl;
        std::cout << "Guessed threshold  " << 0.5*(xmax+xmin)   << std::endl;
        std::cout << "Guessed resolution " << 0.5*(xmax-xmin)+1 << std::endl;
        std::cout << "Guessed pT > " << cut << "  " << ymax30            << std::endl;
    }
}
void doOneFit(TString name, TString mainDir, int etabin) {
    TDirectory *  dir = (TDirectory *) gFile->Get("tpTreeMuMu/"+mainDir+"/fit_eff_plots");
    if (dir == 0) dir = (TDirectory *) gFile->Get("tpTreeElElBDT/"+mainDir+"/fit_eff_plots");
    if (dir == 0) return;
    TCanvas *plot = (TCanvas *) getFromPrefix(dir, Form("pt_PLOT_abseta_bin%d_", etabin));
    if (plot == 0) return;
    TGraphAsymmErrors *eff = (TGraphAsymmErrors*) plot->FindObject("hxy_fit_eff");
    if (eff == 0 || eff->GetN() == 0) return;
    sanityCheck(eff);
    eff->SetLineWidth(2);
    eff->SetMarkerStyle(20);
    eff->SetMarkerSize(1.2);
    int cut = name.Contains("2011A") ? 30 : 40; double xmin = 9.99, xmax = 70.;
    if (name.Contains("Mu8")) {
        if (name.Contains("Mu8leg") || name.Contains("TkMu8") ) cut = 8;
        else cut = 999;
    } else if (name.Contains("SingleEl")) {
        cut = 999; xmin = 19.9; xmax = 110;
    } else if (name.Contains("Ele8")) {
        if (name.Contains("Ele8Leg")) cut = 999;
        else cut = 999;
    }
    TCanvas *c1 = new TCanvas(name,name); //c1->SetGridx(1); c1->SetGridy(1);
    c1->cd(); c1->SetBottomMargin(0.15);
    //TH1F *frame = new TH1F("frame_"+name,name+" turn-on",1,0,30); 
    TH1F *frame = new TH1F("frame_"+name,name+" turn-on",1,xmin,xmax); c1->SetLogx(1);  
    frame->GetYaxis()->SetTitle("Efficiency");
    frame->GetYaxis()->SetDecimals(1);
    frame->GetXaxis()->SetTitleOffset(1.2);
    TLine line; line.SetLineWidth(2.0); line.SetLineStyle(2); line.SetLineColor(202);
    if (name.Contains("SingleMu")) {
        frame->GetXaxis()->SetTitle("muon p_{T} (GeV/c)");
    } else if (name.Contains("Ele")) {
        frame->GetXaxis()->SetTitle("electron p_{T} (GeV/c)");
    }
    c1->cd(); 
    frame->Draw(); 
    frame->GetYaxis()->SetRangeUser(0,1.14);
    frame->GetXaxis()->SetMoreLogLabels(1);
    eff->Draw("P SAME");

    setFunc(eff,userf1, cut);
    eff->Fit(userf1,"B","",xmin,xmax);
    eff->Draw("P SAME");
    c1->SetTitle(name);
    if (name.Contains("SingleMu")) {
        if (name.Contains("175973")) {
            line.DrawLine(24,0,24,1.0);
            line.DrawLine(40,0,40,1.0);
        }
    } else if (name.Contains("Ele17Leg")) {
        line.DrawLine(17,0,17,1.0);
    }
    cmsprelim();
    c1->Print(name+".png");
    c1->Print(name+".pdf");
    fprintf(logfile,"%-20s   %5.3f   %5.2f   %5.3f  %5.3f  %5.3f\n", name.Data(), 
            userf1->GetParameter(0), userf1->GetParameter(1), 
            userf1->GetParameter(2), userf1->GetParameter(4), userf1->GetParameter(3));
    FILE *thistext = fopen((name+".txt").Data(), "w");
    fprintf(thistext, "%.4f/(1 + exp(-(x-%.4f)/%.4f))*(x<%d)+%.4f*(x>%d)\n", 
            userf1->GetParameter(0), userf1->GetParameter(1), 
            userf1->GetParameter(2), 
            cut, userf1->GetParameter(3), cut);
    fclose(thistext);
}

void doFitTriggerPlot() {
    logfile = fopen("fit_results.txt","w");
    fprintf(logfile,"%-20s %5s  %5s   %5s  %5s  %5s\n", "Efficiency", "Plateau", "Turnon", "Resol", "PtCut", "Plat2");
    userf1 = new TF1("userf1","[0]/(1 + exp(-(x-[1])/[2]))*(x<[4])+[3]*(x>[4])", 5, 45);
    userf1->SetParNames("plateu","threshold","resol", "plat2", "cut");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

//#define ALL
#ifdef ALL
    const int   nRunRangesA = 7;
    const char *runRangesA[nRunRangesA] = { "1-163261", "163262-164237", "164238-167043", "167044-167913", "167914-173235", "173236-175910" };
    for (int i = 0; i < nRunRangesA; ++i) {
        TFile *f = TFile::Open(Form("MuonsTrigger_vpvPlusExpo_TnP_Z_DATA_2011A_SingleMu_PtEtaBins_MuFromTk_GiovanniBins_RunRange%s.root", runRangesA[i]));
        if (f == 0) continue;
        doOneFit(Form("2011A_%s_SingleMu_barrel", runRangesA[i]),  "PASSING_all", 0);
        doOneFit(Form("2011A_%s_SingleMu_endcap", runRangesA[i]),  "PASSING_all", 1);
        doOneFit(Form("2011A_%s_SingleMu_forward", runRangesA[i]), "PASSING_all", 2);
        f->Close();
    }

    const int   nRunRangesB = 5;
    const char *runRangesB[nRunRangesB] = {"173236-175910", "175911-175921", "175922-175972", "175973-176544", "176545-999999" };
    for (int i = 0; i < nRunRangesB; ++i) {
        TFile *f = TFile::Open(Form("MuonsTrigger_vpvPlusExpo_TnP_Z_DATA_2011B_SingleMu_PtEtaBins_MuFromTk_GiovanniBins_RunRange%s.root", runRangesB[i]));
        if (f == 0) continue;
        doOneFit(Form("2011B_%s_SingleMu_barrel", runRangesB[i]),  "PASSING_all", 0);
        doOneFit(Form("2011B_%s_SingleMu_endcap", runRangesB[i]),  "PASSING_all", 1);
        doOneFit(Form("2011B_%s_SingleMu_forward", runRangesB[i]), "PASSING_all", 2);
        f->Close();
    }

    const int   nRunRangesDouble = 3;
    const char *runRangesDouble[nRunRangesDouble] = {"165085-175769", "175832-178419", "178420-999999" };
    const char *erasDouble[nRunRangesDouble]      = {"2011A",         "2011B",         "2011B" };
    const char *legsDouble[nRunRangesDouble*2]  = { "Mu13_Mu8leg", "Mu13leg_Mu8", "Mu13_Mu8leg", "Mu13leg_Mu8", "Mu17_Mu8leg", "Mu17leg_Mu8" };
    for (int i = 0; i < nRunRangesDouble*2; ++i) {
        TFile *f = TFile::Open(Form("MuonsTrigger_vpvPlusExpo_TnP_Z_HLT_%s_PtEtaBins_MuFromTk_GiovanniBins_RunRange%s.root", legsDouble[i], runRangesDouble[i/2]));
        if (f == 0) continue;
        doOneFit(Form("%s_%s_%s_barrel",  erasDouble[i/2], runRangesDouble[i/2], legsDouble[i]),  "PASSING_all", 0);
        doOneFit(Form("%s_%s_%s_endcap",  erasDouble[i/2], runRangesDouble[i/2], legsDouble[i]),  "PASSING_all", 1);
        doOneFit(Form("%s_%s_%s_forward", erasDouble[i/2], runRangesDouble[i/2], legsDouble[i]),  "PASSING_all", 2);
        f->Close();
    }

    {
        TFile *f = TFile::Open("MuonsTrigger_vpvPlusExpo_TnP_Z_HLT_Mu17_TkMu8_PtEtaBins_MuFromTk_GiovanniBins_RunRange178420-999999.root");
        if (f == 0) continue;
        doOneFit("2011B_178420-999999_Mu17_TkMu8_Inclusive_barrel",  "PASSING_all", 0);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_Inclusive_endcap",  "PASSING_all", 1);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_Inclusive_forward", "PASSING_all", 2);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_Exclusive_barrel",  "PASSING_this", 0);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_Exclusive_endcap",  "PASSING_this", 1);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_Exclusive_forward", "PASSING_this", 2);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_OR_barrel",  "PASSING_either", 0);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_OR_endcap",  "PASSING_either", 1);
        doOneFit("2011B_178420-999999_Mu17_TkMu8_OR_forward", "PASSING_either", 2);
        f->Close();
    }

#endif
    const int   nRunRangesDoubleE = 3;
    const char *runRangesDoubleE[nRunRangesDoubleE] = { "100001-170052", "170053-175769", "175832-999999" };
    const char *erasDoubleE[nRunRangesDoubleE]      = { "2011A",         "2011A",         "2011B" };
    const char *leg1DoubleE[nRunRangesDoubleE]      = { "Ele17_Ele8_Ele17Leg", "Ele17_Ele8_Tighter_Ele17Leg", "Ele17_Ele8_Tighter_Ele17Leg" };
    const char *leg2DoubleE[nRunRangesDoubleE]      = { "Ele17_Ele8_Ele8Leg",  "Ele17_Ele8_Tighter_Ele8Leg",  "Ele17_Ele8_Tighter_Ele8Leg" };
    const int   nDenominatorDoubleE = 3;
    const char  *denominatorDoubleE[nDenominatorDoubleE] = { "SC17", "SC8", "Ele8" };
    for (int i = 0; i < nRunRangesDoubleE; ++i) {
        for (int j = 0; j < nDenominatorDoubleE; ++j) {
            if (i != 2 || j != 2) continue;
            //TFile *f1 = TFile::Open(Form("ElectronTrigger_vpvPlusExpo_TnP_Z_HLT_%s_from%s_PtEtaBins_GiovanniBins_RunRange%s.root", 
            //                             leg1DoubleE[i], denominatorDoubleE[j], runRangesDoubleE[i]));
            TFile *f1 = TFile::Open(Form("ElectronTriggerV2_vpvPlusExpo_TnP_Z_HLT_%s_from%s_PtEtaBins_GiovanniBins_RunRange%s.root", 
                                         leg1DoubleE[i], denominatorDoubleE[j], runRangesDoubleE[i]));
            if (f1 == 0) continue;
            f1->cd();
            doOneFit(Form("%s_%s_%s_from%s_barrel",  erasDoubleE[i], runRangesDoubleE[i], leg1DoubleE[i], denominatorDoubleE[j]), "PASSING_all", 0);
            doOneFit(Form("%s_%s_%s_from%s_endcap",  erasDoubleE[i], runRangesDoubleE[i], leg1DoubleE[i], denominatorDoubleE[j]), "PASSING_all", 1);
            TFile *f2 = TFile::Open(Form("ElectronTrigger_vpvPlusExpo_TnP_Z_HLT_%s_from%s_PtEtaBins_GiovanniBins_RunRange%s.root", 
                                         leg2DoubleE[i], denominatorDoubleE[j], runRangesDoubleE[i]));
            f1->Close();
            break;
            if (f2 == 0) continue;
            f2->cd();
            doOneFit(Form("%s_%s_%s_from%s_barrel",  erasDoubleE[i], runRangesDoubleE[i], leg2DoubleE[i], denominatorDoubleE[j]), "PASSING_all", 0);
            doOneFit(Form("%s_%s_%s_from%s_endcap",  erasDoubleE[i], runRangesDoubleE[i], leg2DoubleE[i], denominatorDoubleE[j]), "PASSING_all", 1);
            f2->Close();
        }
    }
#ifdef ALL
    const int   nRunRangesSingleE = 6;
    const char *runRangesSingleE[nRunRangesSingleE] = { "100001-164237", "165085-166967", "166968-170901", "170902-175769", "175832-178419", "178420-999999" };
    const char *erasSingleE[nRunRangesSingleE]      = { "2011A",         "2011A",         "2011A",         "2011A",         "2011B",         "2011B" };
    for (int i = 0; i < nRunRangesSingleE; ++i) {
        TFile *f = TFile::Open(Form("ElectronTrigger_vpvPlusExpo_TnP_Z_SingleEl_fromSC17_PtEtaBins_GiovanniBins_RunRange%s.root", runRangesSingleE[i]));
        if (f == 0) continue;
        doOneFit(Form("%s_%s_SingleEl_barrel", erasSingleE[i], runRangesSingleE[i]),  "PASSING_all", 0);
        doOneFit(Form("%s_%s_SingleEl_endcap", erasSingleE[i], runRangesSingleE[i]),  "PASSING_all", 1);
        f->Close();
    }

    fclose(logfile);
#endif
}

#ifdef STANDALONE
#include <TROOT.h>
int main(int argc, char **argv) {
    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    doFitTriggerPlot();
}
#endif
