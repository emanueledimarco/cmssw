#include "WWAnalysis/Misc/scripts/BasePlot.h"
#include "WWAnalysis/Misc/scripts/LatinoStyle.h"


class MWLPlot : public BasePlot{

    public: 
        MWLPlot() { 
            _legendTextSize   = 0.048;
            _axisLabelSize    = 0.08;
            _titleOffset      = 1.4;

            TStyle *ls = LatinoStyle();
            ls->SetCanvasDefH      (794);
            ls->SetTitleFontSize  (0.045);
            ls->SetTitleY         (0.830);

        }


        void Draw(const int &rebin=1, const bool &showErr=true) {
            Draw(new TCanvas(),rebin, showErr);
        }

        void Draw(TCanvas *c1, const int &rebin=1, const bool &showErr=true) {

            gStyle->SetOptStat(0);
            c1->cd();
            c1->Clear();

            TPad *pad1;
            pad1 = new TPad("pad1","pad1",0,1-0.614609572,1,1);
            pad1->SetTopMargin(0.0983606557);
            pad1->SetBottomMargin(0.025);
            pad1->Draw();
            pad1->cd();

            RebinHists(rebin);
            THStack *hstack = GetStack(c1->GetLogy());
            TH1F *signal = GetSignalHist();
            TH1F *data   = GetDataHist();

            if(c1->GetLogy()) gPad->SetLogy();
            hstack->GetHistogram()->SetLabelSize(0.00,"X");
            hstack->GetHistogram()->SetLabelSize(0.06,"Y");
            hstack->GetHistogram()->SetTitleSize(0.06,"XY");
            hstack->Draw("hist");
            if(signal && ! _stackSignal) signal->Draw("hist,same");
            if(data)     data->Draw("ep,same");
            DrawLabels();
            pad1->GetFrame()->DrawClone();


            TH1F *summed = GetSummedMCHist();

            TH1F *rdat = (TH1F*)data->Clone("rdat");   
            if(gROOT->FindObject("rref")) gROOT->FindObject("rref")->Delete();
            TH1F *rref = new TH1F("rref","rref",
                summed->GetNbinsX(),
                summed->GetBinLowEdge(1),
                summed->GetBinLowEdge(summed->GetNbinsX()+1)
            );
            for (int i = 1, n = rref->GetNbinsX(); i <= n+1; ++i) {
                rref->SetBinContent(i,summed->GetBinContent(i));
                rref->SetBinError(i,summed->GetBinError(i));
            }
            rref->SetTitle("");
            rref->SetLineWidth(0);
            rref->SetFillColor(kGray+1);
            rref->SetFillStyle(1001);
            double absmax = 0;
            for (int i = 0, n = rdat->GetNbinsX(); i <= n+1; ++i) {
                double scale = rref->GetBinContent(i);
                if (scale == 0) {
                    rdat->SetBinContent(i, 0);
                    rref->SetBinContent(i, 10000);
                    rdat->SetBinError(i, 0);
                    rref->SetBinError(i, 0);
                } else {
                    rdat->SetBinContent(i, rdat->GetBinContent(i)/scale);
                    rref->SetBinContent(i, rref->GetBinContent(i)/scale);
                    rdat->SetBinError(i, rdat->GetBinError(i)/scale);
                    rref->SetBinError(i, rref->GetBinError(i)/scale);
                    double mymax = TMath::Max(1.2*fabs(rdat->GetBinContent(i)-1)+1.4*rdat->GetBinError(i), 2.0*rref->GetBinError(i));
                    absmax = TMath::Max(mymax, absmax);
                }
            }

            c1->cd();
            TPad *pad2 = new TPad("pad2","pad2",0,0,1,1-0.614609572);
            pad2->SetTopMargin(0.0261437908);
            pad2->SetBottomMargin(0.392156863);
            pad2->Draw();
            pad2->cd();

            TLine *line = new TLine(rref->GetXaxis()->GetXmin(), 1.0, rref->GetXaxis()->GetXmax(), 1.0);
            line->SetLineColor(kBlack);
            line->SetLineWidth(1);
            line->SetLineStyle(1);

            if(showErr) {
                rref->GetYaxis()->SetRangeUser(TMath::Max(-1.,1.-absmax), TMath::Min(3.,absmax+1.));
                AxisFonts(rref->GetXaxis(), "x", hstack->GetXaxis()->GetTitle());
                rref->GetYaxis()->SetTitle("data/mc");
                rref->GetYaxis()->SetLabelSize(0.10);
                rref->GetYaxis()->SetTitleSize(0.10);
                rref->GetYaxis()->SetTitleOffset(0.85);
                rref->GetXaxis()->SetLabelSize(0.10);
                rref->GetXaxis()->SetTitleSize(0.10);
                rref->GetXaxis()->SetTitleOffset(1.5);
                rref->Draw("E2"); 
                rdat->SetMarkerStyle(20);
                rdat->Draw("E SAME p");
                line->Draw("SAME"); 
                c1->Update();
                pad2->GetFrame()->DrawClone();

            } else {
                rdat->GetYaxis()->SetRangeUser(TMath::Max(-1.,1.-absmax), TMath::Min(3.,absmax+1.));
                rdat->GetYaxis()->SetTitle("data/mc");
                rdat->GetYaxis()->SetLabelSize(0.10);
                rdat->GetYaxis()->SetTitleSize(0.10);
                rdat->GetYaxis()->SetTitleOffset(0.85);
                AxisFonts(rdat->GetXaxis(), "x", hstack->GetXaxis()->GetTitle());
                rdat->GetXaxis()->SetLabelSize(0.10);
                rdat->GetXaxis()->SetTitleSize(0.10);
                rdat->GetXaxis()->SetTitleOffset(1.5);
                rdat->SetTitle("");
                rdat->SetMarkerStyle(20);
                rdat->Draw("E p");
                line->Draw("SAME"); 
                c1->Update();
                pad2->GetFrame()->DrawClone();
            }

        }

//         ClassDef(MWLPlot,0);


};

// ClassImp(MWLPlot);

