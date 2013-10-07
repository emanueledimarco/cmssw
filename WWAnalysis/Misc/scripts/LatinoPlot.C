#include "WWAnalysis/Misc/scripts/BasePlot.h"
#include "WWAnalysis/Misc/scripts/LatinoStyle.h"

class LatinoPlot : public BasePlot {

    public: 
        LatinoPlot() { 
            _yoffset          = 0.05;
            _axisLabelSize    = 0.055;
            _legendStart      = 0.84;
        }

        void Draw(const int &rebin=1) {

            Draw(new TCanvas(),rebin);

        }

        void Draw(TCanvas *c1, const int &rebin=1) {

            gStyle->SetOptStat(0);
            c1->cd();
            c1->Clear();

            TPad *pad1;
            pad1 = new TPad("pad1","pad1",0,0,1,1);
            pad1->Draw();
            pad1->cd();

            RebinHists(rebin);
            THStack *hstack = GetStack(c1->GetLogy());
            TH1F *signal = GetSignalHist();
            TH1F *data   = GetDataHist();

            if(c1->GetLogy()) gPad->SetLogy();
            hstack->Draw("hist");
            if(signal) signal->Draw("hist,same");
            if(data)     data->Draw("ep,same");
            DrawLabels();
            gPad->Update();
            pad1->GetFrame()->DrawClone();



        }

        ClassDef(LatinoPlot,0);



};

// ClassImp(LatinoPlot);
