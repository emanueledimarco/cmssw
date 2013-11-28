#include<vector>

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "THStack.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TStyle.h"
#endif

#include <iostream>
#include <algorithm>
#include <utility>

enum samp                       {  Undefined ,  iHWW ,  iWW    ,  iZJets ,  iTop   ,  iVV   ,  iWJets ,  iWZ    ,  iZZ   ,  iFakes ,  nSamples  };
const std::string sampNames[] = { "Undefined",  "HWW",  "WW"   ,  "ZJets",  "Top"  ,  "VV"  ,  "WJets",  "WZ"   ,  "ZZ"  ,  "Fakes", "nSamples" };
Color_t sampleColor[]         = { kWhite     , kRed+1, kAzure-9, kGreen+2, kYellow , kAzure-2, kGray+1, kAzure-2, kAzure-9, kGray+1,  kBlack    };

Color_t autreColors[]         = {kMagenta-9, kBlue-3, kSpring-9, kOrange-9, kRed-9, kViolet-9, kCyan-6, kGreen-6, kYellow-7, kOrange-6, kRed+3, kMagenta+2 }; 

std::string nameBySamp(samp a){
    if(int(a) < int(nSamples) && int(a)>0) return sampNames[int(a)];
    else                                   return "undefined";
}

samp sampByName(const std::string &name){
    int index = std::find(sampNames, sampNames+nSamples, name)-sampNames;
    if(index == nSamples) return Undefined; 
    return samp(index);
}


//              1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
//              1    2    3    4    5    6    7    8    9    10             13   14        16   17   18
float xPos[] = {0.22,0.22,0.22,0.39,0.39,0.39,0.22,0.39,0.56,0.73,0.22,0.39,0.56,0.73,0.22,0.39,0.56,0.73}; 
float yOff[] = {0   ,1   ,2   ,0   ,1   ,2   ,3   ,3   ,3   ,3   ,4   ,4   ,4   ,4   ,5   ,5   ,5   ,5   };

//              1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
//                                                                     12             15
float xPosA[] = {0.22,0.22,0.22,0.39,0.39,0.39,0.22,0.39,0.56,0.22,0.39,0.56,0.22,0.39,0.56}; 
float yOffA[] = {0   ,1   ,2   ,0   ,1   ,2   ,3   ,3   ,3   ,4   ,4   ,4   ,5   ,5   ,5   };

//              1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
//                                                                11   
float xPosB[] = {0.22,0.22,0.22,0.39,0.39,0.39,0.22,0.39,0.22,0.39,0.56}; 
float yOffB[] = {0   ,1   ,2   ,0   ,1   ,2   ,3   ,3   ,4   ,4   ,4   };



//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1F* h)
{
    Float_t maxWithErrors = 0;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
}





class LatinoPlot {

    public: 
        LatinoPlot() { 
            _hist.resize(nSamples,0); 
            _data = 0; 
            _breakdown = false; 
            _mass = 0; 
            _nbins=_low=_high=-1;
            _labelFont        = 42;
            _legendTextSize   = 0.03;
            _xoffset          = 0.20;
            _yoffset          = 0.05;
            _labelOffset      = 0.015;
            _axisLabelSize    = 0.050;
            _titleOffset      = 1.6;
            _sf               = 1.0;
        }

        void setDataHist (TH1F * h)                        { _data          = h;        } 
        void setHWWHist  (TH1F * h)                        { setMCHist(iHWW  ,h);       } 
        void setWWHist   (TH1F * h)                        { setMCHist(iWW   ,h);       } 
        void setZJetsHist(TH1F * h)                        { setMCHist(iZJets,h);       } 
        void setTopHist  (TH1F * h)                        { setMCHist(iTop  ,h);       } 
        void setVVHist   (TH1F * h)                        { setMCHist(iVV   ,h);       } 
        void setWZHist   (TH1F * h)                        { setMCHist(iWZ   ,h);       } 
        void setZZHist   (TH1F * h)                        { setMCHist(iZZ   ,h);       } 
        void setFakesHist(TH1F * h)                        { setMCHist(iFakes,h);       } 
        void setWJetsHist(TH1F * h)                        { setMCHist(iWJets,h);       }

        void setMCHist   (const samp &s,        TH1F * h)  { 
            if ( _nbins == -1 ) {
                _nbins  = h->GetNbinsX();
                _low    = h->GetXaxis()->GetBinLowEdge(1);
                _high   = h->GetBinLowEdge(h->GetNbinsX()+1);
            }
            if (h->GetNbinsX() != _nbins || fabs(h->GetXaxis()->GetBinLowEdge(1)-_low) > 0.001 || fabs(h->GetBinLowEdge(h->GetNbinsX()+1)-_high) > 0.001) {
                std::cout << "Warning, we have a bin mismatch problem" << std::endl;
            }
            _hist[s]       = h;        
        } 

        void setMCHist   (const std::string &s, TH1F * h)  { 
            if (sampByName(s) == Undefined) {
                _autreHists.push_back(make_pair(s,h));
            } else {
                setMCHist(sampByName(s),h);
            }
        } 

        void setMass(const int &m) {_mass=m;}

        void Draw(const int &rebin=1,const bool &div=false) {

            _sf = div?1.6:1.0;
            Draw(new TCanvas(),rebin,div);

        }

        void Draw(TCanvas *c1, const int &rebin=1, const bool &div=false) {

            _sf = div?1.3:1.0;
            gStyle->SetOptStat(0);
            c1->cd();
            c1->Clear();

            TPad *pad1;
            if(div) {
                pad1 = new TPad("pad1","pad1",0,1-0.614609572,1,1);
                pad1->SetTopMargin(0.0983606557);
                pad1->SetBottomMargin(0.0163934426);
            } else {
                pad1 = new TPad("pad1","pad1",0,0,1,1);
            }
            pad1->Draw();
            pad1->cd();

            RebinHists(rebin);
            THStack *hstack = GetStack(c1->GetLogy());
            TH1F *signal = GetSignalHist();
            TH1F *data   = GetDataHist();

            cout << "Here" << endl;

            if(c1->GetLogy()) gPad->SetLogy();
            if(div) hstack->GetHistogram()->SetLabelSize(0.00,"X");
            if(div) hstack->GetHistogram()->SetLabelSize(0.06,"Y");
            if(div) hstack->GetHistogram()->SetTitleSize(0.06,"XY");
            hstack->Draw("hist");
            if(signal) signal->Draw("hist,same");
            if(data)     data->Draw("ep,same");
            DrawLabels();
            gPad->Update();

            if(div) {

                TH1F *summed = GetSummedMCHist();

                TH1F *rdat = (TH1F*)data->Clone("rdat");   
                TH1F *rref = (TH1F*)summed->Clone("rref"); 
                double absmax = 0;
                for (int i = 0, n = rdat->GetNbinsX(); i <= n+1; ++i) {
                    double scale = rref->GetBinContent(i);
                    if (scale == 0) {
                        rdat->SetBinContent(i, 0);
                        rref->SetBinContent(i, 0);
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
                line->SetLineColor(kRed);
                line->SetLineWidth(2);

                rref->GetYaxis()->SetRangeUser(TMath::Max(0.,1.-absmax), absmax+1.);
                AxisFonts(rref->GetYaxis(), "y", "ratio");
                AxisFonts(rref->GetXaxis(), "x", hstack->GetXaxis()->GetTitle());
                rref->SetLineWidth(0);
                rref->SetFillColor(kYellow);
                rref->GetYaxis()->SetTitle("data/mc");
                rref->SetLabelSize(0.06*_sf,"XY");
                rref->SetTitleSize(0.06*_sf,"XY");
                rref->SetTitle("");
                pad2->cd();
                rref->Draw("E2");
                line->Draw();
                rdat->SetMarkerStyle(20);
                rdat->Draw("E SAME p");


            }


        }


        TH1F* GetDataHist() { 
        
            if(_data) _data->SetLineColor  (kBlack);
            if(_data) _data->SetMarkerStyle(kFullCircle);
            return _data; 
        }

        TH1F *GetSignalHist() {

            if( _hist[iHWW] ) {
                _hist[iHWW]->SetLineColor(sampleColor[iHWW]);
                _hist[iHWW]->SetLineWidth(3);
            }
            return _hist[iHWW];

        }

        TH1F *GetSummedMCHist() {

            if(gROOT->FindObject("hMC")) gROOT->FindObject("hMC")->Delete();
            TH1F* hMC = new TH1F("hMC","hMC",_nbins,_low,_high);
            hMC->Sumw2();

            for (int i=0; i<nSamples; i++) if( _hist[i] && i != iHWW) hMC->Add(_hist[i]);
            for (size_t i=0; i<_autreHists.size(); i++)               hMC->Add(_autreHists[i].second);

            return hMC;

        }

        THStack* GetStack(bool isLog) {
            THStack* hstack = new THStack();

            float binWidth = 0;
            for (int i=0; i<nSamples; i++) if( _hist[i] && i != iHWW) {

                _hist[i]->SetLineColor(sampleColor[i]);
                _hist[i]->SetFillColor(sampleColor[i]);
                _hist[i]->SetFillStyle(1001);
                binWidth = _hist[i]->GetBinWidth(1);

                hstack->Add(_hist[i]);
            }

            for (size_t i=0; i<_autreHists.size(); i++) {

                _autreHists[i].second->SetLineColor(autreColors[i]);
                _autreHists[i].second->SetFillColor(autreColors[i]);
                _autreHists[i].second->SetFillStyle(1001);

                hstack->Add(_autreHists[i].second);
            }


            hstack->Draw("GOFF");
            hstack->SetTitle("CMS preliminary");

            Float_t theMax = hstack->GetMaximum();
            Float_t theMin = hstack->GetMinimum();

            if (_hist[iHWW]) {
                if (_hist[iHWW]->GetMaximum() > theMax) theMax = _hist[iHWW]->GetMaximum();
                if (_hist[iHWW]->GetMinimum() < theMin) theMin = _hist[iHWW]->GetMinimum();
            }

            if (_data) {
                Float_t dataMax = GetMaximumIncludingErrors(_data);
                if (dataMax > theMax) theMax = dataMax;
            }

            int sampCount = GetSampCount();
            float scaleBy = 1.35 + 0.2*(sampCount>6) + 0.2*(sampCount>10) + 0.2*(sampCount>14);

            if (isLog) {
                theMin = theMin==0?0.1:theMin/10;
                hstack->SetMinimum(theMin);
                hstack->SetMaximum(pow(10,(log(theMax)/log(10)-log(theMin)/log(10)+1)*scaleBy+log(theMin)/log(10)-1));
            } else {
                hstack->SetMaximum(scaleBy * theMax);
            }

            if(_breakdown) {
                THStackAxisFonts(hstack, "y", "entries");
                hstack->GetHistogram()->LabelsOption("v");
            } else {
                THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                if(_units.Sizeof() == 1) {
                    THStackAxisFonts(hstack, "x", _xLabel.Data());
                    THStackAxisFonts(hstack, "y", "entries");
                } else {
                    THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                    THStackAxisFonts(hstack, "y", TString::Format("entries / %.0f %s", binWidth,_units.Data()));
                }
            }
            return hstack;
        }

        void setLumi(const float &l) { _lumi = l; }
        void setLabel(const TString &s) { _xLabel = s; }
        void setUnits(const TString &s) { _units = s; }
        void setBreakdown(const bool &b = true) { _breakdown = b; }
        void addLabel(const std::string &s) {
            _extraLabel = new TLatex(0.9, 0.77, TString(s));
            _extraLabel->SetNDC();
            _extraLabel->SetTextAlign(32);
            _extraLabel->SetTextFont(42);
            _extraLabel->SetTextSize(_legendTextSize*_sf);
        }

    private: 
        int GetSampCount() {
            int sampCount = 0;
            for (int i=0; i<nSamples; i++) if( _hist[i] ) sampCount++;
            for (size_t i=0; i<_autreHists.size(); i++)   sampCount++;
            if  (_data)                                   sampCount++;
            return sampCount;
        }


        void RebinHists(const int &rebin) {

            for (int i=0; i<nSamples; i++) if( _hist[i] )            _hist[i]->Rebin(rebin);
            for (size_t i=0; i<_autreHists.size(); i++) _autreHists[i].second->Rebin(rebin);
            if(_data) _data->Rebin(rebin);

        }


        void DrawLabels() {

            // total mess to get it nice, should be redone
            size_t j=0;
            TString higgsLabel = " HWW";
            if(_mass != 0) higgsLabel.Form(" m_{H}=%d",_mass);

            float *pos,*off;
            int sampCount = GetSampCount();
            if(sampCount == 12 || sampCount == 15) { pos = xPosA; off = yOffA; }
            else if(sampCount == 11 )              { pos = xPosB; off = yOffB; }
            else                                   { pos = xPos;  off = yOff;  }
            if(_data        ) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _data,                  " data",                "lp"); j++; }
            if(_hist[iHWW  ]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iHWW  ]         , higgsLabel,             "l" ); j++; }
            if(_hist[iWW   ]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iWW   ]         , " WW",                  "f" ); j++; }
            if(_hist[iZJets]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iZJets]         , " Z+jets",              "f" ); j++; }
            if(_hist[iTop  ]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iTop  ]         , " top",                 "f" ); j++; }
            if(_hist[iVV   ]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iVV   ]         , " WZ/ZZ",               "f" ); j++; }
            if(_hist[iWJets]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iWJets]         , " W+jets",              "f" ); j++; }
            if(_hist[iWZ   ]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iWZ   ]         , " WZ",                  "f" ); j++; }
            if(_hist[iZZ   ]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iZZ   ]         , " ZZ",                  "f" ); j++; }
            if(_hist[iFakes]) { DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _hist[iFakes]         , " fakes",               "f" ); j++; }
            for(size_t i=0;i<_autreHists.size();++i) {
                                DrawLegend(pos[j], 0.84 - off[j]*_yoffset*_sf, _autreHists[i].second , _autreHists[i].first,   "f" ); j++; 
            }

            TLatex* luminosity = new TLatex(0.9, 0.815, TString::Format("L = %.1f fb^{-1}",_lumi));
            luminosity->SetNDC();
            luminosity->SetTextAlign(32);
            luminosity->SetTextFont(42);
            luminosity->SetTextSize(_legendTextSize*_sf);
            luminosity->Draw("same");
            if(_extraLabel) _extraLabel->Draw("same");
        }

        //------------------------------------------------------------------------------
        // AxisFonts
        //------------------------------------------------------------------------------
        void AxisFonts(TAxis*  axis,
                TString coordinate,
                TString title)
        {
            axis->SetLabelFont  (_labelFont  );
            axis->SetLabelOffset(_labelOffset*_sf);
            axis->SetLabelSize  (_axisLabelSize*_sf);
            axis->SetNdivisions (  505);
            axis->SetTitleFont  (_labelFont);
            axis->SetTitleOffset(  1.5);
            axis->SetTitleSize  (_axisLabelSize*_sf);
        
            if (coordinate == "y") axis->SetTitleOffset(_titleOffset*_sf);
        
            axis->SetTitle(title);
        }
        
        //------------------------------------------------------------------------------
        // DrawLegend
        //------------------------------------------------------------------------------
        void DrawLegend(Float_t x1,
                Float_t y1,
                TH1F*   hist,
                TString label,
                TString option)
        {
            TLegend* legend = new TLegend(x1,
                    y1,
                    x1 + _xoffset*_sf,
                    y1 + _yoffset*_sf);
        
            legend->SetBorderSize(     0);
            legend->SetFillColor (     0);
            legend->SetTextAlign (    12);
            legend->SetTextFont  (_labelFont);
            legend->SetTextSize  (_legendTextSize*_sf);
        
            legend->AddEntry(hist, label.Data(), option.Data());
        
            legend->Draw();
        }


        //------------------------------------------------------------------------------
        // THStackAxisFonts
        //------------------------------------------------------------------------------
        void THStackAxisFonts(THStack* h,
                TString  coordinate,
                TString  title)
        {
            TAxis* axis = NULL;
        
            if (coordinate.Contains("x")) axis = h->GetHistogram()->GetXaxis();
            if (coordinate.Contains("y")) axis = h->GetHistogram()->GetYaxis();
        
            AxisFonts(axis, coordinate, title);
        }
        

        std::vector<TH1F*> _hist;
        std::vector<std::pair<std::string,TH1F*> > _autreHists;
        TH1F* _data;

        //MWL
        float    _lumi;
        TString  _xLabel;
        TString  _units;
        TLatex * _extraLabel;
        bool     _breakdown;
        int      _mass;
        int      _nbins;
        float    _low;
        float    _high;
        Int_t   _labelFont      ;
        Float_t _legendTextSize ;
        Float_t _xoffset        ;
        Float_t _yoffset        ;
        Float_t _labelOffset    ;
        Float_t _axisLabelSize  ;
        Float_t _titleOffset    ;
        Float_t _sf    ;


};


