#ifndef LATINO_BASE_PLOT_H
#define LATINO_BASE_PLOT_H

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

#include "WWAnalysis/Misc/scripts/LatinoStyle.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <math.h>

enum samp                       {  Undefined ,  iHWW ,  iWW    ,  iZJets ,  iTop   , iTT       , iTW      ,  iVV   ,  iWJets ,  iWZ    ,  iZZ   ,  iFakes ,  nSamples  };
const std::string sampNames[] = { "Undefined",  "HWW",  "WW"   ,  "ZJets",  "Top"  , "TT"      , "TW"     ,  "VV"  ,  "WJets",  "WZ"   ,  "ZZ"  ,  "Fakes", "nSamples" };
Color_t sampleColor[]         = { kWhite     , kRed+1, kAzure-9, kGreen+2, kYellow , kYellow   , kYellow-7, kAzure-2, kGray+1, kAzure-2, kAzure-9, kGray+1,  kBlack    };

Color_t autreColors[]         = {kMagenta-9, kBlue-3, kSpring-9, kRed-9, kViolet-9, kCyan-6, kGreen-6, kYellow-7, kOrange-6, kRed+3, kMagenta+2 }; 

std::string nameBySamp(samp a){
    if(int(a) < int(nSamples) && int(a)>0) return sampNames[int(a)];
    else                                   return "undefined";
}

samp sampByName(const std::string &name){
    int index = std::find(sampNames, sampNames+nSamples, name)-sampNames;
    if(index == nSamples) return Undefined; 
    return samp(index);
}


//              1 2 3 4 5 6 7 8 9 101112131415161718
//              1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
//              1    2    3    4    5    6    7    8    9    10             13   14        16   17   18
float xPos[] = {0.22,0.22,0.22,0.42,0.42,0.42,0.22,0.42,0.58,0.76,0.22,0.42,0.58,0.76,0.22,0.42,0.58,0.76};
float yOff[] = {0   ,1   ,2   ,0   ,1   ,2   ,3   ,3   ,3   ,3   ,4   ,4   ,4   ,4   ,5   ,5   ,5   ,5   };

//              1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
//                                                                     12             15
float xPosA[] = {0.22,0.22,0.22,0.42,0.42,0.42,0.22,0.42,0.58,0.22,0.42,0.58,0.22,0.42,0.58};
float yOffA[] = {0   ,1   ,2   ,0   ,1   ,2   ,3   ,3   ,3   ,4   ,4   ,4   ,5   ,5   ,5   };

//              1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
//                                                                11
float xPosB[] = {0.22,0.22,0.22,0.42,0.42,0.42,0.22,0.42,0.22,0.42,0.58};
float yOffB[] = {0   ,1   ,2   ,0   ,1   ,2   ,3   ,3   ,4   ,4   ,4   };




// Retruns the max(bin + err) for all bins
Float_t GetMaximumIncludingErrors(TH1F* h) {
    Float_t maxWithErrors = 0;

    Float_t binHeight;
    for (Int_t i=1; i<=h->GetNbinsX(); i++) {
        binHeight = h->GetBinContent(i) + h->GetBinError(i);
        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
}


class BasePlot {

    public: 
        BasePlot();
        void setMCHist   (const samp &s,        TH1F * h);
        void setMCHist   (const std::string &s, TH1F * h) ;

        void setDataHist (TH1F * h)                        { _data          = h;        } 
        void setHWWHist  (TH1F * h)                        { setMCHist(iHWW  ,h);       } 
        void setWWHist   (TH1F * h)                        { setMCHist(iWW   ,h);       } 
        void setZJetsHist(TH1F * h)                        { setMCHist(iZJets,h);       } 
        void setTopHist  (TH1F * h)                        { setMCHist(iTop  ,h);       } 
        void setTTHist   (TH1F * h)                        { setMCHist(iTT   ,h);       } 
        void setTWHist   (TH1F * h)                        { setMCHist(iTW   ,h);       } 
        void setVVHist   (TH1F * h)                        { setMCHist(iVV   ,h);       } 
        void setWZHist   (TH1F * h)                        { setMCHist(iWZ   ,h);       } 
        void setZZHist   (TH1F * h)                        { setMCHist(iZZ   ,h);       } 
        void setFakesHist(TH1F * h)                        { setMCHist(iFakes,h);       } 
        void setWJetsHist(TH1F * h)                        { setMCHist(iWJets,h);       }

        void setMass(const int &m) {_mass=m;}
        void setMax(const float &m) {_myMax=m;}
        void setLumi(const float &l) { _lumi = l; }
        void setLabel(const TString &s) { _xLabel = s; }
        void setUnits(const TString &s) { _units = s; }
        void setBreakdown(const bool &b = true) { _breakdown = b; }
        void notPrelim() {_prelim = false;};

        void addLabel(const std::string &s);
        TH1F* GetDataHist();
        TH1F *GetSignalHist();
        TH1F *GetSummedMCHist();
        THStack* GetStack(bool isLog);

//         void Draw(const int &rebin=1) { Draw(new TCanvas(),rebin); }
//         void Draw(TCanvas *c1, const int &rebin=1) { c1 = c1; rebin;}

    private: 

        int GetSampCount();
        void RebinHists(const int &rebin);
        void DrawLabels();
        void AxisFonts(TAxis*  axis, TString coordinate, TString title);
        void DrawLegend(Float_t x1, Float_t y1, TH1F*   hist, TString label, TString option);
        void THStackAxisFonts(THStack* h, TString  coordinate, TString  title);
        

        // all of the histograms
        std::vector<TH1F*> _hist;
        std::vector<std::pair<std::string,TH1F*> > _autreHists;
        TH1F* _data;

        // to possibly display
        float    _lumi          ;
        TString  _xLabel        ;
        TString  _units         ;
        TLatex * _extraLabel    ;
        int      _mass          ;

        // plot options
        bool     _breakdown     ;
        bool     _showErr       ;
        bool     _stackSignal   ;

        // storing binning info comes in handy
        int      _nbins         ;
        float    _low           ;
        float    _high          ;

        // defaults which might change based on the type of plot
        Int_t   _labelFont      ;
        Float_t _legendTextSize ;
        Float_t _xoffset        ;
        Float_t _yoffset        ;
        Float_t _labelOffset    ;
        Float_t _axisLabelSize  ;
        Float_t _titleOffset    ;
        Float_t _myMax          ;
        Float_t _legendStart    ;
        bool _prelim    ;

        ClassDef(BasePlot,0);

};

#endif
