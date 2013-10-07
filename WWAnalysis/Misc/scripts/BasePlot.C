#include "WWAnalysis/Misc/scripts/LatinoStyle.h"
#include "WWAnalysis/Misc/scripts/BasePlot.h"


BasePlot::BasePlot() { 
    _hist.resize(nSamples,0); 
    _data = 0; 
    _breakdown = false; 
    _mass = 0; 
    _nbins=-1;
    _low=_high=-1;
    _labelFont        = 42;
    _legendTextSize   = 0.04;
    _xoffset          = 0.20;
    _yoffset          = 0.06;
    _labelOffset      = 0.015;
    _axisLabelSize    = 0.050;
    _titleOffset      = 1.6;
    _myMax = -1.0;
    _stackSignal      = false;
    _prelim           = true;
    _legendStart      = 0.8;

    LatinoStyle();
}

void BasePlot::setMCHist   (const samp &s,        TH1F * h)  { 
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

void BasePlot::setMCHist   (const std::string &s, TH1F * h)  { 
    if (sampByName(s) == Undefined) {
        _autreHists.push_back(make_pair(s,h));
    } else {
        setMCHist(sampByName(s),h);
    }
} 

void BasePlot::addLabel(const std::string &s) {
    _extraLabel = new TLatex(0.89, _legendStart-1*_yoffset, TString(s));
    _extraLabel->SetNDC();
    _extraLabel->SetTextAlign(33);
    _extraLabel->SetTextFont(42);
    _extraLabel->SetTextSize(_legendTextSize*0.9);
}


TH1F* BasePlot::GetDataHist() { 

    if(_data) _data->SetLineColor  (kBlack);
    if(_data) _data->SetMarkerStyle(kFullCircle);
    return _data; 
}

TH1F *BasePlot::GetSignalHist() {

    if( _hist[iHWW] ) {
        _hist[iHWW]->SetLineColor(sampleColor[iHWW]);
        _hist[iHWW]->SetLineWidth(3);
    }
    return _hist[iHWW];

}

TH1F *BasePlot::GetSummedMCHist() {

    if(gROOT->FindObject("hMC")) gROOT->FindObject("hMC")->Delete();
    TH1F* hMC = new TH1F("hMC","hMC",_nbins,_low,_high);
    hMC->Sumw2();

    for (int i=0; i<nSamples; i++) if( _hist[i] && i != iHWW) hMC->Add(_hist[i]);
    for (size_t i=0; i<_autreHists.size(); i++)               hMC->Add(_autreHists[i].second);

    return hMC;

}

THStack* BasePlot::GetStack(bool isLog)  {
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
    if(_prelim) hstack->SetTitle("CMS preliminary");
    else        hstack->SetTitle("CMS, #sqrt{s} = 7 TeV");

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

    // Min --- only need to change if log
    theMin = theMin==0?0.1:theMin/10;
    if(isLog) hstack->SetMinimum(theMin);

    // Max
    if (_myMax != -1) {
        hstack->SetMaximum(_myMax);
    } else if (isLog) {
        hstack->SetMaximum(pow(10,(log(theMax)/log(10)-log(theMin)/log(10)+1)*scaleBy+log(theMin)/log(10)-1));
    } else {
        hstack->SetMaximum(scaleBy * theMax);
    }


    if(_breakdown) {
        THStackAxisFonts(hstack, "y", "entries");
        hstack->GetHistogram()->LabelsOption("v");
    } else {
        THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
        THStackAxisFonts(hstack, "y", TString::Format("entries / %.1f %s", binWidth,_units.Data()));
        if(_units.Sizeof() == 1) {
            THStackAxisFonts(hstack, "x", _xLabel.Data());
//            THStackAxisFonts(hstack, "y", "entries");
        }
    }
    return hstack;
}


int BasePlot::GetSampCount() {
    int sampCount = 0;
    for (int i=0; i<nSamples; i++) if( _hist[i] ) sampCount++;
    for (size_t i=0; i<_autreHists.size(); i++)   sampCount++;
    if  (_data)                                   sampCount++;
    return sampCount;
}

void BasePlot::RebinHists(const int &rebin) {
    for (int i=0; i<nSamples; i++) if( _hist[i] )            _hist[i]->Rebin(rebin);
    for (size_t i=0; i<_autreHists.size(); i++) _autreHists[i].second->Rebin(rebin);
    if(_data) _data->Rebin(rebin);

}

void BasePlot::DrawLabels() {
    // total mess to get it nice, should be redone
    size_t j=0;
    TString higgsLabel = " HWW";
    if(_mass != 0) higgsLabel.Form(" m_{H}=%d",_mass);

    float *pos,*off;
    int sampCount = GetSampCount();
    if(sampCount == 12 || sampCount == 15) { pos = xPosA; off = yOffA; }
    else if(sampCount == 11 )              { pos = xPosB; off = yOffB; }
    else                                   { pos = xPos;  off = yOff;  }
    if(_data        ) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _data,                  " data",                "lp"); j++; }
    if(_hist[iHWW  ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iHWW  ]         , higgsLabel,             "l" ); j++; }
    if(_hist[iWW   ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iWW   ]         , " WW",                  "f" ); j++; }
    if(_hist[iZJets]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iZJets]         , " DY",                  "f" ); j++; }
    if(_hist[iTop  ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iTop  ]         , " top",                 "f" ); j++; }
    if(_hist[iTT   ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iTT   ]         , " t#bar{t}",            "f" ); j++; }
    if(_hist[iTW   ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iTW   ]         , " tW",                  "f" ); j++; }
    if(_hist[iVV   ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iVV   ]         , " WZ/ZZ",               "f" ); j++; }
    if(_hist[iWJets]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iWJets]         , " W+jets",              "f" ); j++; }
    if(_hist[iWZ   ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iWZ   ]         , " WZ",                  "f" ); j++; }
    if(_hist[iZZ   ]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iZZ   ]         , " ZZ",                  "f" ); j++; }
    if(_hist[iFakes]) { DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _hist[iFakes]         , " fakes",               "f" ); j++; }
    for(size_t i=0;i<_autreHists.size();++i) {
                        DrawLegend(pos[j], _legendStart - off[j]*_yoffset, _autreHists[i].second , _autreHists[i].first,   "f" ); j++; 
    }

    TLatex* luminosity = new TLatex(0.896, _legendStart-0.4*_yoffset, TString::Format("L = %.1f fb^{-1}",_lumi).Data());
//     _extraLabel = new TLatex(0.9, , TString(s));
    luminosity->SetNDC();
    luminosity->SetTextAlign(32);
    luminosity->SetTextFont(42);
    luminosity->SetTextSize(_legendTextSize*0.95);
    luminosity->Draw("same");
    if(_extraLabel) _extraLabel->Draw("same");
}

void BasePlot::AxisFonts(TAxis*  axis, TString coordinate, TString title) {
    axis->SetLabelFont  (_labelFont  );
    axis->SetLabelOffset(_labelOffset);
    axis->SetLabelSize  (_axisLabelSize);
    axis->SetNdivisions (  505);
    axis->SetTitleFont  (_labelFont);
    axis->SetTitleOffset(  1.5);
    axis->SetTitleSize  (_axisLabelSize);

    if (coordinate == "y") axis->SetTitleOffset(_titleOffset);

    axis->SetTitle(title);
}

void BasePlot::DrawLegend(Float_t x1, Float_t y1, TH1F*   hist, TString label, TString option) {
    TLegend* legend = new TLegend(x1, y1, x1 + _xoffset, y1 + _yoffset);

    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (_labelFont);
    legend->SetTextSize  (_legendTextSize);

    legend->AddEntry(hist, label.Data(), option.Data());

    legend->Draw();
}

void BasePlot::THStackAxisFonts(THStack* h, TString  coordinate, TString  title) {
    TAxis* axis = NULL;

    if (coordinate.Contains("x")) axis = h->GetHistogram()->GetXaxis();
    if (coordinate.Contains("y")) axis = h->GetHistogram()->GetYaxis();

    AxisFonts(axis, coordinate, title);
}
        
ClassImp(BasePlot);

