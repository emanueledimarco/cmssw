#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "YieldMaker.h"

using namespace RooFit;

float melacut           = -1.0;
float z1min             = 40.;
float z2min             = 12.;
float massLow           = 100.;
float massHigh          = 1600.;
int   nBinsMass2D       = 750;
int   nBinsMELA2D       = 30.; 
bool  do7TeV            = true; 
bool  combine2e2mu      = false; 
std::string treeFolder7 = "/home/avartak/CMS/Higgs/Moriond/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";
std::string treeFolder8 = "/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/";

void fillMass(std::string rootfile, float xsecsf, bool isLowMass, int hm, bool i7, TH2F* h2D_sig_mm, TH2F* h2D_sig_ee, TH2F* h2D_sig_em, TH2F* h2D_sig_me) {
    std::string base_folder = do7TeV ? treeFolder7 : treeFolder8;
    ZZYieldMaker   ymaker;

    std::cout << "Filling " << rootfile <<  " with scale factor " << xsecsf << std::endl;

    ymaker.fill(base_folder+rootfile, xsecsf, 0.0, true, 1, hm, i7 );

    if (isLowMass) {
    ymaker.get2DHist(0, z1min, z2min, massLow, 180., melacut, h2D_sig_mm);
    ymaker.get2DHist(1, z1min, z2min, massLow, 180., melacut, h2D_sig_ee);
    ymaker.get2DHist(2, z1min, z2min, massLow, 180., melacut, h2D_sig_em);
    ymaker.get2DHist(3, z1min, z2min, massLow, 180., melacut, h2D_sig_me);
    }


    ymaker.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_sig_mm);
    ymaker.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_sig_mm);
    ymaker.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_sig_mm);
    ymaker.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_sig_mm);

    ymaker.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_sig_ee);
    ymaker.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_sig_ee);
    ymaker.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_sig_ee);
    ymaker.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_sig_ee);

    ymaker.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_sig_em);
    ymaker.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_sig_em);
    ymaker.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_sig_em);
    ymaker.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_sig_em);

    if (!combine2e2mu) {
    ymaker.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_sig_me);
    ymaker.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_sig_me);
    ymaker.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_sig_me);
    ymaker.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_sig_me);
    }

}

void fillBkg(TH2F* h2D_bkg_mm, TH2F* h2D_bkg_ee, TH2F* h2D_bkg_em, TH2F* h2D_bkg_me, bool isGG) {
    ZZYieldMaker   ymaker_ggzz;
    ZZYieldMaker   ymaker_qqzz;

    if (do7TeV) {
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id102.root" ,  getBkgXsec(102)/evt_7TeV(102), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id103.root" ,  getBkgXsec(103)/evt_7TeV(103), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id104.root" ,  getBkgXsec(104)/evt_7TeV(104), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id105.root" ,  getBkgXsec(105)/evt_7TeV(105), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id106.root" ,  getBkgXsec(106)/evt_7TeV(106), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id107.root" ,  getBkgXsec(107)/evt_7TeV(107), 0.0, false);
        ymaker_ggzz.fill(treeFolder7+"hzzTree_id101.root" ,  getBkgXsec(101)/evt_7TeV(101), 0.0, false);
        ymaker_ggzz.fill(treeFolder7+"hzzTree_id100.root" ,  getBkgXsec(100)/evt_7TeV(100), 0.0, false);
    }
    else {
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id102.root" ,  getBkgXsec(102)/evt_8TeV(102), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id103.root" ,  getBkgXsec(103)/evt_8TeV(103), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id104.root" ,  getBkgXsec(104)/evt_8TeV(104), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id105.root" ,  getBkgXsec(105)/evt_8TeV(105), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id106.root" ,  getBkgXsec(106)/evt_8TeV(106), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id107.root" ,  getBkgXsec(107)/evt_8TeV(107), 0.0, false);
        ymaker_ggzz.fill(treeFolder8+"hzzTree_id101.root" ,  getBkgXsec(101)/evt_8TeV(101), 0.0, false);
        ymaker_ggzz.fill(treeFolder8+"hzzTree_id100.root" ,  getBkgXsec(100)/evt_8TeV(100), 0.0, false);
    }

    if (isGG) {
        ymaker_ggzz.get2DHist(0, z1min, z2min, massLow, 180., melacut, h2D_bkg_mm);
        ymaker_ggzz.get2DHist(1, z1min, z2min, massLow, 180., melacut, h2D_bkg_ee);
        ymaker_ggzz.get2DHist(2, z1min, z2min, massLow, 180., melacut, h2D_bkg_em);
        ymaker_ggzz.get2DHist(3, z1min, z2min, massLow, 180., melacut, h2D_bkg_me);
        
        ymaker_ggzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        ymaker_ggzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        ymaker_ggzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        ymaker_ggzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        
        ymaker_ggzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        ymaker_ggzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        ymaker_ggzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        ymaker_ggzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        
        ymaker_ggzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
        ymaker_ggzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
        ymaker_ggzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
        ymaker_ggzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);

        if (!combine2e2mu) {
        ymaker_ggzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        ymaker_ggzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        ymaker_ggzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        ymaker_ggzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        }
    }

    else {
        ymaker_qqzz.get2DHist(0, z1min, z2min, massLow, 180., melacut, h2D_bkg_mm);
        ymaker_qqzz.get2DHist(1, z1min, z2min, massLow, 180., melacut, h2D_bkg_ee);
        ymaker_qqzz.get2DHist(2, z1min, z2min, massLow, 180., melacut, h2D_bkg_em);
        ymaker_qqzz.get2DHist(3, z1min, z2min, massLow, 180., melacut, h2D_bkg_me);
        
        ymaker_qqzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        ymaker_qqzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        ymaker_qqzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        ymaker_qqzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
        
        ymaker_qqzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        ymaker_qqzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        ymaker_qqzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        ymaker_qqzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
        
        ymaker_qqzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
        ymaker_qqzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
        ymaker_qqzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
        ymaker_qqzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);

        if (!combine2e2mu) {
        ymaker_qqzz.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        ymaker_qqzz.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        ymaker_qqzz.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        ymaker_qqzz.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
        }
    }

}

void fillZjets(TH2F* h2D_bkg_mm, TH2F* h2D_bkg_ee, TH2F* h2D_bkg_em, TH2F* h2D_bkg_me, bool doSS) {
    std::string base_folder = do7TeV ? treeFolder7 : treeFolder8;

    FakeRateCalculator FR(base_folder+"data.root", do7TeV, 40, 120, 0.0, 0.0, true);

    ZXYieldMaker   ymaker_zxss;
    ymaker_zxss.fill(base_folder+"data.root", 1.0, FR, doSS);

    ymaker_zxss.get2DHist(0, z1min, z2min, massLow, 180., melacut, h2D_bkg_mm);
    ymaker_zxss.get2DHist(1, z1min, z2min, massLow, 180., melacut, h2D_bkg_ee);
    ymaker_zxss.get2DHist(2, z1min, z2min, massLow, 180., melacut, h2D_bkg_em);
    ymaker_zxss.get2DHist(3, z1min, z2min, massLow, 180., melacut, h2D_bkg_me);

    ymaker_zxss.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
    ymaker_zxss.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
    ymaker_zxss.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);
    ymaker_zxss.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_mm);

    ymaker_zxss.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
    ymaker_zxss.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
    ymaker_zxss.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);
    ymaker_zxss.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_ee);

    ymaker_zxss.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
    ymaker_zxss.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
    ymaker_zxss.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);
    ymaker_zxss.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_em);

    if (!combine2e2mu) {
    ymaker_zxss.get2DHist(0, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
    ymaker_zxss.get2DHist(1, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
    ymaker_zxss.get2DHist(2, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
    ymaker_zxss.get2DHist(3, z1min, z2min, 180., massHigh, melacut, h2D_bkg_me);
    }

}


///// This code is taken from the HZZ4L_Combination package, since we need to rely on the MELA experts for this.

void reweightForInterference(TH2F* hist){
    TF1* reweightFunc =new TF1("reweightFunc","gaus",100,1600);
      
    reweightFunc->SetParameter(0,0.4053);
    reweightFunc->SetParameter(1,110.6);
    reweightFunc->SetParameter(2,22.57);
    
    double oldTempValue=0;
    double newTempValue=0;
    double slope;
    
    for(int i=1; i<=hist->GetNbinsX(); i++){
        for(int j=1; j<=hist->GetNbinsY(); j++){
            slope=reweightFunc->Eval((double)((i-1)*2+101));
            oldTempValue = hist->GetBinContent(i,j);
            newTempValue = oldTempValue*(1+slope*((double)j/30.-.5));
            hist->SetBinContent(i,j,newTempValue);
        }
    }

    delete reweightFunc;
}


void smooth(TH2F* h, int binmin, int binmax, float arraysize) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
    for (int i = binmin; i <= hist->GetNbinsX() && i <= binmax; ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float count = 0.;
            float val = 0;
            for (int m = i-arraysize; m <= i+arraysize; m++) {
                //for (int n = j-arraysize; n <= j+arraysize; n++) {
                //for (int n = j-1; n <= j+1; n++) {
                    int n = j;
                    if (m >= binmin && m <= hist->GetNbinsX() && m <= binmax && n > 0 && n <= hist->GetNbinsY()) {
                        count += 1.0;
                        val += hist->GetBinContent(m, n);
                    }
                //}
            }
            val /= count;
            if (val > 0.0) h->SetBinContent(i,j,val);
            else h->SetBinContent(i,j,0.00001);
        }
    }
}

void removeZeroes(TH2F* hist) {
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            if (hist->GetBinContent(i, j) == 0) hist->SetBinContent(i,j,0.00001);
        }
    }
}


void smoothSwiss(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) {
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.0;;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.0;
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.0;;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.0;
                val /= (9.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}

void smoothVertical(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_vert").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) {
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.0;;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.0;
                val /= (3.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}

void smoothHorizontal(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_hori").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) { 
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.0;;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.0;
                val /= (3.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}


void normalize(TH2F* hist) {
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double histsum = 0;
        for(int j = 1; j <= hist->GetNbinsY(); ++j) histsum += hist->GetBinContent(i,j);
        if (histsum > 0) histsum = 1.0/histsum;
        for(int j = 1; j <= hist->GetNbinsY(); ++j) hist->SetBinContent(i,j, hist->GetBinContent(i,j) * histsum);
    }
}

void normalize(TH2F* hist1, TH2F* hist2) {
    for (int i = 1; i <= hist1->GetNbinsX(); ++i) {
        double hist1sum = 0;
        double hist2sum = 0;
        for(int j = 1; j <= hist1->GetNbinsY(); ++j) hist1sum += hist1->GetBinContent(i,j);
        for(int j = 1; j <= hist2->GetNbinsY(); ++j) hist2sum += hist2->GetBinContent(i,j);
        if (hist1sum > 0) hist1sum = hist2sum/hist1sum;
        for(int j = 1; j <= hist1->GetNbinsY(); ++j) hist1->SetBinContent(i,j, hist1->GetBinContent(i,j) * hist1sum);
    }
}

void create2DHistogram() {

    TH2F* h2D_sig_em;
    TH2F* h2D_sig_me;
    TH2F* h2D_sig_mm;
    TH2F* h2D_sig_ee;

    TH2F* h2D_s2g_em;
    TH2F* h2D_s2g_me;
    TH2F* h2D_s2g_mm;
    TH2F* h2D_s2g_ee;

    TH2F* h2D_bkg_em;
    TH2F* h2D_bkg_me;
    TH2F* h2D_bkg_mm;
    TH2F* h2D_bkg_ee;

    TH2F* h2D_b2g_em;
    TH2F* h2D_b2g_me;
    TH2F* h2D_b2g_mm;
    TH2F* h2D_b2g_ee;

    h2D_sig_em = new TH2F("hist2D_sig_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    if (!combine2e2mu) {
    h2D_sig_me = new TH2F("hist2D_sig_me", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    }
    else {
    h2D_sig_me = h2D_sig_em;    
    }
    h2D_sig_mm = new TH2F("hist2D_sig_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    h2D_sig_ee = new TH2F("hist2D_sig_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);




    h2D_s2g_em = new TH2F("hist2D_s2g_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    if (!combine2e2mu) {
    h2D_s2g_me = new TH2F("hist2D_s2g_me", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    }
    else {
    h2D_s2g_me = h2D_s2g_em;    
    }
    h2D_s2g_mm = new TH2F("hist2D_s2g_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    h2D_s2g_ee = new TH2F("hist2D_s2g_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);




    h2D_bkg_em = new TH2F("hist2D_bkg_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    if (!combine2e2mu) {
    h2D_bkg_me = new TH2F("hist2D_bkg_me", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    }
    else {
    h2D_bkg_me = h2D_bkg_em;    
    }
    h2D_bkg_mm = new TH2F("hist2D_bkg_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    h2D_bkg_ee = new TH2F("hist2D_bkg_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);




    h2D_b2g_em = new TH2F("hist2D_b2g_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    if (!combine2e2mu) {
    h2D_b2g_me = new TH2F("hist2D_b2g_me", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    }
    else {
    h2D_b2g_me = h2D_b2g_em;    
    }
    h2D_b2g_mm = new TH2F("hist2D_b2g_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    h2D_b2g_ee = new TH2F("hist2D_b2g_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);


    h2D_sig_em->Sumw2();
    if (!combine2e2mu) {
    h2D_sig_me->Sumw2();
    }
    h2D_sig_mm->Sumw2();
    h2D_sig_ee->Sumw2();
    
    h2D_s2g_em->Sumw2();
    if (!combine2e2mu) {
    h2D_s2g_me->Sumw2();
    }
    h2D_s2g_mm->Sumw2();
    h2D_s2g_ee->Sumw2();
    
    h2D_bkg_em->Sumw2();
    if (!combine2e2mu) {
    h2D_bkg_me->Sumw2();
    }
    h2D_bkg_mm->Sumw2();
    h2D_bkg_ee->Sumw2();
    
    h2D_b2g_em->Sumw2();
    if (!combine2e2mu) {
    h2D_b2g_me->Sumw2();
    }
    h2D_b2g_mm->Sumw2();
    h2D_b2g_ee->Sumw2();

    do7TeV = true;
    init(do7TeV);

    fillMass("hzzTree_id1115.root" , 1.0/evt_7TeV(1115) , true , 115, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1120.root" , 1.0/evt_7TeV(1120) , true , 120, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1130.root" , 1.0/evt_7TeV(1130) , true , 130, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1140.root" , 1.0/evt_7TeV(1140) , true , 140, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1150.root" , 1.0/evt_7TeV(1150) , true , 150, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1160.root" , 1.0/evt_7TeV(1160) , true , 160, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1170.root" , 1.0/evt_7TeV(1170) , true , 170, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1180.root" , 1.0/evt_7TeV(1180) , true , 180, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1190.root" , 1.0/evt_7TeV(1190) , true , 190, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1200.root" , 1.0/evt_7TeV(1200) , true , 200, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1210.root" , 1.0/evt_7TeV(1210) , true , 210, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1220.root" , 1.0/evt_7TeV(1220) , false, 220, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1230.root" , 1.0/evt_7TeV(1230) , false, 230, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1250.root" , 1.0/evt_7TeV(1250) , false, 250, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1275.root" , 1.0/evt_7TeV(1275) , false, 275, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1300.root" , 1.0/evt_7TeV(1300) , false, 300, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1325.root" , 1.0/evt_7TeV(1325) , false, 325, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1350.root" , 1.0/evt_7TeV(1350) , false, 350, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1375.root" , 1.0/evt_7TeV(1375) , false, 375, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1400.root" , 1.0/evt_7TeV(1400) , false, 400, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1450.root" , 1.0/evt_7TeV(1450) , false, 450, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1500.root" , 1.0/evt_7TeV(1500) , false, 500, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1550.root" , 1.0/evt_7TeV(1550) , false, 550, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1600.root" , 1.0/evt_7TeV(1600) , false, 600, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1650.root" , 1.0/evt_7TeV(1650) , false, 650, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1700.root" , 1.0/evt_7TeV(1700) , false, 700, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1750.root" , 1.0/evt_7TeV(1750) , false, 750, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1800.root" , 1.0/evt_7TeV(1800) , false, 800, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1900.root" , 1.0/evt_7TeV(1900) , false, 900, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1950.root" , 1.0/evt_7TeV(1950) , false, 950, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id11000.root", 1.0/evt_7TeV(11000), false,1000, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);

    std::cout << "Filling background" << std::endl;
    fillBkg (h2D_bkg_mm, h2D_bkg_ee, h2D_bkg_em, h2D_bkg_me, false);
    fillBkg (h2D_b2g_mm, h2D_b2g_ee, h2D_bkg_em, h2D_bkg_me, true);

    do7TeV = false;
    init(do7TeV);

    fillMass("hzzTree_id1115.root" , 1.0/evt_8TeV(1115),  true , 115, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1116.root" , 1.0/evt_8TeV(1116),  true , 116, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1117.root" , 1.0/evt_8TeV(1117),  true , 117, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1118.root" , 1.0/evt_8TeV(1118),  true , 118, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1119.root" , 1.0/evt_8TeV(1119),  true , 119, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1120.root" , 1.0/evt_8TeV(1120),  true , 120, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1121.root" , 1.0/evt_8TeV(1121),  true , 121, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1122.root" , 1.0/evt_8TeV(1122),  true , 122, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1123.root" , 1.0/evt_8TeV(1123),  true , 123, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1124.root" , 1.0/evt_8TeV(1124),  true , 124, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1125.root" , 1.0/evt_8TeV(1125),  true , 125, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1126.root" , 1.0/evt_8TeV(1126),  true , 126, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1127.root" , 1.0/evt_8TeV(1127),  true , 127, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1128.root" , 1.0/evt_8TeV(1128),  true , 128, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1129.root" , 1.0/evt_8TeV(1129),  true , 129, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1130.root" , 1.0/evt_8TeV(1130),  true , 130, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1135.root" , 1.0/evt_8TeV(1135),  true , 135, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1140.root" , 1.0/evt_8TeV(1140),  true , 140, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1145.root" , 1.0/evt_8TeV(1145),  true , 145, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1150.root" , 1.0/evt_8TeV(1150),  true , 150, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1160.root" , 1.0/evt_8TeV(1160),  true , 160, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1170.root" , 1.0/evt_8TeV(1170),  true , 170, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1180.root" , 1.0/evt_8TeV(1180),  true , 180, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1190.root" , 1.0/evt_8TeV(1190),  true , 190, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1200.root" , 1.0/evt_8TeV(1200),  true , 200, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1220.root" , 1.0/evt_8TeV(1220),  false, 220, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1250.root" , 1.0/evt_8TeV(1250),  false, 250, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1275.root" , 1.0/evt_8TeV(1275),  false, 275, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1300.root" , 1.0/evt_8TeV(1300),  false, 300, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1325.root" , 1.0/evt_8TeV(1325),  false, 325, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1350.root" , 1.0/evt_8TeV(1350),  false, 350, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1375.root" , 1.0/evt_8TeV(1375),  false, 375, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1400.root" , 1.0/evt_8TeV(1400),  false, 400, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1450.root" , 1.0/evt_8TeV(1450),  false, 450, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1500.root" , 1.0/evt_8TeV(1500),  false, 500, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1550.root" , 1.0/evt_8TeV(1550),  false, 550, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1600.root" , 1.0/evt_8TeV(1600),  false, 600, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1650.root" , 1.0/evt_8TeV(1650),  false, 650, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1700.root" , 1.0/evt_8TeV(1700),  false, 700, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1750.root" , 1.0/evt_8TeV(1750),  false, 750, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1800.root" , 1.0/evt_8TeV(1800),  false, 800, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1850.root" , 1.0/evt_8TeV(1850),  false, 850, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1900.root" , 1.0/evt_8TeV(1900),  false, 900, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id1950.root" , 1.0/evt_8TeV(1950),  false, 950, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);
    fillMass("hzzTree_id11000.root", 1.0/evt_8TeV(11000), false,1000, do7TeV, h2D_sig_mm, h2D_sig_ee, h2D_sig_em, h2D_sig_me);


    std::cout << "Filling background" << std::endl;
    fillBkg (h2D_bkg_mm, h2D_bkg_ee, h2D_bkg_em, h2D_bkg_me, false);
    fillBkg (h2D_b2g_mm, h2D_b2g_ee, h2D_b2g_em, h2D_bkg_me, true);

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating alternate shapes ...." << std::endl;
    TH2F* h2D_bkg_em_up = ((TH2F*)h2D_bkg_em->Clone("hist2D_bkg_em_up"));
    TH2F* h2D_bkg_me_up = ((TH2F*)h2D_bkg_me->Clone("hist2D_bkg_me_up"));
    TH2F* h2D_bkg_mm_up = ((TH2F*)h2D_bkg_mm->Clone("hist2D_bkg_mm_up"));
    TH2F* h2D_bkg_ee_up = ((TH2F*)h2D_bkg_ee->Clone("hist2D_bkg_ee_up"));

    TH2F* h2D_bkg_em_dn = ((TH2F*)h2D_bkg_em->Clone("hist2D_bkg_em_dn"));
    TH2F* h2D_bkg_me_dn = ((TH2F*)h2D_bkg_me->Clone("hist2D_bkg_me_dn"));
    TH2F* h2D_bkg_mm_dn = ((TH2F*)h2D_bkg_mm->Clone("hist2D_bkg_mm_dn"));
    TH2F* h2D_bkg_ee_dn = ((TH2F*)h2D_bkg_ee->Clone("hist2D_bkg_ee_dn"));

    nBinsMass2D = 6;
    Double_t systbins[] = {100., 120., 140., 160., 180., 300., 1600.};

    TH2F* h2D_zz_em;
    TH2F* h2D_zz_me;
    TH2F* h2D_zz_mm;
    TH2F* h2D_zz_ee;
    TH2F* h2D_zx_em;
    TH2F* h2D_zx_me;
    TH2F* h2D_zx_mm;
    TH2F* h2D_zx_ee;

    h2D_zz_em = new TH2F("hist2D_zz_em", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    if (!combine2e2mu) {
    h2D_zz_me = new TH2F("hist2D_zz_me", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    }
    else {
    h2D_zz_me = h2D_zz_em;
    }
    h2D_zz_mm = new TH2F("hist2D_zz_mm", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    h2D_zz_ee = new TH2F("hist2D_zz_ee", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);

    h2D_zx_em = new TH2F("hist2D_zx_em", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    if (!combine2e2mu) {
    h2D_zx_me = new TH2F("hist2D_zx_me", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    }
    else {
    h2D_zx_me = h2D_zx_em;
    }
    h2D_zx_mm = new TH2F("hist2D_zx_mm", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    h2D_zx_ee = new TH2F("hist2D_zx_ee", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);

    h2D_zx_em->Sumw2();
    if (!combine2e2mu) {
    h2D_zx_me->Sumw2();
    }
    h2D_zx_mm->Sumw2();
    h2D_zx_ee->Sumw2();

    do7TeV = true;
    init(do7TeV);
    fillBkg  (h2D_zz_mm, h2D_zz_ee, h2D_zz_em, h2D_zz_me, false);
    fillZjets(h2D_zx_mm, h2D_zx_ee, h2D_zx_em, h2D_zz_me, true);

    do7TeV = false;
    init(do7TeV);
    fillBkg  (h2D_zz_mm, h2D_zz_ee, h2D_zz_em, h2D_zz_me, false);
    fillZjets(h2D_zx_mm, h2D_zx_ee, h2D_zx_em, h2D_zz_me, true);

    h2D_zz_em->Add(h2D_zz_mm);
    h2D_zz_em->Add(h2D_zz_ee);
    if (!combine2e2mu) {
    h2D_zz_em->Add(h2D_zz_me);
    }
    h2D_zx_em->Add(h2D_zx_mm);
    h2D_zx_em->Add(h2D_zx_ee);
    if (!combine2e2mu) {
    h2D_zx_em->Add(h2D_zx_me);
    }
    normalize(h2D_zz_em, h2D_zx_em);

    std::vector<float> fitp0;
    std::vector<float> fitp1;
    for (int i = 1; i <= nBinsMass2D; i++) {
        stringstream   histname;
        histname   << "hist";
        histname   << i;
        TH1F hist(histname.str().c_str()  , "", nBinsMELA2D, 0., 1.);
        TF1 f("f", "[0]+[1]*x", 0, 1600.);
        for (int j = 1; j <= nBinsMELA2D; j++) {
            hist.SetBinContent  (j, h2D_zx_em->GetBinContent(i,j)/h2D_zz_em->GetBinContent(i,j));
            hist.SetBinError    (j, h2D_zx_em->GetBinError  (i,j)/h2D_zz_em->GetBinContent(i,j));
        }
        hist.Fit(&f);
        fitp0.push_back(f.GetParameter(0));
        fitp1.push_back(f.GetParameter(1));
    }

    for (int i = 1; i <= h2D_bkg_em->GetNbinsX(); i++) {
        for (int j = 1; j <= h2D_bkg_em->GetNbinsY(); j++) {
            int k = nBinsMass2D;
            for (int m = 0; m < nBinsMass2D; m++) {
                if ((massLow + i*2) > systbins[m] && (massLow + i*2) <= systbins[m+1]) {
                    k = m;   
                }       
            }

            if (k < nBinsMass2D) {

                h2D_bkg_me_up->SetBinContent(i, j, h2D_bkg_me->GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));
                h2D_bkg_em_up->SetBinContent(i, j, h2D_bkg_em->GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));
                h2D_bkg_ee_up->SetBinContent(i, j, h2D_bkg_ee->GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));
                h2D_bkg_mm_up->SetBinContent(i, j, h2D_bkg_mm->GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));

                h2D_bkg_me_dn->SetBinContent(i, j, h2D_bkg_me->GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));
                h2D_bkg_em_dn->SetBinContent(i, j, h2D_bkg_em->GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));
                h2D_bkg_ee_dn->SetBinContent(i, j, h2D_bkg_ee->GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));
                h2D_bkg_mm_dn->SetBinContent(i, j, h2D_bkg_mm->GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em->GetNbinsY())));
            }
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Smoothing ...." << std::endl;

    removeZeroes(h2D_sig_em);
    if (!combine2e2mu) {
    removeZeroes(h2D_sig_me);
    }
    removeZeroes(h2D_sig_mm);
    removeZeroes(h2D_sig_ee);

    removeZeroes(h2D_bkg_em);
    if (!combine2e2mu) {
    removeZeroes(h2D_bkg_me);
    }
    removeZeroes(h2D_bkg_mm);
    removeZeroes(h2D_bkg_ee);
    
    removeZeroes(h2D_b2g_em);
    if (!combine2e2mu) {
    removeZeroes(h2D_b2g_me);
    }
    removeZeroes(h2D_b2g_mm);
    removeZeroes(h2D_b2g_ee);
    
    removeZeroes(h2D_bkg_em_up);
    if (!combine2e2mu) {
    removeZeroes(h2D_bkg_me_up);
    }
    removeZeroes(h2D_bkg_mm_up);
    removeZeroes(h2D_bkg_ee_up);

    removeZeroes(h2D_bkg_em_dn);
    if (!combine2e2mu) {
    removeZeroes(h2D_bkg_me_dn);
    }
    removeZeroes(h2D_bkg_mm_dn);
    removeZeroes(h2D_bkg_ee_dn);



    normalize(h2D_sig_em);
    if (!combine2e2mu) {
    normalize(h2D_sig_me);
    }
    normalize(h2D_sig_mm);
    normalize(h2D_sig_ee);

    normalize(h2D_bkg_em);
    if (!combine2e2mu) {
    normalize(h2D_bkg_me);
    }
    normalize(h2D_bkg_mm);
    normalize(h2D_bkg_ee);

    normalize(h2D_b2g_em);
    if (!combine2e2mu) {
    normalize(h2D_b2g_me);
    }
    normalize(h2D_b2g_mm);
    normalize(h2D_b2g_ee);

    normalize(h2D_bkg_em_up);
    if (!combine2e2mu) {
    normalize(h2D_bkg_me_up);
    }
    normalize(h2D_bkg_mm_up);
    normalize(h2D_bkg_ee_up);

    normalize(h2D_bkg_em_dn);
    if (!combine2e2mu) {
    normalize(h2D_bkg_me_dn);
    }
    normalize(h2D_bkg_mm_dn);
    normalize(h2D_bkg_ee_dn);


    int binmin1 = 100;
    int binmin2 = 125;
    int binmin3 = 200;
    int binmin4 = 250;
    int binmin5 = 350;
    int binmin6 = 450;
    int binmin7 = 550;

    int binmax1 = 124;
    int binmax2 = 199;
    int binmax3 = 249;
    int binmax4 = 349;
    int binmax5 = 449;
    int binmax6 = 549;
    int binmax7 = 750;

    smooth(h2D_sig_em,    binmin1, binmax1, 1);
    if (!combine2e2mu) {
    smooth(h2D_sig_me,    binmin1, binmax1, 1);
    }
    smooth(h2D_sig_mm,    binmin1, binmax1, 1);
    smooth(h2D_sig_ee,    binmin1, binmax1, 1);
    
    smooth(h2D_bkg_em,    binmin1, binmax1, 1);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me,    binmin1, binmax1, 1);
    }
    smooth(h2D_bkg_mm,    binmin1, binmax1, 1);
    smooth(h2D_bkg_ee,    binmin1, binmax1, 1);

    smooth(h2D_b2g_em,    binmin1, binmax1, 1);
    if (!combine2e2mu) {
    smooth(h2D_b2g_me,    binmin1, binmax1, 1);
    }
    smooth(h2D_b2g_mm,    binmin1, binmax1, 1);
    smooth(h2D_b2g_ee,    binmin1, binmax1, 1);

    smooth(h2D_bkg_em_up, binmin1, binmax1, 1);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_up, binmin1, binmax1, 1);
    }
    smooth(h2D_bkg_mm_up, binmin1, binmax1, 1);
    smooth(h2D_bkg_ee_up, binmin1, binmax1, 1);

    smooth(h2D_bkg_em_dn, binmin1, binmax1, 1);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_dn, binmin1, binmax1, 1);
    }
    smooth(h2D_bkg_mm_dn, binmin1, binmax1, 1);
    smooth(h2D_bkg_ee_dn, binmin1, binmax1, 1);

    smooth(h2D_sig_em,    binmin2, binmax2, 3);
    if (!combine2e2mu) {
    smooth(h2D_sig_me,    binmin2, binmax2, 3);
    }
    smooth(h2D_sig_mm,    binmin2, binmax2, 3);
    smooth(h2D_sig_ee,    binmin2, binmax2, 3);
    
    smooth(h2D_bkg_em,    binmin2, binmax2, 3);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me,    binmin2, binmax2, 3);
    }
    smooth(h2D_bkg_mm,    binmin2, binmax2, 3);
    smooth(h2D_bkg_ee,    binmin2, binmax2, 3);

    smooth(h2D_b2g_em,    binmin2, binmax2, 3);
    if (!combine2e2mu) {
    smooth(h2D_b2g_me,    binmin2, binmax2, 3);
    }
    smooth(h2D_b2g_mm,    binmin2, binmax2, 3);
    smooth(h2D_b2g_ee,    binmin2, binmax2, 3);

    smooth(h2D_bkg_em_up, binmin2, binmax2, 3);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_up, binmin2, binmax2, 3);
    }
    smooth(h2D_bkg_mm_up, binmin2, binmax2, 3);
    smooth(h2D_bkg_ee_up, binmin2, binmax2, 3);

    smooth(h2D_bkg_em_dn, binmin2, binmax2, 3);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_dn, binmin2, binmax2, 3);
    }
    smooth(h2D_bkg_mm_dn, binmin2, binmax2, 3);
    smooth(h2D_bkg_ee_dn, binmin2, binmax2, 3);

    smooth(h2D_sig_em,    binmin3, binmax3, 5);
    if (!combine2e2mu) {
    smooth(h2D_sig_me,    binmin3, binmax3, 5);
    }
    smooth(h2D_sig_mm,    binmin3, binmax3, 5);
    smooth(h2D_sig_ee,    binmin3, binmax3, 5);
    
    smooth(h2D_bkg_em,    binmin3, binmax3, 5);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me,    binmin3, binmax3, 5);
    }
    smooth(h2D_bkg_mm,    binmin3, binmax3, 5);
    smooth(h2D_bkg_ee,    binmin3, binmax3, 5);

    smooth(h2D_b2g_em,    binmin3, binmax3, 5);
    if (!combine2e2mu) {
    smooth(h2D_b2g_me,    binmin3, binmax3, 5);
    }
    smooth(h2D_b2g_mm,    binmin3, binmax3, 5);
    smooth(h2D_b2g_ee,    binmin3, binmax3, 5);

    smooth(h2D_bkg_em_up, binmin3, binmax3, 5);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_up, binmin3, binmax3, 5);
    }
    smooth(h2D_bkg_mm_up, binmin3, binmax3, 5);
    smooth(h2D_bkg_ee_up, binmin3, binmax3, 5);

    smooth(h2D_bkg_em_dn, binmin3, binmax3, 5);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_dn, binmin3, binmax3, 5);
    }
    smooth(h2D_bkg_mm_dn, binmin3, binmax3, 5);
    smooth(h2D_bkg_ee_dn, binmin3, binmax3, 5);

    smooth(h2D_sig_em,    binmin4, binmax4, 7);
    if (!combine2e2mu) {
    smooth(h2D_sig_me,    binmin4, binmax4, 7);
    }
    smooth(h2D_sig_mm,    binmin4, binmax4, 7);
    smooth(h2D_sig_ee,    binmin4, binmax4, 7);
    
    smooth(h2D_bkg_em,    binmin4, binmax4, 7);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me,    binmin4, binmax4, 7);
    }
    smooth(h2D_bkg_mm,    binmin4, binmax4, 7);
    smooth(h2D_bkg_ee,    binmin4, binmax4, 7);

    smooth(h2D_b2g_em,    binmin4, binmax4, 7);
    if (!combine2e2mu) {
    smooth(h2D_b2g_me,    binmin4, binmax4, 7);
    }
    smooth(h2D_b2g_mm,    binmin4, binmax4, 7);
    smooth(h2D_b2g_ee,    binmin4, binmax4, 7);

    smooth(h2D_bkg_em_up, binmin4, binmax4, 7);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_up, binmin4, binmax4, 7);
    }
    smooth(h2D_bkg_mm_up, binmin4, binmax4, 7);
    smooth(h2D_bkg_ee_up, binmin4, binmax4, 7);

    smooth(h2D_bkg_em_dn, binmin4, binmax4, 7);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_dn, binmin4, binmax4, 7);
    }
    smooth(h2D_bkg_mm_dn, binmin4, binmax4, 7);
    smooth(h2D_bkg_ee_dn, binmin4, binmax4, 7);

    smooth(h2D_sig_em,    binmin5, binmax5, 11);
    if (!combine2e2mu) {
    smooth(h2D_sig_me,    binmin5, binmax5, 11);
    }
    smooth(h2D_sig_mm,    binmin5, binmax5, 11);
    smooth(h2D_sig_ee,    binmin5, binmax5, 11);
    
    smooth(h2D_bkg_em,    binmin5, binmax5, 11);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me,    binmin5, binmax5, 11);
    }
    smooth(h2D_bkg_mm,    binmin5, binmax5, 11);
    smooth(h2D_bkg_ee,    binmin5, binmax5, 11);

    smooth(h2D_b2g_em,    binmin5, binmax5, 11);
    if (!combine2e2mu) {
    smooth(h2D_b2g_me,    binmin5, binmax5, 11);
    }
    smooth(h2D_b2g_mm,    binmin5, binmax5, 11);
    smooth(h2D_b2g_ee,    binmin5, binmax5, 11);

    smooth(h2D_bkg_em_up, binmin5, binmax5, 11);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_up, binmin5, binmax5, 11);
    }
    smooth(h2D_bkg_mm_up, binmin5, binmax5, 11);
    smooth(h2D_bkg_ee_up, binmin5, binmax5, 11);

    smooth(h2D_bkg_em_dn, binmin5, binmax5, 11);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_dn, binmin5, binmax5, 11);
    }
    smooth(h2D_bkg_mm_dn, binmin5, binmax5, 11);
    smooth(h2D_bkg_ee_dn, binmin5, binmax5, 11);

    smooth(h2D_sig_em,    binmin6, binmax6, 15);
    if (!combine2e2mu) {
    smooth(h2D_sig_me,    binmin6, binmax6, 15);
    }
    smooth(h2D_sig_mm,    binmin6, binmax6, 15);
    smooth(h2D_sig_ee,    binmin6, binmax6, 15);
    
    smooth(h2D_bkg_em,    binmin6, binmax6, 15);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me,    binmin6, binmax6, 15);
    }
    smooth(h2D_bkg_mm,    binmin6, binmax6, 15);
    smooth(h2D_bkg_ee,    binmin6, binmax6, 15);

    smooth(h2D_b2g_em,    binmin6, binmax6, 15);
    if (!combine2e2mu) {
    smooth(h2D_b2g_me,    binmin6, binmax6, 15);
    }
    smooth(h2D_b2g_mm,    binmin6, binmax6, 15);
    smooth(h2D_b2g_ee,    binmin6, binmax6, 15);

    smooth(h2D_bkg_em_up, binmin6, binmax6, 15);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_up, binmin6, binmax6, 15);
    }
    smooth(h2D_bkg_mm_up, binmin6, binmax6, 15);
    smooth(h2D_bkg_ee_up, binmin6, binmax6, 15);

    smooth(h2D_bkg_em_dn, binmin6, binmax6, 15);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_dn, binmin6, binmax6, 15);
    }
    smooth(h2D_bkg_mm_dn, binmin6, binmax6, 15);
    smooth(h2D_bkg_ee_dn, binmin6, binmax6, 15);

    smooth(h2D_sig_em,    binmin7, binmax7, 25);
    if (!combine2e2mu) {
    smooth(h2D_sig_me,    binmin7, binmax7, 25);
    }
    smooth(h2D_sig_mm,    binmin7, binmax7, 25);
    smooth(h2D_sig_ee,    binmin7, binmax7, 25);
    
    smooth(h2D_bkg_em,    binmin7, binmax7, 25);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me,    binmin7, binmax7, 25);
    }
    smooth(h2D_bkg_mm,    binmin7, binmax7, 25);
    smooth(h2D_bkg_ee,    binmin7, binmax7, 25);

    smooth(h2D_b2g_em,    binmin7, binmax7, 25);
    if (!combine2e2mu) {
    smooth(h2D_b2g_me,    binmin7, binmax7, 25);
    }
    smooth(h2D_b2g_mm,    binmin7, binmax7, 25);
    smooth(h2D_b2g_ee,    binmin7, binmax7, 25);

    smooth(h2D_bkg_em_up, binmin7, binmax7, 25);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_up, binmin7, binmax7, 25);
    }
    smooth(h2D_bkg_mm_up, binmin7, binmax7, 25);
    smooth(h2D_bkg_ee_up, binmin7, binmax7, 25);

    smooth(h2D_bkg_em_dn, binmin7, binmax7, 25);
    if (!combine2e2mu) {
    smooth(h2D_bkg_me_dn, binmin7, binmax7, 25);
    }
    smooth(h2D_bkg_mm_dn, binmin7, binmax7, 25);
    smooth(h2D_bkg_ee_dn, binmin7, binmax7, 25);


    h2D_sig_em->Smooth();
    if (!combine2e2mu) {
    h2D_sig_me->Smooth();
    }
    h2D_sig_mm->Smooth();
    h2D_sig_ee->Smooth();
    
    h2D_bkg_em->Smooth();
    if (!combine2e2mu) {
    h2D_bkg_me->Smooth();
    }
    h2D_bkg_mm->Smooth();
    h2D_bkg_ee->Smooth();
    
    h2D_b2g_em->Smooth();
    if (!combine2e2mu) {
    h2D_b2g_me->Smooth();
    }
    h2D_b2g_mm->Smooth();
    h2D_b2g_ee->Smooth();
    
    h2D_bkg_em_up->Smooth();
    if (!combine2e2mu) {
    h2D_bkg_me_up->Smooth();
    }
    h2D_bkg_mm_up->Smooth();
    h2D_bkg_ee_up->Smooth();
    
    h2D_bkg_em_dn->Smooth();
    if (!combine2e2mu) {
    h2D_bkg_me_dn->Smooth();
    }
    h2D_bkg_mm_dn->Smooth();
    h2D_bkg_ee_dn->Smooth();

    h2D_sig_em->GetXaxis()->SetRange(41,750);
    h2D_sig_me->GetXaxis()->SetRange(41,750);
    h2D_sig_mm->GetXaxis()->SetRange(41,750);
    h2D_sig_ee->GetXaxis()->SetRange(41,750);

    h2D_bkg_em->GetXaxis()->SetRange(41,750);
    h2D_bkg_me->GetXaxis()->SetRange(41,750);
    h2D_bkg_mm->GetXaxis()->SetRange(41,750);
    h2D_bkg_ee->GetXaxis()->SetRange(41,750);

    h2D_b2g_em->GetXaxis()->SetRange(41,750);
    h2D_b2g_me->GetXaxis()->SetRange(41,750);
    h2D_b2g_mm->GetXaxis()->SetRange(41,750);
    h2D_b2g_ee->GetXaxis()->SetRange(41,750);

    h2D_bkg_em_up->GetXaxis()->SetRange(41,750);
    h2D_bkg_me_up->GetXaxis()->SetRange(41,750);
    h2D_bkg_mm_up->GetXaxis()->SetRange(41,750);
    h2D_bkg_ee_up->GetXaxis()->SetRange(41,750);

    h2D_bkg_em_dn->GetXaxis()->SetRange(41,750);
    h2D_bkg_me_dn->GetXaxis()->SetRange(41,750);
    h2D_bkg_mm_dn->GetXaxis()->SetRange(41,750);
    h2D_bkg_ee_dn->GetXaxis()->SetRange(41,750);

    h2D_sig_em->Smooth();
    if (!combine2e2mu) {
    h2D_sig_me->Smooth();
    }
    h2D_sig_mm->Smooth();
    h2D_sig_ee->Smooth();
    
    h2D_bkg_em->Smooth();
    if (!combine2e2mu) {
    h2D_bkg_me->Smooth();
    }
    h2D_bkg_mm->Smooth();
    h2D_bkg_ee->Smooth();
    
    h2D_b2g_em->Smooth();
    if (!combine2e2mu) {
    h2D_b2g_me->Smooth();
    }
    h2D_b2g_mm->Smooth();
    h2D_b2g_ee->Smooth();
    
    h2D_bkg_em_up->Smooth();
    if (!combine2e2mu) {
    h2D_bkg_me_up->Smooth();
    }
    h2D_bkg_mm_up->Smooth();
    h2D_bkg_ee_up->Smooth();
    
    h2D_bkg_em_dn->Smooth();
    if (!combine2e2mu) {
    h2D_bkg_me_dn->Smooth();
    }
    h2D_bkg_mm_dn->Smooth();
    h2D_bkg_ee_dn->Smooth();

    h2D_sig_em->GetXaxis()->SetRange(1,750);
    h2D_sig_me->GetXaxis()->SetRange(1,750);
    h2D_sig_mm->GetXaxis()->SetRange(1,750);
    h2D_sig_ee->GetXaxis()->SetRange(1,750);

    h2D_bkg_em->GetXaxis()->SetRange(1,750);
    h2D_bkg_me->GetXaxis()->SetRange(1,750);
    h2D_bkg_mm->GetXaxis()->SetRange(1,750);
    h2D_bkg_ee->GetXaxis()->SetRange(1,750);

    h2D_b2g_em->GetXaxis()->SetRange(1,750);
    h2D_b2g_me->GetXaxis()->SetRange(1,750);
    h2D_b2g_mm->GetXaxis()->SetRange(1,750);
    h2D_b2g_ee->GetXaxis()->SetRange(1,750);

    h2D_bkg_em_up->GetXaxis()->SetRange(1,750);
    h2D_bkg_me_up->GetXaxis()->SetRange(1,750);
    h2D_bkg_mm_up->GetXaxis()->SetRange(1,750);
    h2D_bkg_ee_up->GetXaxis()->SetRange(1,750);

    h2D_bkg_em_dn->GetXaxis()->SetRange(1,750);
    h2D_bkg_me_dn->GetXaxis()->SetRange(1,750);
    h2D_bkg_mm_dn->GetXaxis()->SetRange(1,750);
    h2D_bkg_ee_dn->GetXaxis()->SetRange(1,750);

    std::cout << "Applying interference reweighting to signal ...." << std::endl;
    reweightForInterference(h2D_sig_mm);
    reweightForInterference(h2D_sig_ee);

    reweightForInterference(h2D_s2g_mm);
    reweightForInterference(h2D_s2g_ee);

    normalize(h2D_sig_em);
    if (!combine2e2mu) {
    normalize(h2D_sig_me);
    }
    normalize(h2D_sig_mm);
    normalize(h2D_sig_ee);
    
    normalize(h2D_bkg_em);
    if (!combine2e2mu) {
    normalize(h2D_bkg_me);
    }
    normalize(h2D_bkg_mm);
    normalize(h2D_bkg_ee);
    
    normalize(h2D_b2g_em);
    if (!combine2e2mu) {
    normalize(h2D_b2g_me);
    }
    normalize(h2D_b2g_mm);
    normalize(h2D_b2g_ee);
    
    normalize(h2D_bkg_em_up);
    if (!combine2e2mu) {
    normalize(h2D_bkg_me_up);
    }
    normalize(h2D_bkg_mm_up);
    normalize(h2D_bkg_ee_up);
    
    normalize(h2D_bkg_em_dn);
    if (!combine2e2mu) {
    normalize(h2D_bkg_me_dn);
    }
    normalize(h2D_bkg_mm_dn);
    normalize(h2D_bkg_ee_dn);


//////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = 1; i <= h2D_sig_em->GetNbinsX(); i++) {
        for (int j = 1; j <= h2D_sig_em->GetNbinsY(); j++) {
            h2D_s2g_em->SetBinContent(i,j,h2D_sig_em->GetBinContent(i, j));
            h2D_s2g_me->SetBinContent(i,j,h2D_sig_me->GetBinContent(i, j));
            h2D_s2g_ee->SetBinContent(i,j,h2D_sig_ee->GetBinContent(i, j));
            h2D_s2g_mm->SetBinContent(i,j,h2D_sig_mm->GetBinContent(i, j));
        }
    }

    std::cout << "Writing to file ...." << std::endl;
    TFile histout("mela2DShapes.root", "RECREATE");
    h2D_sig_em->Write();
    if (!combine2e2mu) h2D_sig_me->Write();
    h2D_sig_mm->Write();
    h2D_sig_ee->Write();
    h2D_s2g_em->Write();
    if (!combine2e2mu) h2D_s2g_me->Write();
    h2D_s2g_mm->Write();
    h2D_s2g_ee->Write();
    h2D_bkg_em->Write();
    if (!combine2e2mu) h2D_bkg_me->Write();
    h2D_bkg_mm->Write();
    h2D_bkg_ee->Write();
    h2D_b2g_em->Write();
    if (!combine2e2mu) h2D_b2g_me->Write();
    h2D_b2g_mm->Write();
    h2D_b2g_ee->Write();
    h2D_bkg_em_up->Write();
    if (!combine2e2mu) h2D_bkg_me_up->Write();
    h2D_bkg_mm_up->Write();
    h2D_bkg_ee_up->Write();
    h2D_bkg_em_dn->Write();
    if (!combine2e2mu) h2D_bkg_me_dn->Write();
    h2D_bkg_mm_dn->Write();
    h2D_bkg_ee_dn->Write();
    histout.Close();
}
