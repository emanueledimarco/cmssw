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

float melacut           = 0.5;
float z1min             = 40.;
float z2min             = 12.;
float massLow           = 100.;
float massHigh          = 180.;
int   nBinsMass2D       = 40; //2GeV bins
int   nBinsMELA2D       = 30.; 
bool  do7TeV            = true; 
std::string treeFolder7 = "/data1/sixie/ntuples/step3new/hcp/2011/";
std::string treeFolder8 = "/data1/sixie/ntuples/step3new/hcp/2012/";

void fillMass(std::string rootfile, float xsecsf, int hm, bool i7, TH2F& h2D_sig_em, TH2F& h2D_sig_mm, TH2F& h2D_sig_ee) {
    std::string base_folder = i7 ? treeFolder7 : treeFolder8;
    ZZYieldMaker   ymaker;

    std::cout << "7TeV? : " << do7TeV << " : " << base_folder << std::endl;
    std::cout << "Filling " << base_folder+rootfile <<  " with scale factor " << xsecsf << std::endl;

    ymaker.fill(base_folder+rootfile, xsecsf, 0.0, true, 1, hm, i7 );

    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_mm);
    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_ee);
    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);
    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);

}

void fillBkg(TH2F& h2D_bkg_em, TH2F& h2D_bkg_mm, TH2F& h2D_bkg_ee, bool isGG) {
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
        ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
        ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
        ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
        ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);       
    }
    else {
        ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
        ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
        ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
        ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);        
    }

}

void fillZjets(TH2F& h2D_bkg_em, TH2F& h2D_bkg_mm, TH2F& h2D_bkg_ee, bool doSS) {
    std::string base_folder = do7TeV ? treeFolder7 : treeFolder8;

    FakeRateCalculator FR(base_folder+"hcp.root", do7TeV, 40, 120, 0.0, 0.0, true);

    ZXYieldMaker   ymaker_zxss;
    ymaker_zxss.fill(base_folder+"hcp.root", 1.0, FR, doSS);

    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);

    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);

}


///// This code is taken from the HZZ4L_Combination package, since we need to rely on the MELA experts for this.

void reweightForInterference(TH2F& hist){
    TF1* reweightFunc =new TF1("reweightFunc","gaus",100,1600);
      
    reweightFunc->SetParameter(0,0.4053);
    reweightFunc->SetParameter(1,110.6);
    reweightFunc->SetParameter(2,22.57);
    
    double oldTempValue=0;
    double newTempValue=0;
    double slope;
    
    for(int i=1; i<=hist.GetNbinsX(); i++){
        for(int j=1; j<=hist.GetNbinsY(); j++){
            slope=reweightFunc->Eval((double)((i-1)*2+101));
            oldTempValue = hist.GetBinContent(i,j);
            newTempValue = oldTempValue*(1+slope*((double)j/30.-.5));
            hist.SetBinContent(i,j,newTempValue);
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

void createSpinTwoMinimalMelaVsMassHistogram() {
    TH2F h2D_sig_em("hist2D_sig_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_sig_mm("hist2D_sig_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_sig_ee("hist2D_sig_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    TH2F h2D_spin2minimalsig_em("hist2D_spin2minimalsig_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_spin2minimalsig_mm("hist2D_spin2minimalsig_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_spin2minimalsig_ee("hist2D_spin2minimalsig_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    TH2F h2D_bkg_em("hist2D_bkg_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_mm("hist2D_bkg_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_ee("hist2D_bkg_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    TH2F h2D_b2g_em("hist2D_b2g_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_b2g_mm("hist2D_b2g_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_b2g_ee("hist2D_b2g_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);


    h2D_sig_em.Sumw2();
    h2D_sig_mm.Sumw2();
    h2D_sig_ee.Sumw2();
    
    h2D_spin2minimalsig_em.Sumw2();
    h2D_spin2minimalsig_mm.Sumw2();
    h2D_spin2minimalsig_ee.Sumw2();
    
    h2D_bkg_em.Sumw2();
    h2D_bkg_mm.Sumw2();
    h2D_bkg_ee.Sumw2();
    
    h2D_b2g_em.Sumw2();
    h2D_b2g_mm.Sumw2();
    h2D_b2g_ee.Sumw2();

    //***********************************************************************************************
    //Fill histogram for scalar higgs
    //***********************************************************************************************

    bool do7TeV = true;
    init(do7TeV);

    std::cout << "Filling background" << std::endl;
    fillBkg (h2D_bkg_em, h2D_bkg_mm, h2D_bkg_ee, false);
    fillBkg (h2D_b2g_em, h2D_b2g_mm, h2D_b2g_ee, true);

    do7TeV = false;
    init(do7TeV);

    fillMass("hzzTree_id1125.root" , 1.0,  125, kFALSE, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id9125.root" , 1.0 , 125, kFALSE, h2D_spin2minimalsig_em, h2D_spin2minimalsig_mm, h2D_spin2minimalsig_ee);

    std::cout << "Filling background" << std::endl;
    fillBkg (h2D_bkg_em, h2D_bkg_mm, h2D_bkg_ee, false);
    fillBkg (h2D_b2g_em, h2D_b2g_mm, h2D_b2g_ee, true);

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating alternate shapes ...." << std::endl;
    TH2F h2D_bkg_em_up = *((TH2F*)h2D_bkg_em.Clone("hist2D_bkg_em_up"));
    TH2F h2D_bkg_mm_up = *((TH2F*)h2D_bkg_mm.Clone("hist2D_bkg_mm_up"));
    TH2F h2D_bkg_ee_up = *((TH2F*)h2D_bkg_ee.Clone("hist2D_bkg_ee_up"));

    TH2F h2D_bkg_em_dn = *((TH2F*)h2D_bkg_em.Clone("hist2D_bkg_em_dn"));
    TH2F h2D_bkg_mm_dn = *((TH2F*)h2D_bkg_mm.Clone("hist2D_bkg_mm_dn"));
    TH2F h2D_bkg_ee_dn = *((TH2F*)h2D_bkg_ee.Clone("hist2D_bkg_ee_dn"));

    nBinsMass2D = 6;
    Double_t systbins[] = {100., 120., 140., 160., 180., 300., 1600.};
    TH2F h2D_zz_em("hist2D_zz_em", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zz_mm("hist2D_zz_mm", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zz_ee("hist2D_zz_ee", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);

    TH2F h2D_zx_em("hist2D_zx_em", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zx_mm("hist2D_zx_mm", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zx_ee("hist2D_zx_ee", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);

    h2D_zx_em.Sumw2();
    h2D_zx_mm.Sumw2();
    h2D_zx_ee.Sumw2();

    do7TeV = true;
    init(do7TeV);
    fillBkg  (h2D_zz_em, h2D_zz_mm, h2D_zz_ee, false);
    fillZjets(h2D_zx_em, h2D_zx_mm, h2D_zx_ee, true);

    do7TeV = false;
    init(do7TeV);
    fillBkg  (h2D_zz_em, h2D_zz_mm, h2D_zz_ee, false);
    fillZjets(h2D_zx_em, h2D_zx_mm, h2D_zx_ee, true);

    h2D_zz_em.Add(&h2D_zz_mm);
    h2D_zz_em.Add(&h2D_zz_ee);
    h2D_zx_em.Add(&h2D_zx_mm);
    h2D_zx_em.Add(&h2D_zx_ee);
    normalize(&h2D_zz_em, &h2D_zx_em);

    std::vector<float> fitp0;
    std::vector<float> fitp1;
    for (int i = 1; i <= nBinsMass2D; i++) {
        stringstream   histname;
        histname   << "hist";
        histname   << i;
        TH1F hist(histname.str().c_str()  , "", nBinsMELA2D, 0., 1.);
        TF1 f("f", "[0]+[1]*x", 0, 1600.);
        for (int j = 1; j <= nBinsMELA2D; j++) {
            hist.SetBinContent  (j, h2D_zx_em.GetBinContent(i,j)/h2D_zz_em.GetBinContent(i,j));
            hist.SetBinError    (j, h2D_zx_em.GetBinError  (i,j)/h2D_zz_em.GetBinContent(i,j));
        }
        hist.Fit(&f);
        fitp0.push_back(f.GetParameter(0));
        fitp1.push_back(f.GetParameter(1));
    }

    for (int i = 1; i <= h2D_bkg_em.GetNbinsX(); i++) {
        for (int j = 1; j <= h2D_bkg_em.GetNbinsY(); j++) {
            int k = nBinsMass2D;
            for (int m = 0; m < nBinsMass2D; m++) {
                if ((massLow + i*2) > systbins[m] && (massLow + i*2) <= systbins[m+1]) {
                    k = m;   
                }       
            }

            if (k < nBinsMass2D) {

                h2D_bkg_em_up.SetBinContent(i, j, h2D_bkg_em.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_ee_up.SetBinContent(i, j, h2D_bkg_ee.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_mm_up.SetBinContent(i, j, h2D_bkg_mm.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));

                h2D_bkg_em_dn.SetBinContent(i, j, h2D_bkg_em.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_ee_dn.SetBinContent(i, j, h2D_bkg_ee.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_mm_dn.SetBinContent(i, j, h2D_bkg_mm.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
            }
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Smoothing ...." << std::endl;
    normalize(&h2D_sig_em);
    normalize(&h2D_sig_mm);
    normalize(&h2D_sig_ee);

    normalize(&h2D_spin2minimalsig_em);
    normalize(&h2D_spin2minimalsig_mm);
    normalize(&h2D_spin2minimalsig_ee);

    normalize(&h2D_bkg_em);
    normalize(&h2D_bkg_mm);
    normalize(&h2D_bkg_ee);

    normalize(&h2D_b2g_em);
    normalize(&h2D_b2g_mm);
    normalize(&h2D_b2g_ee);

    normalize(&h2D_bkg_em_up);
    normalize(&h2D_bkg_mm_up);
    normalize(&h2D_bkg_ee_up);

    normalize(&h2D_bkg_em_dn);
    normalize(&h2D_bkg_mm_dn);
    normalize(&h2D_bkg_ee_dn);

    //just use smoothing of TH2
    h2D_sig_em.Smooth();
    h2D_sig_mm.Smooth();
    h2D_sig_ee.Smooth();
    
    h2D_spin2minimalsig_em.Smooth();
    h2D_spin2minimalsig_mm.Smooth();
    h2D_spin2minimalsig_ee.Smooth();
    
    h2D_bkg_em.Smooth();
    h2D_bkg_mm.Smooth();
    h2D_bkg_ee.Smooth();
    
    h2D_b2g_em.Smooth();
    h2D_b2g_mm.Smooth();
    h2D_b2g_ee.Smooth();
    
    h2D_bkg_em_up.Smooth();
    h2D_bkg_mm_up.Smooth();
    h2D_bkg_ee_up.Smooth();
    
    h2D_bkg_em_dn.Smooth();
    h2D_bkg_mm_dn.Smooth();
    h2D_bkg_ee_dn.Smooth();

//No Reweighting for Spin2
//     reweightForInterference(h2D_sig_em);
//     reweightForInterference(h2D_sig_mm);
//     reweightForInterference(h2D_sig_ee);

//     reweightForInterference(h2D_spin2minimalsig_em);
//     reweightForInterference(h2D_spin2minimalsig_mm);
//     reweightForInterference(h2D_spin2minimalsig_ee);

    normalize(&h2D_sig_em);
    normalize(&h2D_sig_mm);
    normalize(&h2D_sig_ee);
    
    normalize(&h2D_spin2minimalsig_em);
    normalize(&h2D_spin2minimalsig_mm);
    normalize(&h2D_spin2minimalsig_ee);
    
    normalize(&h2D_bkg_em);
    normalize(&h2D_bkg_mm);
    normalize(&h2D_bkg_ee);
    
    normalize(&h2D_b2g_em);
    normalize(&h2D_b2g_mm);
    normalize(&h2D_b2g_ee);
    
    normalize(&h2D_bkg_em_up);
    normalize(&h2D_bkg_mm_up);
    normalize(&h2D_bkg_ee_up);
    
    normalize(&h2D_bkg_em_dn);
    normalize(&h2D_bkg_mm_dn);
    normalize(&h2D_bkg_ee_dn);


//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Writing to file ...." << std::endl;
    TFile histout("SpinTwoMinimalMelaVsMass2DShapes.root", "RECREATE");
    h2D_sig_em.Write();
    h2D_sig_mm.Write();
    h2D_sig_ee.Write();
    h2D_spin2minimalsig_em.Write();
    h2D_spin2minimalsig_mm.Write();
    h2D_spin2minimalsig_ee.Write();
    h2D_bkg_em.Write();
    h2D_bkg_mm.Write();
    h2D_bkg_ee.Write();
    h2D_b2g_em.Write();
    h2D_b2g_mm.Write();
    h2D_b2g_ee.Write();
    h2D_bkg_em_up.Write();
    h2D_bkg_mm_up.Write();
    h2D_bkg_ee_up.Write();
    h2D_bkg_em_dn.Write();
    h2D_bkg_mm_dn.Write();
    h2D_bkg_ee_dn.Write();
    histout.Close();
}
