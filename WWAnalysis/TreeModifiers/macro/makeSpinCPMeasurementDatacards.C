#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooLognormal.h> 
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>
#include <RooHistFunc.h>

#include "CardTemplate.h"
#include "FitMaker.h"
#include "YieldMaker.h"
#include "findAndReplace.h"
#include "scales2.h"
#include "SignalInterpolationStrings.h"

using namespace RooFit;

struct HiggsMassPointInfo {

    float lumi;
    float base_lumi;
    float melacut;
    float z1min;
    float z2min;
    float massLowBkgFit;
    float massHighBkgFit;
    bool  do1D;
    bool  do7TeV;
    bool  doFFT;
    bool  doMassError;
    int   alternativeSignalType;
    std::string treeFolder;
    std::string melafilename;
    std::string alternativeMELAfilename;
    
    DataYieldMaker ymaker_data;


    void createCard(float mass, float massLow, float massHigh, int ch) {

        std::string chstr;
        if (ch == 0) chstr = "4mu";
        if (ch == 1) chstr = "4e";
        if (ch == 2) chstr = "2e2mu";

        stringstream mass_str_ss;
        mass_str_ss << mass;
        std::string mass_str = mass_str_ss.str();
        
        std::cout << "Creating datacard for " << mass_str << " GeV mass point and channel " << chstr << " ... " << std::endl;
        
        std::string card_name = "";
        if (alternativeSignalType == 1) card_name   = (std::string("card_pseudoscalar_2D_m")+mass_str+"_"+(do7TeV ? "7TeV_" : "8TeV_"));
        else if (alternativeSignalType == 2) card_name   = (std::string("card_spinTwoMinimal_2D_m")+mass_str+"_"+(do7TeV ? "7TeV_" : "8TeV_"));
        else { cout << "Error. signal type = " <<  alternativeSignalType << " not supported\n"; return; }
        card_name += chstr;
        std::string workspace = card_name+"_workspace.root";
        
        float yield_dt  = ymaker_data.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
        if (ch == 2) yield_dt += ymaker_data.getYield(3, z1min, z2min, massLow, massHigh, melacut);
         
        RooWorkspace w("w", "");
        
        RooRealVar CMS_zz4l_melaLD ("CMS_zz4l_melaLD" , "MELA"       ,   0    , 1       , "");
        RooRealVar CMS_zz4l_AlternativeMELA ("CMS_zz4l_AlternativeMELA" , "AlternativeMELA"       ,   0    , 1       , "");
        RooRealVar CMS_zz4l_mass_1D("CMS_zz4l_mass_1D", "M(4l)"      , massLow, massHigh, "GeV/c^{2}");
        RooRealVar CMS_zz4l_massErr("CMS_zz4l_massErr", "M(4l) Error",   0    , 5       , "GeV/c^{2}");
        CMS_zz4l_massErr.setBins(25);
        
        if (doFFT) CMS_zz4l_mass_1D.setBins(100000, "fft");
        
        RooArgSet argset_obs(CMS_zz4l_mass_1D, CMS_zz4l_melaLD, CMS_zz4l_AlternativeMELA, "argset_obs");
        RooDataSet data_obs("data_obs", "data_obs", argset_obs);
        
        if (alternativeSignalType == 1) {
          ymaker_data.getDataSetMassVsPSMela(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_melaLD, CMS_zz4l_AlternativeMELA);
          if (ch == 2) ymaker_data.getDataSetMassVsPSMela(3, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_melaLD, CMS_zz4l_AlternativeMELA);
        } else if (alternativeSignalType == 2) {
          ymaker_data.getDataSetMassVsSpinTwoMinimalMela(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_melaLD, CMS_zz4l_AlternativeMELA);
          if (ch == 2) ymaker_data.getDataSetMassVsSpinTwoMinimalMela(3, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_melaLD, CMS_zz4l_AlternativeMELA);
        }
        
        w.import(data_obs);
         
        ///////////////////// Define parameters //////////////////////////////////

        float qqa0  = 0.;
        float qqa1  = 0.;
        float qqa2  = 0.;
        float qqa3  = 0.;
        float qqa4  = 0.;
        float qqa5  = 0.;
        float qqa6  = 0.;
        float qqa7  = 0.;
        float qqa8  = 0.;
        float qqa9  = 0.;
        float qqa10 = 0.;
        float qqa11 = 0.;
        float qqa12 = 0.;
        float qqa13 = 0.;

        float gga0  = 0.;
        float gga1  = 0.;
        float gga2  = 0.;
        float gga3  = 0.;
        float gga4  = 0.;
        float gga5  = 0.;
        float gga6  = 0.;
        float gga7  = 0.;
        float gga8  = 0.;
        float gga9  = 0.;

        float zxme  = 0.;
        float zxsi  = 0.;

        if (do7TeV) {

            if (ch == 0) {
                zxme  = 141.2;
                zxsi = 19.1;


                qqa0  = 107.997;
                qqa1  = 5.30551;
                qqa2  = 116.325;
                qqa3  = 0.0435405;
                qqa4  = 185.792;
                qqa5  = 8.42431;
                qqa6  = 35.3956;
                qqa7  = 0.0992744;
                qqa8  = 57.4238;
                qqa9  = 0.0498783;
                qqa10 = 100.267;
                qqa11 = -8.48533;
                qqa12 = 1240.17;
                qqa13 = 0.0863839;


                gga0  = 145.053;
                gga1  = 59.4322;
                gga2  = 120.157;
                gga3  = 0.0448841;
                gga4  = 183.875;
                gga5  = 9.01533;
                gga6  = 41.0342;
                gga7  = 0.466347;
                gga8  = 37.616;
                gga9  = -0.216485;

            }

            else if (ch == 1) {
                zxme  = 156.1;
                zxsi = 23.5;


                qqa0  = 115.07;
                qqa1  = 24.2583;
                qqa2  = 115.911;
                qqa3  = 0.0423577;
                qqa4  = 183.922;
                qqa5  = 11.0038;
                qqa6  = 35.0579;
                qqa7  = 0.104589;
                qqa8  = 55.9744;
                qqa9  = 0.0509771;
                qqa10 = 94.338;
                qqa11 = -6.95537;
                qqa12 = 1169.59;
                qqa13 = 0.138667;


                gga0  = 148.362;
                gga1  = 49.0322;
                gga2  = 124.348;
                gga3  = 0.0445723;
                gga4  = 183.049;
                gga5  = 11.4454;
                gga6  = 41.7409;
                gga7  = 0.495309;
                gga8  = 40.5164;
                gga9  = -0.198228;

            }

            else {
                zxme  = 153.7;
                zxsi = 23.3;


                qqa0  = 109.684;
                qqa1  = 17.3937;
                qqa2  = 116.859;
                qqa3  = 0.04159;
                qqa4  = 184.657;
                qqa5  = 9.8782;
                qqa6  = 36.0073;
                qqa7  = 0.105264;
                qqa8  = 57.6245;
                qqa9  = 0.0541577;
                qqa10 = 95.7085;
                qqa11 = -5.85807;
                qqa12 = 619.854;
                qqa13 = 0.0613803;


                gga0  = 155.967;
                gga1  = 51.3515;
                gga2  = 120.248;
                gga3  = 0.050175;
                gga4  = 184.576;
                gga5  = 11.5281;
                gga6  = 42.8414;
                gga7  = 0.478245;
                gga8  = 44.629;
                gga9  = -0.180257;

            }

        }

        else {
            if (ch == 0) {
                zxme  = 147.1;
                zxsi = 32.8;


                qqa0  = 107.085;
                qqa1  = 10.4552;
                qqa2  = 118.199;
                qqa3  = 0.045443;
                qqa4  = 185.568;
                qqa5  = 8.86547;
                qqa6  = 38.1948;
                qqa7  = 0.0986066;
                qqa8  = 52.993;
                qqa9  = 0.047083;
                qqa10 = 100.71;
                qqa11 = -9.74852;
                qqa12 = 2359.79;
                qqa13 = 0.0566742;


                gga0  = 129.229;
                gga1  = 101.025;
                gga2  = 133.144;
                gga3  = 0.0488151;
                gga4  = 185.503;
                gga5  = 9.5927;
                gga6  = 44.0617;
                gga7  = 0.471199;
                gga8  = 42.2465;
                gga9  = -0.175023;

            }

            else if (ch == 1) {
                zxme  = 158.6;
                zxsi = 23.5;


                qqa0  = 114.809;
                qqa1  = 22.0457;
                qqa2  = 121.635;
                qqa3  = 0.0383133;
                qqa4  = 184.887;
                qqa5  = 10.87;
                qqa6  = 34.6625;
                qqa7  = 0.115359;
                qqa8  = 61.2978;
                qqa9  = 0.0631035;
                qqa10 = 97.7677;
                qqa11 = -8.25902;
                qqa12 = 2772.12;
                qqa13 = 0.069981;


                gga0  = 140.519;
                gga1  = 57.6559;
                gga2  = 134.049;
                gga3  = 0.0462216;
                gga4  = 184.591;
                gga5  = 9.97182;
                gga6  = 42.9434;
                gga7  = 0.476164;
                gga8  = 39.892;
                gga9  = -0.173694;

            }
        
            else {
                zxme  = 154.8;
                zxsi = 24.0;


                qqa0  = 108.602;
                qqa1  = 15.2791;
                qqa2  = 122.349;
                qqa3  = 0.0400475;
                qqa4  = 185.232;
                qqa5  = 9.46401;
                qqa6  = 34.4605;
                qqa7  = 0.110499;
                qqa8  = 61.0546;
                qqa9  = 0.059151;
                qqa10 = 95.2315;
                qqa11 = -6.48344;
                qqa12 = 1665.9;
                qqa13 = 0.095991;


                gga0  = 164.834;
                gga1  = 60.5107;
                gga2  = 121.166;
                gga3  = 0.0506449;
                gga4  = 185.414;
                gga5  = 9.31951;
                gga6  = 43.4125;
                gga7  = 0.462093;
                gga8  = 43.9754;
                gga9  = -0.18577;


            }
        }



        std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
        stringstream lumiss;
        lumiss << lumi;
        std::string lumistr = lumiss.str(); 
        RooRealVar qqzz_a0 (("bkg_qqzz_"+chstr+tevstr+"_a0" ).c_str(), "", qqa0 );
        RooRealVar qqzz_a1 (("bkg_qqzz_"+chstr+tevstr+"_a1" ).c_str(), "", qqa1 );
        RooRealVar qqzz_a2 (("bkg_qqzz_"+chstr+tevstr+"_a2" ).c_str(), "", qqa2 );
        RooRealVar qqzz_a3 (("bkg_qqzz_"+chstr+tevstr+"_a3" ).c_str(), "", qqa3 );
        RooRealVar qqzz_a4 (("bkg_qqzz_"+chstr+tevstr+"_a4" ).c_str(), "", qqa4 );
        RooRealVar qqzz_a5 (("bkg_qqzz_"+chstr+tevstr+"_a5" ).c_str(), "", qqa5 );
        RooRealVar qqzz_a6 (("bkg_qqzz_"+chstr+tevstr+"_a6" ).c_str(), "", qqa6 );
        RooRealVar qqzz_a7 (("bkg_qqzz_"+chstr+tevstr+"_a7" ).c_str(), "", qqa7 );
        RooRealVar qqzz_a8 (("bkg_qqzz_"+chstr+tevstr+"_a8" ).c_str(), "", qqa8 );
        RooRealVar qqzz_a9 (("bkg_qqzz_"+chstr+tevstr+"_a9" ).c_str(), "", qqa9 );
        RooRealVar qqzz_a10(("bkg_qqzz_"+chstr+tevstr+"_a10").c_str(), "", qqa10);
        RooRealVar qqzz_a11(("bkg_qqzz_"+chstr+tevstr+"_a11").c_str(), "", qqa11);
        RooRealVar qqzz_a12(("bkg_qqzz_"+chstr+tevstr+"_a12").c_str(), "", qqa12);
        RooRealVar qqzz_a13(("bkg_qqzz_"+chstr+tevstr+"_a13").c_str(), "", qqa13);
        
        RooRealVar ggzz_a0 (("bkg_ggzz_"+chstr+tevstr+"_a0" ).c_str(), "", gga0 );
        RooRealVar ggzz_a1 (("bkg_ggzz_"+chstr+tevstr+"_a1" ).c_str(), "", gga1 );
        RooRealVar ggzz_a2 (("bkg_ggzz_"+chstr+tevstr+"_a2" ).c_str(), "", gga2 );
        RooRealVar ggzz_a3 (("bkg_ggzz_"+chstr+tevstr+"_a3" ).c_str(), "", gga3 );
        RooRealVar ggzz_a4 (("bkg_ggzz_"+chstr+tevstr+"_a4" ).c_str(), "", gga4 );
        RooRealVar ggzz_a5 (("bkg_ggzz_"+chstr+tevstr+"_a5" ).c_str(), "", gga5 );
        RooRealVar ggzz_a6 (("bkg_ggzz_"+chstr+tevstr+"_a6" ).c_str(), "", gga6 );
        RooRealVar ggzz_a7 (("bkg_ggzz_"+chstr+tevstr+"_a7" ).c_str(), "", gga7 );
        RooRealVar ggzz_a8 (("bkg_ggzz_"+chstr+tevstr+"_a8" ).c_str(), "", gga8 );
        RooRealVar ggzz_a9 (("bkg_ggzz_"+chstr+tevstr+"_a9" ).c_str(), "", gga9 );
        
        RooRealVar zx_mean (("bkg_zjets_"+chstr+tevstr+"_mean_zx" ).c_str(), "", zxme);
        RooRealVar zx_sigma(("bkg_zjets_"+chstr+tevstr+"_sigma_zx").c_str(), "", zxsi);
        
        RooRealVar sig_mean_err    (("sig_"+chstr+"_mean_err"  +tevstr).c_str()  , "", 0., -10., 10.);
        RooRealVar sig_gamma_err   ("sig_gamma_err"                              , "", 0., -10., 10.);
        RooRealVar sig_sigma_err   (("sig_"+chstr+"_sigma_err" +tevstr).c_str()  , "", 0., -10., 10.);

        RooRealVar qqzz_sigma_EBK(("bkg_qqzz_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);
        RooRealVar ggzz_sigma_EBK(("bkg_ggzz_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);
        RooRealVar zjet_sigma_EBK(("bkg_zjet_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);

        RooFormulaVar qqzz_sigma_EBM(("bkg_qqzz_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(CMS_zz4l_mass_1D));
        RooFormulaVar ggzz_sigma_EBM(("bkg_ggzz_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(CMS_zz4l_mass_1D));
        RooFormulaVar zjet_sigma_EBM(("bkg_zjet_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(CMS_zz4l_mass_1D));

        RooRealVar masshiggs       ("MH", "", mass);

        std::string cs_scale_str = "";
        if (do7TeV) cs_scale_str += "0.5+0.5*TMath::Erf((@0 - 80.85)/50.42)";
        else cs_scale_str += "1";    
        RooFormulaVar cs_scale_z2_ggh  (("cs_scale_z2_ggH"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_vbf  (("cs_scale_z2_VBF"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_whi  (("cs_scale_z2_WHi"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_zhi  (("cs_scale_z2_ZHi"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_tth  (("cs_scale_z2_ttH"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar yield_var_ggh    (("yield_eff_ggH_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_vbf    (("yield_eff_VBF_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_whi    (("yield_eff_WHi_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_zhi    (("yield_eff_ZHi_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_tth    (("yield_eff_ttH_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));

        RooFormulaVar cs_scale_z2_gghA (("cs_scale_z2_ggH_ALT"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_vbfA (("cs_scale_z2_VBF_ALT"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_whiA (("cs_scale_z2_WHi_ALT"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_zhiA (("cs_scale_z2_ZHi_ALT"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_tthA (("cs_scale_z2_ttH_ALT"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar yield_var_gghA   (("yield_eff_ggH_ALT_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_vbfA   (("yield_eff_VBF_ALT_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_whiA   (("yield_eff_WHi_ALT_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_zhiA   (("yield_eff_ZHi_ALT_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_tthA   (("yield_eff_ttH_ALT_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));

        TH1F* histxsecbrggh;
        TH1F* histxsecbrvbf;
        TH1F* histxsecbrwhi;
        TH1F* histxsecbrzhi;
        TH1F* histxsecbrtth;

        if (mass >= 110 && mass < 160) {
            Float_t gghxsecbry[50];
            Float_t vbfxsecbry[50];
            Float_t whixsecbry[50];
            Float_t zhixsecbry[50];
            Float_t tthxsecbry[50];
        
        
            for (float i = 110.; i < 160.; i = i+1.) {
                gghxsecbry[(int)(i-110.)] = getXsecggHByChannel(i, ch);
                vbfxsecbry[(int)(i-110.)] = getXsecVBFByChannel(i, ch);
                whixsecbry[(int)(i-110.)] = getXsecWHiByChannel(i, ch);
                zhixsecbry[(int)(i-110.)] = getXsecZHiByChannel(i, ch);
                tthxsecbry[(int)(i-110.)] = getXsecttHByChannel(i, ch);
            }
        
            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 50, 110., 160.);
        
            for (int i = 1; i <= 50; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 50; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 50; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 50; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 50; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }


        else if (mass >= 160 && mass < 290) {
            Float_t gghxsecbry[65];
            Float_t vbfxsecbry[65];
            Float_t whixsecbry[65];
            Float_t zhixsecbry[65];
            Float_t tthxsecbry[65];


            for (float i = 160.; i < 290.; i = i+2.) {
                gghxsecbry[((int)(i-160.))/2] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-160.))/2] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-160.))/2] = getXsecWHiByChannel(i, ch);
                zhixsecbry[((int)(i-160.))/2] = getXsecZHiByChannel(i, ch);
                tthxsecbry[((int)(i-160.))/2] = getXsecttHByChannel(i, ch);
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 65, 160., 290.);

            for (int i = 1; i <= 65; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }

        else if (mass >= 290 && mass < 350) {
            Float_t gghxsecbry[12];
            Float_t vbfxsecbry[12];
            Float_t whixsecbry[12];
            Float_t zhixsecbry[12];
            Float_t tthxsecbry[12];


            for (float i = 290.; i < 350.; i = i+5.) {
                gghxsecbry[((int)(i-290.))/5] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-290.))/5] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecWHiByChannel(i, ch) : 0.0;
                zhixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecZHiByChannel(i, ch) : 0.0;
                tthxsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecttHByChannel(i, ch) : 0.0;
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 12, 290., 350);

            for (int i = 1; i <= 12; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }

        else if (mass >= 350 && mass < 400) {
            Float_t gghxsecbry[5];
            Float_t vbfxsecbry[5];
            Float_t whixsecbry[5];
            Float_t zhixsecbry[5];
            Float_t tthxsecbry[5];


            for (float i = 350.; i < 400.; i = i+10.) {
                gghxsecbry[((int)(i-350.))/10] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-350.))/10] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-350.))/10] = 0.0;
                zhixsecbry[((int)(i-350.))/10] = 0.0;
                tthxsecbry[((int)(i-350.))/10] = 0.0;
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 5, 350., 400.);

            for (int i = 1; i <= 5; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }

        else  {
            Float_t gghxsecbry[30];
            Float_t vbfxsecbry[30];
            Float_t whixsecbry[30];
            Float_t zhixsecbry[30];
            Float_t tthxsecbry[30];


            for (float i = 400.; i < 1000.; i = i+20.) {
                gghxsecbry[((int)(i-400.))/20] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-400.))/20] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-400.))/20] = 0.;
                zhixsecbry[((int)(i-400.))/20] = 0.;
                tthxsecbry[((int)(i-400.))/20] = 0.;
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);

            for (int i = 1; i <= 30; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }


        RooDataHist dhistxsecbrggh(("rdhistxsecbr_ggH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrggh);
        RooDataHist dhistxsecbrvbf(("rdhistxsecbr_VBF_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrvbf);
        RooDataHist dhistxsecbrwhi(("rdhistxsecbr_WHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrwhi);
        RooDataHist dhistxsecbrzhi(("rdhistxsecbr_ZHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrzhi);
        RooDataHist dhistxsecbrtth(("rdhistxsecbr_ttH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrtth);

        RooHistFunc ggh_xsecbr(("xsecbr_ggH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrggh, 1);
        RooHistFunc vbf_xsecbr(("xsecbr_VBF_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrvbf, 1);
        RooHistFunc whi_xsecbr(("xsecbr_WHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrwhi, 1);
        RooHistFunc zhi_xsecbr(("xsecbr_ZHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrzhi, 1);
        RooHistFunc tth_xsecbr(("xsecbr_ttH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrtth, 1);

        RooHistFunc gghA_xsecbr(("xsecbr_ggH_ALT_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrggh, 1);
        RooHistFunc vbfA_xsecbr(("xsecbr_VBF_ALT_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrvbf, 1);
        RooHistFunc whiA_xsecbr(("xsecbr_WHi_ALT_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrwhi, 1);
        RooHistFunc zhiA_xsecbr(("xsecbr_ZHi_ALT_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrzhi, 1);
        RooHistFunc tthA_xsecbr(("xsecbr_ttH_ALT_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrtth, 1);

        RooFormulaVar ggh_mean_CB  (("sig_ggh_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch, do7TeV, doFFT).c_str()           , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar ggh_sigma_CB (("sig_ggh_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str()                  , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar ggh_alphaL   (("sig_ggh_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str()                            , RooArgList(masshiggs));
        RooFormulaVar ggh_alphaR   (("sig_ggh_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str()                            , RooArgList(masshiggs));
        RooFormulaVar ggh_nL       (("sig_ggh_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()                                , RooArgList(masshiggs));
        RooFormulaVar ggh_nR       (("sig_ggh_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()                                , RooArgList(masshiggs));
        RooFormulaVar ggh_gamma_BW (("sig_ggh_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str()                             , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar ggh_sigma_EBE(("sig_ggh_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                                                          , RooArgList(CMS_zz4l_massErr, ggh_sigma_CB));
        RooFormulaVar ggh_sigma_EBM(("sig_ggh_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar ggh_norm     ("ggH_norm"    , ("@0*@1*@2*"+lumistr).c_str(), RooArgList(cs_scale_z2_ggh, ggh_xsecbr, yield_var_ggh));

        RooFormulaVar gghA_mean_CB  (("sig_ggh_ALT_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch, do7TeV, doFFT).c_str()      , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar gghA_sigma_CB (("sig_ggh_ALT_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str()             , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar gghA_alphaL   (("sig_ggh_ALT_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str()                       , RooArgList(masshiggs));
        RooFormulaVar gghA_alphaR   (("sig_ggh_ALT_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str()                       , RooArgList(masshiggs));
        RooFormulaVar gghA_nL       (("sig_ggh_ALT_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()                           , RooArgList(masshiggs));
        RooFormulaVar gghA_nR       (("sig_ggh_ALT_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()                           , RooArgList(masshiggs));
        RooFormulaVar gghA_gamma_BW (("sig_ggh_ALT_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str()                        , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar gghA_sigma_EBE(("sig_ggh_ALT_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                                                          , RooArgList(CMS_zz4l_massErr, ggh_sigma_CB));
        RooFormulaVar gghA_sigma_EBM(("sig_ggh_ALT_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar gghA_norm     ("ggH_ALT_norm", ("@0*@1*@2*"+lumistr).c_str(), RooArgList(cs_scale_z2_gghA, gghA_xsecbr, yield_var_gghA));
        
        RooFormulaVar vbf_mean_CB  (("sig_vbf_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar vbf_sigma_CB (("sig_vbf_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar vbf_alphaL   (("sig_vbf_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar vbf_alphaR   (("sig_vbf_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar vbf_nL       (("sig_vbf_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar vbf_nR       (("sig_vbf_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar vbf_gamma_BW (("sig_vbf_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar vbf_sigma_EBE(("sig_vbf_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, vbf_sigma_CB));
        RooFormulaVar vbf_sigma_EBM(("sig_vbf_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar vbf_norm     ("qqH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_vbf, vbf_xsecbr, yield_var_vbf));

        RooFormulaVar vbfA_mean_CB  (("sig_vbf_ALT_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar vbfA_sigma_CB (("sig_vbf_ALT_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar vbfA_alphaL   (("sig_vbf_ALT_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar vbfA_alphaR   (("sig_vbf_ALT_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar vbfA_nL       (("sig_vbf_ALT_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar vbfA_nR       (("sig_vbf_ALT_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar vbfA_gamma_BW (("sig_vbf_ALT_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar vbfA_sigma_EBE(("sig_vbf_ALT_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, vbf_sigma_CB));
        RooFormulaVar vbfA_sigma_EBM(("sig_vbf_ALT_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar vbfA_norm     ("qqH_ALT_norm"                            , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_vbfA, vbfA_xsecbr, yield_var_vbfA));
      
        RooFormulaVar whi_mean_CB  (("sig_whi_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch, do7TeV, doFFT).c_str() , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar whi_sigma_CB (("sig_whi_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar whi_alphaL   (("sig_whi_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar whi_alphaR   (("sig_whi_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar whi_nL       (("sig_whi_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar whi_nR       (("sig_whi_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar whi_gamma_BW (("sig_whi_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar whi_sigma_EBE(("sig_whi_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, whi_sigma_CB));
        RooFormulaVar whi_sigma_EBM(("sig_whi_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar whi_norm     ("WH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_whi, whi_xsecbr, yield_var_whi));

        RooFormulaVar whiA_mean_CB  (("sig_whi_ALT_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch, do7TeV, doFFT).c_str() , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar whiA_sigma_CB (("sig_whi_ALT_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar whiA_alphaL   (("sig_whi_ALT_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar whiA_alphaR   (("sig_whi_ALT_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar whiA_nL       (("sig_whi_ALT_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar whiA_nR       (("sig_whi_ALT_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar whiA_gamma_BW (("sig_whi_ALT_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar whiA_sigma_EBE(("sig_whi_ALT_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, whi_sigma_CB));
        RooFormulaVar whiA_sigma_EBM(("sig_whi_ALT_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar whiA_norm     ("WH_ALT_norm"                            , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_whiA, whiA_xsecbr, yield_var_whiA));
 
        RooFormulaVar zhi_mean_CB  (("sig_zhi_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar zhi_sigma_CB (("sig_zhi_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar zhi_alphaL   (("sig_zhi_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar zhi_alphaR   (("sig_zhi_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar zhi_nL       (("sig_zhi_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar zhi_nR       (("sig_zhi_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar zhi_gamma_BW (("sig_zhi_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar zhi_sigma_EBE(("sig_zhi_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, zhi_sigma_CB));
        RooFormulaVar zhi_sigma_EBM(("sig_zhi_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar zhi_norm     ("ZH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_zhi, zhi_xsecbr, yield_var_zhi));

        RooFormulaVar zhiA_mean_CB  (("sig_zhi_ALT_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar zhiA_sigma_CB (("sig_zhi_ALT_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar zhiA_alphaL   (("sig_zhi_ALT_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar zhiA_alphaR   (("sig_zhi_ALT_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar zhiA_nL       (("sig_zhi_ALT_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar zhiA_nR       (("sig_zhi_ALT_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar zhiA_gamma_BW (("sig_zhi_ALT_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar zhiA_sigma_EBE(("sig_zhi_ALT_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, zhi_sigma_CB));
        RooFormulaVar zhiA_sigma_EBM(("sig_zhi_ALT_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar zhiA_norm    ("ZH_ALT_norm"                            , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_zhiA, zhiA_xsecbr, yield_var_zhiA));
 
        RooFormulaVar tth_mean_CB  (("sig_tth_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar tth_sigma_CB (("sig_tth_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar tth_alphaL   (("sig_tth_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar tth_alphaR   (("sig_tth_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar tth_nL       (("sig_tth_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar tth_nR       (("sig_tth_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar tth_gamma_BW (("sig_tth_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar tth_sigma_EBE(("sig_tth_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, tth_sigma_CB));
        RooFormulaVar tth_sigma_EBM(("sig_tth_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar tth_norm     ("ttH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_tth, tth_xsecbr, yield_var_tth));

        RooFormulaVar tthA_mean_CB  (("sig_tth_ALT_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar tthA_sigma_CB (("sig_tth_ALT_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar tthA_alphaL   (("sig_tth_ALT_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar tthA_alphaR   (("sig_tth_ALT_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar tthA_nL       (("sig_tth_ALT_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar tthA_nR       (("sig_tth_ALT_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar tthA_gamma_BW (("sig_tth_ALT_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar tthA_sigma_EBE(("sig_tth_ALT_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, tth_sigma_CB));
        RooFormulaVar tthA_sigma_EBM(("sig_tth_ALT_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar tthA_norm    ("ttH_ALT_norm"                            , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_tthA, tthA_xsecbr, yield_var_tthA));

        RooRealVar ggh_sigma_EBK(("sig_ggh_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);
        RooRealVar vbf_sigma_EBK(("sig_vbf_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);
        RooRealVar whi_sigma_EBK(("sig_whi_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);
        RooRealVar zhi_sigma_EBK(("sig_zhi_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);
        RooRealVar tth_sigma_EBK(("sig_tth_"+chstr+tevstr+"_sigma_EBK").c_str(), "EBE Kappa", 1.4);

        RooRealVar ggh_scale_BW (("sig_ggh_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar vbf_scale_BW (("sig_vbf_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar whi_scale_BW (("sig_whi_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar zhi_scale_BW (("sig_zhi_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar tth_scale_BW (("sig_tth_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);

        RooRealVar gghA_scale_BW (("sig_ggh_ALT_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar vbfA_scale_BW (("sig_vbf_ALT_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar whiA_scale_BW (("sig_whi_ALT_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar zhiA_scale_BW (("sig_zhi_ALT_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar tthA_scale_BW (("sig_tth_ALT_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);

        /////////// Set parameters to constant //////////////////
        qqzz_a0      .setConstant(kTRUE);
        qqzz_a1      .setConstant(kTRUE);
        qqzz_a2      .setConstant(kTRUE);
        qqzz_a3      .setConstant(kTRUE);
        qqzz_a4      .setConstant(kTRUE);
        qqzz_a5      .setConstant(kTRUE);
        qqzz_a6      .setConstant(kTRUE);
        qqzz_a7      .setConstant(kTRUE);
        qqzz_a8      .setConstant(kTRUE);
        qqzz_a9      .setConstant(kTRUE);
        qqzz_a10     .setConstant(kTRUE);
        qqzz_a11     .setConstant(kTRUE);
        qqzz_a12     .setConstant(kTRUE);
        qqzz_a13     .setConstant(kTRUE);
        
        ggzz_a0      .setConstant(kTRUE);
        ggzz_a1      .setConstant(kTRUE);
        ggzz_a2      .setConstant(kTRUE);
        ggzz_a3      .setConstant(kTRUE);
        ggzz_a4      .setConstant(kTRUE);
        ggzz_a5      .setConstant(kTRUE);
        ggzz_a6      .setConstant(kTRUE);
        ggzz_a7      .setConstant(kTRUE);
        ggzz_a8      .setConstant(kTRUE);
        ggzz_a9      .setConstant(kTRUE);
        
        zx_mean      .setConstant(kTRUE);
        zx_sigma     .setConstant(kTRUE);
      
        qqzz_sigma_EBK.setConstant(kTRUE);
        ggzz_sigma_EBK.setConstant(kTRUE);
        zjet_sigma_EBK.setConstant(kTRUE);
 
        masshiggs    .setConstant(kTRUE);
        sig_mean_err .setConstant(kTRUE);
        sig_sigma_err.setConstant(kTRUE);
        sig_gamma_err.setConstant(kTRUE);
        
        ggh_sigma_EBK.setConstant(kTRUE);
        vbf_sigma_EBK.setConstant(kTRUE);
        whi_sigma_EBK.setConstant(kTRUE);
        zhi_sigma_EBK.setConstant(kTRUE);
        tth_sigma_EBK.setConstant(kTRUE);
        
        ggh_scale_BW.setConstant(kTRUE);
        vbf_scale_BW.setConstant(kTRUE);
        whi_scale_BW.setConstant(kTRUE);
        zhi_scale_BW.setConstant(kTRUE);
        tth_scale_BW.setConstant(kTRUE);
        gghA_scale_BW.setConstant(kTRUE);
        vbfA_scale_BW.setConstant(kTRUE);
        whiA_scale_BW.setConstant(kTRUE);
        zhiA_scale_BW.setConstant(kTRUE);
        tthA_scale_BW.setConstant(kTRUE);
        
        ////////////////// Define the PDFs /////////////////////////////////
        
        const char* bkg_qqzz_pdf_name  = (do1D && !doMassError) ? "bkg_qqzz"  : "bkg_qqzz_1D" ;
        const char* bkg_ggzz_pdf_name  = (do1D && !doMassError) ? "bkg_ggzz"  : "bkg_ggzz_1D" ;
        const char* bkg_zjets_pdf_name = (do1D && !doMassError) ? "bkg_zjets" : "bkg_zjets_1D";
        const char* sig_ggH_pdf_name   = (do1D && !doMassError) ? "ggH"       : "ggH_1D"  ;
        const char* sig_VBF_pdf_name   = (do1D && !doMassError) ? "qqH"       : "qqH_1D"  ;
        const char* sig_WHi_pdf_name   = (do1D && !doMassError) ? "WH"        : "WH_1D"  ;
        const char* sig_ZHi_pdf_name   = (do1D && !doMassError) ? "ZH"        : "ZH_1D"  ;
        const char* sig_ttH_pdf_name   = (do1D && !doMassError) ? "ttH"       : "ttH_1D"  ;
       
        const char* sig_ggHA_pdf_name   = (do1D && !doMassError) ? "ggH_ALT"   : "ggH_1D_ALT"  ;
        const char* sig_VBFA_pdf_name   = (do1D && !doMassError) ? "qqH_ALT"   : "qqH_1D_ALT"  ;
        const char* sig_WHiA_pdf_name   = (do1D && !doMassError) ? "WH_ALT"    : "WH_1D_ALT"  ;
        const char* sig_ZHiA_pdf_name   = (do1D && !doMassError) ? "ZH_ALT"    : "ZH_1D_ALT"  ;
        const char* sig_ttHA_pdf_name   = (do1D && !doMassError) ? "ttH_ALT"   : "ttH_1D_ALT"  ;

        RooqqZZPdf_v2 bkg_qqzz_pdf(bkg_qqzz_pdf_name,"",CMS_zz4l_mass_1D,
                           qqzz_a0,
                           qqzz_a1,
                           qqzz_a2,
                           qqzz_a3,
                           qqzz_a4,
                           qqzz_a5,
                           qqzz_a6,
                           qqzz_a7,
                           qqzz_a8,
                           qqzz_a9,
                           qqzz_a10,
                           qqzz_a11,
                           qqzz_a12,
                           qqzz_a13);
        
        RooggZZPdf_v2 bkg_ggzz_pdf(bkg_ggzz_pdf_name,"",CMS_zz4l_mass_1D,
                           ggzz_a0,
                           ggzz_a1,
                           ggzz_a2,
                           ggzz_a3,
                           ggzz_a4,
                           ggzz_a5,
                           ggzz_a6,
                           ggzz_a7,
                           ggzz_a8,
                           ggzz_a9);
        
        RooLandau bkg_zjets_pdf(bkg_zjets_pdf_name, "",CMS_zz4l_mass_1D,zx_mean,zx_sigma);
      
        RooLognormal  bkgLN_qqzz_EBE("bkgLN_qqzz_EBE", "", CMS_zz4l_massErr, qqzz_sigma_EBM, qqzz_sigma_EBK);
        RooLognormal  bkgLN_ggzz_EBE("bkgLN_ggzz_EBE", "", CMS_zz4l_massErr, ggzz_sigma_EBM, ggzz_sigma_EBK);
        RooLognormal  bkgLN_zjet_EBE("bkgLN_zjet_EBE", "", CMS_zz4l_massErr, zjet_sigma_EBM, zjet_sigma_EBK);
        

        RooRelBWUFParam signalBW_ggH_LM("signalBW_ggH_LM", "", CMS_zz4l_mass_1D, masshiggs,ggh_scale_BW);
        RooRelBWUFParam signalBW_VBF_LM("signalBW_VBF_LM", "", CMS_zz4l_mass_1D, masshiggs,vbf_scale_BW);
        RooRelBWUFParam signalBW_WHi_LM("signalBW_WHi_LM", "", CMS_zz4l_mass_1D, masshiggs,whi_scale_BW);
        RooRelBWUFParam signalBW_ZHi_LM("signalBW_ZHi_LM", "", CMS_zz4l_mass_1D, masshiggs,zhi_scale_BW);
        RooRelBWUFParam signalBW_ttH_LM("signalBW_ttH_LM", "", CMS_zz4l_mass_1D, masshiggs,tth_scale_BW);

        RooRelBWHighMass  signalBW_ggH_HM("signalBW_ggH_HM", "", CMS_zz4l_mass_1D, masshiggs,ggh_gamma_BW);
        RooRelBWHighMass  signalBW_VBF_HM("signalBW_VBF_HM", "", CMS_zz4l_mass_1D, masshiggs,vbf_gamma_BW);
        RooRelBWHighMass  signalBW_WHi_HM("signalBW_WHi_HM", "", CMS_zz4l_mass_1D, masshiggs,whi_gamma_BW);
        RooRelBWHighMass  signalBW_ZHi_HM("signalBW_ZHi_HM", "", CMS_zz4l_mass_1D, masshiggs,zhi_gamma_BW);
        RooRelBWHighMass  signalBW_ttH_HM("signalBW_ttH_HM", "", CMS_zz4l_mass_1D, masshiggs,tth_gamma_BW);

        RooRelBWUFParam signalBW_ggHA_LM("signalBW_ggH_LM_ALT", "", CMS_zz4l_mass_1D, masshiggs,gghA_scale_BW);
        RooRelBWUFParam signalBW_VBFA_LM("signalBW_VBF_LM_ALT", "", CMS_zz4l_mass_1D, masshiggs,vbfA_scale_BW);
        RooRelBWUFParam signalBW_WHiA_LM("signalBW_WHi_LM_ALT", "", CMS_zz4l_mass_1D, masshiggs,whiA_scale_BW);
        RooRelBWUFParam signalBW_ZHiA_LM("signalBW_ZHi_LM_ALT", "", CMS_zz4l_mass_1D, masshiggs,zhiA_scale_BW);
        RooRelBWUFParam signalBW_ttHA_LM("signalBW_ttH_LM_ALT", "", CMS_zz4l_mass_1D, masshiggs,tthA_scale_BW);

        RooRelBWHighMass  signalBW_ggHA_HM("signalBW_ggH_HM_ALT", "", CMS_zz4l_mass_1D, masshiggs,ggh_gamma_BW);
        RooRelBWHighMass  signalBW_VBFA_HM("signalBW_VBF_HM_ALT", "", CMS_zz4l_mass_1D, masshiggs,vbf_gamma_BW);
        RooRelBWHighMass  signalBW_WHiA_HM("signalBW_WHi_HM_ALT", "", CMS_zz4l_mass_1D, masshiggs,whi_gamma_BW);
        RooRelBWHighMass  signalBW_ZHiA_HM("signalBW_ZHi_HM_ALT", "", CMS_zz4l_mass_1D, masshiggs,zhi_gamma_BW);
        RooRelBWHighMass  signalBW_ttHA_HM("signalBW_ttH_HM_ALT", "", CMS_zz4l_mass_1D, masshiggs,tth_gamma_BW);


        RooDoubleCB  signalCB_ggH   (doFFT ? "signalCB_ggH" : sig_ggH_pdf_name, "", CMS_zz4l_mass_1D, ggh_mean_CB,ggh_sigma_CB,ggh_alphaL,ggh_nL,ggh_alphaR,ggh_nR);
        RooDoubleCB  signalCB_VBF   (doFFT ? "signalCB_VBF" : sig_VBF_pdf_name, "", CMS_zz4l_mass_1D, vbf_mean_CB,vbf_sigma_CB,vbf_alphaL,vbf_nL,vbf_alphaR,vbf_nR);
        RooDoubleCB  signalCB_WHi   (doFFT ? "signalCB_WHi" : sig_WHi_pdf_name, "", CMS_zz4l_mass_1D, whi_mean_CB,whi_sigma_CB,whi_alphaL,whi_nL,whi_alphaR,whi_nR);
        RooDoubleCB  signalCB_ZHi   (doFFT ? "signalCB_ZHi" : sig_ZHi_pdf_name, "", CMS_zz4l_mass_1D, zhi_mean_CB,zhi_sigma_CB,zhi_alphaL,zhi_nL,zhi_alphaR,zhi_nR);
        RooDoubleCB  signalCB_ttH   (doFFT ? "signalCB_ttH" : sig_ttH_pdf_name, "", CMS_zz4l_mass_1D, tth_mean_CB,tth_sigma_CB,tth_alphaL,tth_nL,tth_alphaR,tth_nR);

        RooDoubleCB  signalCB_ggHA  (doFFT ? "signalCB_ggH_ALT" : sig_ggHA_pdf_name, "", CMS_zz4l_mass_1D, gghA_mean_CB,gghA_sigma_CB,gghA_alphaL,gghA_nL,gghA_alphaR,gghA_nR);
        RooDoubleCB  signalCB_VBFA  (doFFT ? "signalCB_VBF_ALT" : sig_VBFA_pdf_name, "", CMS_zz4l_mass_1D, vbfA_mean_CB,vbfA_sigma_CB,vbfA_alphaL,vbfA_nL,vbfA_alphaR,vbfA_nR);
        RooDoubleCB  signalCB_WHiA  (doFFT ? "signalCB_WHi_ALT" : sig_WHiA_pdf_name, "", CMS_zz4l_mass_1D, whiA_mean_CB,whiA_sigma_CB,whiA_alphaL,whiA_nL,whiA_alphaR,whiA_nR);
        RooDoubleCB  signalCB_ZHiA  (doFFT ? "signalCB_ZHi_ALT" : sig_ZHiA_pdf_name, "", CMS_zz4l_mass_1D, zhiA_mean_CB,zhiA_sigma_CB,zhiA_alphaL,zhiA_nL,zhiA_alphaR,zhiA_nR);
        RooDoubleCB  signalCB_ttHA  (doFFT ? "signalCB_ttH_ALT" : sig_ttHA_pdf_name, "", CMS_zz4l_mass_1D, tthA_mean_CB,tthA_sigma_CB,tthA_alphaL,tthA_nL,tthA_alphaR,tthA_nR);

        RooFFTConvPdf   sig_ggH_pdf_LM((sig_ggH_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ggH_LM, signalCB_ggH,2);
        RooFFTConvPdf   sig_VBF_pdf_LM((sig_VBF_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_VBF_LM, signalCB_VBF,2);
        RooFFTConvPdf   sig_WHi_pdf_LM((sig_WHi_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_WHi_LM, signalCB_WHi,2);
        RooFFTConvPdf   sig_ZHi_pdf_LM((sig_ZHi_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ZHi_LM, signalCB_ZHi,2);
        RooFFTConvPdf   sig_ttH_pdf_LM((sig_ttH_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ttH_LM, signalCB_ttH,2);

        RooFFTConvPdf   sig_ggH_pdf_HM((sig_ggH_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ggH_HM, signalCB_ggH,2);
        RooFFTConvPdf   sig_VBF_pdf_HM((sig_VBF_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_VBF_HM, signalCB_VBF,2);
        RooFFTConvPdf   sig_WHi_pdf_HM((sig_WHi_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_WHi_HM, signalCB_WHi,2);
        RooFFTConvPdf   sig_ZHi_pdf_HM((sig_ZHi_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ZHi_HM, signalCB_ZHi,2);
        RooFFTConvPdf   sig_ttH_pdf_HM((sig_ttH_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ttH_HM, signalCB_ttH,2);

        RooFFTConvPdf&  sig_ggH_pdf = mass<400 ? sig_ggH_pdf_LM : sig_ggH_pdf_HM;
        RooFFTConvPdf&  sig_VBF_pdf = mass<400 ? sig_VBF_pdf_LM : sig_VBF_pdf_HM;
        RooFFTConvPdf&  sig_WHi_pdf = mass<400 ? sig_WHi_pdf_LM : sig_WHi_pdf_HM;
        RooFFTConvPdf&  sig_ZHi_pdf = mass<400 ? sig_ZHi_pdf_LM : sig_ZHi_pdf_HM;
        RooFFTConvPdf&  sig_ttH_pdf = mass<400 ? sig_ttH_pdf_LM : sig_ttH_pdf_HM;

        if (doFFT) {
            sig_ggH_pdf.SetName(sig_ggH_pdf_name);
            sig_VBF_pdf.SetName(sig_VBF_pdf_name);
            sig_WHi_pdf.SetName(sig_WHi_pdf_name);
            sig_ZHi_pdf.SetName(sig_ZHi_pdf_name);
            sig_ttH_pdf.SetName(sig_ttH_pdf_name);
        }

        sig_ggH_pdf.setBufferFraction(0.2);
        sig_VBF_pdf.setBufferFraction(0.2);
        sig_WHi_pdf.setBufferFraction(0.2);
        sig_ZHi_pdf.setBufferFraction(0.2);
        sig_ttH_pdf.setBufferFraction(0.2);

        RooFFTConvPdf   sig_ggHA_pdf_LM((sig_ggH_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ggHA_LM, signalCB_ggHA,2);
        RooFFTConvPdf   sig_VBFA_pdf_LM((sig_VBF_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_VBFA_LM, signalCB_VBFA,2);
        RooFFTConvPdf   sig_WHiA_pdf_LM((sig_WHi_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_WHiA_LM, signalCB_WHiA,2);
        RooFFTConvPdf   sig_ZHiA_pdf_LM((sig_ZHi_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ZHiA_LM, signalCB_ZHiA,2);
        RooFFTConvPdf   sig_ttHA_pdf_LM((sig_ttH_pdf_name + std::string("_LM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ttHA_LM, signalCB_ttHA,2);

        RooFFTConvPdf   sig_ggHA_pdf_HM((sig_ggH_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ggHA_HM, signalCB_ggHA,2);
        RooFFTConvPdf   sig_VBFA_pdf_HM((sig_VBF_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_VBFA_HM, signalCB_VBFA,2);
        RooFFTConvPdf   sig_WHiA_pdf_HM((sig_WHi_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_WHiA_HM, signalCB_WHiA,2);
        RooFFTConvPdf   sig_ZHiA_pdf_HM((sig_ZHi_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ZHiA_HM, signalCB_ZHiA,2);
        RooFFTConvPdf   sig_ttHA_pdf_HM((sig_ttH_pdf_name + std::string("_HM")).c_str(), "", CMS_zz4l_mass_1D, signalBW_ttHA_HM, signalCB_ttHA,2);

        RooFFTConvPdf&  sig_ggHA_pdf = mass<400 ? sig_ggHA_pdf_LM : sig_ggHA_pdf_HM;
        RooFFTConvPdf&  sig_VBFA_pdf = mass<400 ? sig_VBFA_pdf_LM : sig_VBFA_pdf_HM;
        RooFFTConvPdf&  sig_WHiA_pdf = mass<400 ? sig_WHiA_pdf_LM : sig_WHiA_pdf_HM;
        RooFFTConvPdf&  sig_ZHiA_pdf = mass<400 ? sig_ZHiA_pdf_LM : sig_ZHiA_pdf_HM;
        RooFFTConvPdf&  sig_ttHA_pdf = mass<400 ? sig_ttHA_pdf_LM : sig_ttHA_pdf_HM;

        if (doFFT) {
            sig_ggHA_pdf.SetName(sig_ggHA_pdf_name);
            sig_VBFA_pdf.SetName(sig_VBFA_pdf_name);
            sig_WHiA_pdf.SetName(sig_WHiA_pdf_name);
            sig_ZHiA_pdf.SetName(sig_ZHiA_pdf_name);
            sig_ttHA_pdf.SetName(sig_ttHA_pdf_name);
        }

        sig_ggHA_pdf.setBufferFraction(0.2);
        sig_VBFA_pdf.setBufferFraction(0.2);
        sig_WHiA_pdf.setBufferFraction(0.2);
        sig_ZHiA_pdf.setBufferFraction(0.2);
        sig_ttHA_pdf.setBufferFraction(0.2);


        RooDoubleCB  signalCB_ggH_EBE("signalCB_ggH_EBE", "", CMS_zz4l_mass_1D, ggh_mean_CB,ggh_sigma_EBE,ggh_alphaL,ggh_alphaR,ggh_nL,ggh_nR);
        RooDoubleCB  signalCB_VBF_EBE("signalCB_VBF_EBE", "", CMS_zz4l_mass_1D, vbf_mean_CB,vbf_sigma_EBE,vbf_alphaL,vbf_alphaR,vbf_nL,vbf_nR);
        RooDoubleCB  signalCB_WHi_EBE("signalCB_WHi_EBE", "", CMS_zz4l_mass_1D, whi_mean_CB,whi_sigma_EBE,whi_alphaL,whi_alphaR,whi_nL,whi_nR);
        RooDoubleCB  signalCB_ZHi_EBE("signalCB_ZHi_EBE", "", CMS_zz4l_mass_1D, zhi_mean_CB,zhi_sigma_EBE,zhi_alphaL,zhi_alphaR,zhi_nL,zhi_nR);
        RooDoubleCB  signalCB_ttH_EBE("signalCB_ttH_EBE", "", CMS_zz4l_mass_1D, tth_mean_CB,tth_sigma_EBE,tth_alphaL,tth_alphaR,tth_nL,tth_nR);

        RooFFTConvPdf   sig_ggH_pdf_EBE("signalFFT_ggH_EBE", "", CMS_zz4l_mass_1D, signalBW_ggH_LM, signalCB_ggH_EBE,2);
        RooFFTConvPdf   sig_VBF_pdf_EBE("signalFFT_VBF_EBE", "", CMS_zz4l_mass_1D, signalBW_VBF_LM, signalCB_VBF_EBE,2);
        RooFFTConvPdf   sig_WHi_pdf_EBE("signalFFT_WHi_EBE", "", CMS_zz4l_mass_1D, signalBW_WHi_LM, signalCB_WHi_EBE,2);
        RooFFTConvPdf   sig_ZHi_pdf_EBE("signalFFT_ZHi_EBE", "", CMS_zz4l_mass_1D, signalBW_ZHi_LM, signalCB_ZHi_EBE,2);
        RooFFTConvPdf   sig_ttH_pdf_EBE("signalFFT_ttH_EBE", "", CMS_zz4l_mass_1D, signalBW_ttH_LM, signalCB_ttH_EBE,2);

        sig_ggH_pdf_EBE.setBufferFraction(0.2);
        sig_VBF_pdf_EBE.setBufferFraction(0.2);
        sig_WHi_pdf_EBE.setBufferFraction(0.2);
        sig_ZHi_pdf_EBE.setBufferFraction(0.2);
        sig_ttH_pdf_EBE.setBufferFraction(0.2);

        RooLognormal  signalLN_ggH_EBE("signalLN_ggH_EBE", "", CMS_zz4l_massErr, ggh_sigma_EBM, ggh_sigma_EBK);
        RooLognormal  signalLN_VBF_EBE("signalLN_VBF_EBE", "", CMS_zz4l_massErr, vbf_sigma_EBM, vbf_sigma_EBK);
        RooLognormal  signalLN_WHi_EBE("signalLN_WHi_EBE", "", CMS_zz4l_massErr, whi_sigma_EBM, whi_sigma_EBK);
        RooLognormal  signalLN_ZHi_EBE("signalLN_ZHi_EBE", "", CMS_zz4l_massErr, zhi_sigma_EBM, zhi_sigma_EBK);
        RooLognormal  signalLN_ttH_EBE("signalLN_ttH_EBE", "", CMS_zz4l_massErr, tth_sigma_EBM, tth_sigma_EBK);
        

        //////////////// Yields and printing cards /////////////////////////

        RooRealVar CMS_zz4l_mass_norm("CMS_zz4l_mass_norm", "M(4l)", massLowBkgFit, massHighBkgFit, "GeV/c^{2}");
        CMS_zz4l_mass_norm.setRange("card", massLow      , massHigh);
        CMS_zz4l_mass_norm.setRange("full", massLowBkgFit, massHighBkgFit);

        RooqqZZPdf_v2 bkg_qqzz_norm((std::string(bkg_qqzz_pdf_name)+"_normalization").c_str(),"",CMS_zz4l_mass_norm,
                           qqzz_a0,
                           qqzz_a1,
                           qqzz_a2,
                           qqzz_a3,
                           qqzz_a4,
                           qqzz_a5,
                           qqzz_a6,
                           qqzz_a7,
                           qqzz_a8,
                           qqzz_a9,
                           qqzz_a10,
                           qqzz_a11,
                           qqzz_a12,
                           qqzz_a13);

        RooggZZPdf_v2 bkg_ggzz_norm((std::string(bkg_ggzz_pdf_name)+"_normalization").c_str(),"",CMS_zz4l_mass_norm,
                           ggzz_a0,
                           ggzz_a1,
                           ggzz_a2,
                           ggzz_a3,
                           ggzz_a4,
                           ggzz_a5,
                           ggzz_a6,
                           ggzz_a7,
                           ggzz_a8,
                           ggzz_a9);

        RooLandau bkg_zjet_norm((std::string(bkg_zjets_pdf_name)+"_normalization").c_str(), "",CMS_zz4l_mass_norm,zx_mean,zx_sigma);
        
        double qqzz_fullyield = 0.0;
        double ggzz_fullyield = 0.0;
        double zjet_fullyield = 0.0;

        if (ch == 0) {
            qqzz_fullyield = do7TeV ? 21.0 : 56.5;
            ggzz_fullyield = do7TeV ? 1.30 : 3.08;
            zjet_fullyield = do7TeV ? 0.70 : 3.7;
        }

        else if (ch == 1) {
            qqzz_fullyield = do7TeV ? 13.4 : 37.0;
            ggzz_fullyield = do7TeV ? 0.92 : 2.24;
            zjet_fullyield = do7TeV ? 2.40 : 9.5;
        }

        else {
            qqzz_fullyield = do7TeV ? 32.2 : 90.3;
            ggzz_fullyield = do7TeV ? 2.18 : 6.55;
            zjet_fullyield = do7TeV ? 2.90 : 12.7;
        }

         //**************************************
         //Apply lumi scaling for easy projections
         //**************************************
         qqzz_fullyield = qqzz_fullyield * ( lumi / base_lumi);
         ggzz_fullyield = ggzz_fullyield * ( lumi / base_lumi);
         zjet_fullyield = zjet_fullyield * ( lumi / base_lumi);

        float yield_qq = qqzz_fullyield * ((bkg_qqzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("card")))->getVal()) / ((bkg_qqzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("full")))->getVal());
        float yield_gg = ggzz_fullyield * ((bkg_ggzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("card")))->getVal()) / ((bkg_ggzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("full")))->getVal());
        float yield_zj = zjet_fullyield * ((bkg_zjet_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("card")))->getVal()) / ((bkg_zjet_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("full")))->getVal());


        std::string card   = "";
        if (mass<=300) card += createCardTemplateForSignalHypothesisSingleChannel(do7TeV, ch, do1D, workspace.c_str(), mass<400. ? true : false);
        else card += createCardTemplateForSignalHypothesisSingleChannel(do7TeV, ch, do1D, workspace.c_str(), mass<400. ? true : false);

        std::string binname;
        if (ch == 0) binname = "a1";
        if (ch == 1) binname = "a2";
        if (ch == 2) binname = "a3";

        card = findAndReplace(card, "GGZZ_PDF"       , getGGZZPDFUncertainty7TeV(mass));
        card = findAndReplace(card, "QQZZ_PDF"       , getQQZZPDFUncertainty7TeV(mass));
        card = findAndReplace(card, "GGZZ_QCD"       , getGGZZQCDScaleUncertainty7TeV(mass));
        card = findAndReplace(card, "QQZZ_QCD"       , getQQZZQCDScaleUncertainty7TeV(mass));
        card = findAndReplace(card, "GGH_PDF"        , getggHPDFUncertainty(mass, false), getggHPDFUncertainty(mass, true));
        card = findAndReplace(card, "VBF_PDF"        , getVBFPDFUncertainty(mass, false), getVBFPDFUncertainty(mass, true));
        card = findAndReplace(card, "WHI_PDF"        , getWHiPDFUncertainty(mass, false), getWHiPDFUncertainty(mass, true));
        card = findAndReplace(card, "ZHI_PDF"        , getZHiPDFUncertainty(mass, false), getZHiPDFUncertainty(mass, true));
        card = findAndReplace(card, "TTH_PDF"        , getttHPDFUncertainty(mass, false), getttHPDFUncertainty(mass, true));
        card = findAndReplace(card, "GGH_QCD"        , getggHQCDScaleUncertainty(mass, false), getggHQCDScaleUncertainty(mass, true));
        card = findAndReplace(card, "VBF_QCD"        , getVBFQCDScaleUncertainty(mass, false), getVBFQCDScaleUncertainty(mass, true));
        card = findAndReplace(card, "WHI_QCD"        , getWHiQCDScaleUncertainty(mass, false), getWHiQCDScaleUncertainty(mass, true));
        card = findAndReplace(card, "ZHI_QCD"        , getZHiQCDScaleUncertainty(mass, false), getZHiQCDScaleUncertainty(mass, true));
        card = findAndReplace(card, "TTH_QCD"        , getttHQCDScaleUncertainty(mass, false), getttHQCDScaleUncertainty(mass, true));
        card = findAndReplace(card, "HIGH_MH"        , mass < 200 ? 1 : 1+(mass/1000.));
        card = findAndReplace(card, "SIG_GGH_YIELD"  , (getXsecggHByChannel(mass, ch) + getXsecVBFByChannel(mass, ch) + getXsecZHiByChannel(mass, ch) + getXsecWHiByChannel(mass, ch) + getXsecttHByChannel(mass, ch))/getXsecggHByChannel(mass, ch));
        card = findAndReplace(card, "BKG_QQZZ_YIELD" , yield_qq);
        card = findAndReplace(card, "BKG_GGZZ_YIELD" , yield_gg);
        card = findAndReplace(card, "BKG_ZJETS_YIELD", yield_zj);
        card = findAndReplace(card, "BIN"            , binname);
        card = findAndReplace(card, "OBS"            , yield_dt);
        card = findAndReplace(card, "LOW_EDGE"       , massLow);
        card = findAndReplace(card, "HIGH_EDGE"      , massHigh);

        ofstream file;
        file.open ((card_name +".txt").c_str());
        file << card;
        file.close();



        w.import(ggh_norm);
        w.import(vbf_norm);
        w.import(whi_norm);
        w.import(zhi_norm);
        w.import(tth_norm);
        w.import(gghA_norm);
        w.import(vbfA_norm);
        w.import(whiA_norm);
        w.import(zhiA_norm);
        w.import(tthA_norm);

        //*******************************************************************************************
        //Do 3D Analysis
        //*******************************************************************************************
        cout << "here71\n";
        TFile MELAFile(melafilename.c_str());
        RooArgList v2dList_MassMELA(CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
        RooArgSet  v2dSet_MassMELA (CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
        TFile alternativeMELAFile(alternativeMELAfilename.c_str());
        RooArgList v2dList_MassAlternativeMELA(CMS_zz4l_mass_1D, CMS_zz4l_AlternativeMELA);
        RooArgSet  v2dSet_MassAlternativeMELA (CMS_zz4l_mass_1D, CMS_zz4l_AlternativeMELA);
            
        std::string chstrsmall;
        if (ch == 0) chstrsmall = "mm";
        if (ch == 1) chstrsmall = "ee";
        if (ch == 2) chstrsmall = "em";


        //******************************************************************************************
        //Mass Vs MELA Template PDFs
        //******************************************************************************************
        TH2F* MELAShape_qz = (TH2F*)(MELAFile.Get(("hist2D_bkg_"+chstrsmall).c_str()));
        TH2F* MELAShape_gz = (TH2F*)(MELAFile.Get(("hist2D_b2g_"+chstrsmall).c_str()));
        TH2F* MELAShape_si = (TH2F*)(MELAFile.Get(("hist2D_sig_"+chstrsmall).c_str()));

        TH2F* MELAShape_s2 = 0;
        if (alternativeSignalType == 1) MELAShape_s2 = (TH2F*)(MELAFile.Get(("hist2D_pssig_"+chstrsmall).c_str()));
        else if (alternativeSignalType == 2) MELAShape_s2 = (TH2F*)(MELAFile.Get(("hist2D_spin2minimalsig_"+chstrsmall).c_str()));
        else { cout << "Error. signal type = " <<  alternativeSignalType << " not supported\n"; return;}

        TH2F* MELAShape_zu = (TH2F*)(MELAFile.Get(("hist2D_bkg_"+chstrsmall+"_up").c_str()));
        TH2F* MELAShape_zd = (TH2F*)(MELAFile.Get(("hist2D_bkg_"+chstrsmall+"_dn").c_str()));

        assert(MELAShape_qz);
        assert(MELAShape_gz);
        assert(MELAShape_si);
        assert(MELAShape_s2);
        assert(MELAShape_zu);
        assert(MELAShape_zd);

        RooDataHist rhist_MELA_qqzz (("rhist_MELA_qqzz_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_qz);
        RooDataHist rhist_MELA_ggzz (("rhist_MELA_ggzz_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_gz);
        RooDataHist rhist_MELA_zjets(("rhist_MELA_zjets_"+chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_qz);
        RooDataHist rhist_MELA_ggH  (("rhist_MELA_ggH_"  +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_si);
        RooDataHist rhist_MELA_VBF  (("rhist_MELA_VBF_"  +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_si);
        RooDataHist rhist_MELA_WHi  (("rhist_MELA_WHi_"  +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_si);
        RooDataHist rhist_MELA_ZHi  (("rhist_MELA_ZHi_"  +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_si);
        RooDataHist rhist_MELA_ttH  (("rhist_MELA_ttH_"  +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_si);
        RooDataHist rhist_MELA_ggHA (("rhist_MELA_ggHA_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_s2);
        RooDataHist rhist_MELA_VBFA (("rhist_MELA_VBFA_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_s2);
        RooDataHist rhist_MELA_WHiA (("rhist_MELA_WHiA_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_s2);
        RooDataHist rhist_MELA_ZHiA (("rhist_MELA_ZHiA_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_s2);
        RooDataHist rhist_MELA_ttHA (("rhist_MELA_ttHA_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_s2);
        RooDataHist rhist_MELA_zjup (("rhist_MELA_zjup_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_zu);
        RooDataHist rhist_MELA_zjdn (("rhist_MELA_zjdn_" +chstr+tevstr).c_str(), "", v2dList_MassMELA, MELAShape_zd);
        RooHistPdf rpdf_MELA_qqzz (("bkg_qqzz_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_qqzz);
        RooHistPdf rpdf_MELA_ggzz (("bkg_ggzz_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_ggzz);
        RooHistPdf rpdf_MELA_zjets(("bkg_zjets_MELA2D_pdf_"+chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_zjets);
        RooHistPdf rpdf_MELA_ggH  (("sig_ggH_MELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_ggH);
        RooHistPdf rpdf_MELA_VBF  (("sig_VBF_MELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_VBF);
        RooHistPdf rpdf_MELA_WHi  (("sig_WHi_MELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_WHi);
        RooHistPdf rpdf_MELA_ZHi  (("sig_ZHi_MELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_ZHi);
        RooHistPdf rpdf_MELA_ttH  (("sig_ttH_MELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_ttH);
        RooHistPdf rpdf_MELA_ggHA (("sig_ggHA_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_ggHA);
        RooHistPdf rpdf_MELA_VBFA (("sig_VBFA_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_VBFA);
        RooHistPdf rpdf_MELA_WHiA (("sig_WHiA_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_WHiA);
        RooHistPdf rpdf_MELA_ZHiA (("sig_ZHiA_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_ZHiA);
        RooHistPdf rpdf_MELA_ttHA (("sig_ttHA_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_ttHA);

        RooHistPdf rpdf_MELA_zjup (("bkg_zjup_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_zjup);
        RooHistPdf rpdf_MELA_zjdn (("bkg_zjdn_MELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassMELA , rhist_MELA_zjdn);

        RooArgList list_MELA_zjets;
        list_MELA_zjets.add(rpdf_MELA_zjets);
        list_MELA_zjets.add(rpdf_MELA_zjup);
        list_MELA_zjets.add(rpdf_MELA_zjdn);
        RooRealVar CMS_zz4l_bkgMELA("CMS_zz4l_bkgMELA" ,"" ,0,-10,10); 
 
        FastVerticalInterpHistPdf2D plpdf_MELA_qqzz (("bkg_qqzz_FVIHP_MELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_qqzz) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_ggzz (("bkg_ggzz_FVIHP_MELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_ggzz) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_zjets(("bkg_zjets_FVIHP_MELA_"+chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,list_MELA_zjets            ,RooArgList(CMS_zz4l_bkgMELA),1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_ggH  (("sig_ggH_FVIHP_MELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_ggH)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_VBF  (("sig_VBF_FVIHP_MELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_VBF)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_WHi  (("sig_WHi_FVIHP_MELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_WHi)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_ZHi  (("sig_ZHi_FVIHP_MELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_ZHi)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_ttH  (("sig_ttH_FVIHP_MELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_ttH)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_ggHA (("sig_ggHA_FVIHP_MELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_ggHA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_VBFA (("sig_VBFA_FVIHP_MELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_VBFA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_WHiA (("sig_WHiA_FVIHP_MELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_WHiA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_ZHiA (("sig_ZHiA_FVIHP_MELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_ZHiA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_MELA_ttHA (("sig_ttHA_FVIHP_MELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_MELA_ttHA) ,RooArgList()                ,1.0,1);

        //******************************************************************************************
        //Mass Vs Alternative MELA Template PDFs
        //******************************************************************************************
        TH2F* alternativeMELAShape_qz = (TH2F*)(alternativeMELAFile.Get(("hist2D_bkg_"+chstrsmall).c_str()));
        TH2F* alternativeMELAShape_gz = (TH2F*)(alternativeMELAFile.Get(("hist2D_b2g_"+chstrsmall).c_str()));
        TH2F* alternativeMELAShape_si = (TH2F*)(alternativeMELAFile.Get(("hist2D_sig_"+chstrsmall).c_str()));

        TH2F* alternativeMELAShape_s2 = 0;
        if (alternativeSignalType == 1) alternativeMELAShape_s2 = (TH2F*)(alternativeMELAFile.Get(("hist2D_pssig_"+chstrsmall).c_str()));
        else if (alternativeSignalType == 2) alternativeMELAShape_s2 = (TH2F*)(alternativeMELAFile.Get(("hist2D_spin2minimalsig_"+chstrsmall).c_str()));
        else { cout << "Error. signal type = " <<  alternativeSignalType << " not supported\n"; return;}

        TH2F* alternativeMELAShape_zu = (TH2F*)(alternativeMELAFile.Get(("hist2D_bkg_"+chstrsmall+"_up").c_str()));
        TH2F* alternativeMELAShape_zd = (TH2F*)(alternativeMELAFile.Get(("hist2D_bkg_"+chstrsmall+"_dn").c_str()));

        assert(alternativeMELAShape_qz);
        assert(alternativeMELAShape_gz);
        assert(alternativeMELAShape_si);
        assert(alternativeMELAShape_s2);
        assert(alternativeMELAShape_zu);
        assert(alternativeMELAShape_zd);

        RooDataHist rhist_alternativeMELA_qqzz (("rhist_alternativeMELA_qqzz_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_qz);
        RooDataHist rhist_alternativeMELA_ggzz (("rhist_alternativeMELA_ggzz_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_gz);
        RooDataHist rhist_alternativeMELA_zjets(("rhist_alternativeMELA_zjets_"+chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_qz);
        RooDataHist rhist_alternativeMELA_ggH  (("rhist_alternativeMELA_ggH_"  +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_si);
        RooDataHist rhist_alternativeMELA_VBF  (("rhist_alternativeMELA_VBF_"  +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_si);
        RooDataHist rhist_alternativeMELA_WHi  (("rhist_alternativeMELA_WHi_"  +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_si);
        RooDataHist rhist_alternativeMELA_ZHi  (("rhist_alternativeMELA_ZHi_"  +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_si);
        RooDataHist rhist_alternativeMELA_ttH  (("rhist_alternativeMELA_ttH_"  +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_si);
        RooDataHist rhist_alternativeMELA_ggHA (("rhist_alternativeMELA_ggHA_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_s2);
        RooDataHist rhist_alternativeMELA_VBFA (("rhist_alternativeMELA_VBFA_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_s2);
        RooDataHist rhist_alternativeMELA_WHiA (("rhist_alternativeMELA_WHiA_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_s2);
        RooDataHist rhist_alternativeMELA_ZHiA (("rhist_alternativeMELA_ZHiA_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_s2);
        RooDataHist rhist_alternativeMELA_ttHA (("rhist_alternativeMELA_ttHA_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_s2);
        RooDataHist rhist_alternativeMELA_zjup (("rhist_alternativeMELA_zjup_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_zu);
        RooDataHist rhist_alternativeMELA_zjdn (("rhist_alternativeMELA_zjdn_" +chstr+tevstr).c_str(), "", v2dList_MassAlternativeMELA, alternativeMELAShape_zd);
        RooHistPdf rpdf_alternativeMELA_qqzz (("bkg_qqzz_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_qqzz);
        RooHistPdf rpdf_alternativeMELA_ggzz (("bkg_ggzz_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_ggzz);
        RooHistPdf rpdf_alternativeMELA_zjets(("bkg_zjets_alternativeMELA2D_pdf_"+chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_zjets);
        RooHistPdf rpdf_alternativeMELA_ggH  (("sig_ggH_alternativeMELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_ggH);
        RooHistPdf rpdf_alternativeMELA_VBF  (("sig_VBF_alternativeMELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_VBF);
        RooHistPdf rpdf_alternativeMELA_WHi  (("sig_WHi_alternativeMELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_WHi);
        RooHistPdf rpdf_alternativeMELA_ZHi  (("sig_ZHi_alternativeMELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_ZHi);
        RooHistPdf rpdf_alternativeMELA_ttH  (("sig_ttH_alternativeMELA2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_ttH);
        RooHistPdf rpdf_alternativeMELA_ggHA (("sig_ggHA_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_ggHA);
        RooHistPdf rpdf_alternativeMELA_VBFA (("sig_VBFA_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_VBFA);
        RooHistPdf rpdf_alternativeMELA_WHiA (("sig_WHiA_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_WHiA);
        RooHistPdf rpdf_alternativeMELA_ZHiA (("sig_ZHiA_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_ZHiA);
        RooHistPdf rpdf_alternativeMELA_ttHA (("sig_ttHA_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_ttHA);

        RooHistPdf rpdf_alternativeMELA_zjup (("bkg_zjup_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_zjup);
        RooHistPdf rpdf_alternativeMELA_zjdn (("bkg_zjdn_alternativeMELA2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet_MassAlternativeMELA , rhist_alternativeMELA_zjdn);

        RooArgList list_alternativeMELA_zjets;
        list_alternativeMELA_zjets.add(rpdf_alternativeMELA_zjets);
        list_alternativeMELA_zjets.add(rpdf_alternativeMELA_zjup);
        list_alternativeMELA_zjets.add(rpdf_alternativeMELA_zjdn);
        cout << "here8\n";

        RooRealVar CMS_zz4l_bkgAlternativeMELA("CMS_zz4l_bkgAlternativeMELA" ,"" ,0,-10,10); 
 
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_qqzz (("bkg_qqzz_FVIHP_alternativeMELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_qqzz) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_ggzz (("bkg_ggzz_FVIHP_alternativeMELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_ggzz) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_zjets(("bkg_zjets_FVIHP_alternativeMELA_"+chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,list_alternativeMELA_zjets            ,RooArgList(CMS_zz4l_bkgAlternativeMELA),1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_ggH  (("sig_ggH_FVIHP_alternativeMELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_ggH)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_VBF  (("sig_VBF_FVIHP_alternativeMELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_VBF)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_WHi  (("sig_WHi_FVIHP_alternativeMELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_WHi)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_ZHi  (("sig_ZHi_FVIHP_alternativeMELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_ZHi)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_ttH  (("sig_ttH_FVIHP_alternativeMELA_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_ttH)  ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_ggHA (("sig_ggHA_FVIHP_alternativeMELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_ggHA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_VBFA (("sig_VBFA_FVIHP_alternativeMELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_VBFA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_WHiA (("sig_WHiA_FVIHP_alternativeMELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_WHiA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_ZHiA (("sig_ZHiA_FVIHP_alternativeMELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_ZHiA) ,RooArgList()                ,1.0,1);
        FastVerticalInterpHistPdf2D plpdf_alternativeMELA_ttHA (("sig_ttHA_FVIHP_alternativeMELA_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_AlternativeMELA,true,RooArgList(rpdf_alternativeMELA_ttHA) ,RooArgList()                ,1.0,1);


        //******************************************************************************************
        //Product PDFs : p_mass * p_2dtemplate
        //******************************************************************************************
        if (doFFT) { 
          RooProdPdf sig_ggH_pdf_2D  ("ggH"    ,  "", sig_ggH_pdf   ,Conditional(plpdf_alternativeMELA_ggH  , RooArgSet(CMS_zz4l_AlternativeMELA)), Conditional(plpdf_MELA_ggH  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBF_pdf_2D  ("qqH"    ,  "", sig_VBF_pdf   ,Conditional(plpdf_alternativeMELA_VBF  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_VBF  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_WHi_pdf_2D  ("WH"     ,  "", sig_WHi_pdf   ,Conditional(plpdf_alternativeMELA_WHi  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_WHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ZHi_pdf_2D  ("ZH"     ,  "", sig_ZHi_pdf   ,Conditional(plpdf_alternativeMELA_ZHi  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ZHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ttH_pdf_2D  ("ttH"    ,  "", sig_ttH_pdf   ,Conditional(plpdf_alternativeMELA_ttH  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ttH  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ggHA_pdf_2D ("ggH_ALT",  "", sig_ggHA_pdf  ,Conditional(plpdf_alternativeMELA_ggHA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ggHA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBFA_pdf_2D ("qqH_ALT",  "", sig_VBFA_pdf  ,Conditional(plpdf_alternativeMELA_VBFA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_VBFA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_WHiA_pdf_2D ("WH_ALT" ,  "", sig_WHiA_pdf  ,Conditional(plpdf_alternativeMELA_WHiA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_WHiA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ZHiA_pdf_2D ("ZH_ALT" ,  "", sig_ZHiA_pdf  ,Conditional(plpdf_alternativeMELA_ZHiA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ZHiA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ttHA_pdf_2D ("ttH_ALT",  "", sig_ttHA_pdf  ,Conditional(plpdf_alternativeMELA_ttHA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ttHA , RooArgSet(CMS_zz4l_melaLD))); 
                
            RooProdPdf bkg_qqzz_pdf_2D ("bkg_qqzz", "", bkg_qqzz_pdf  ,Conditional(plpdf_alternativeMELA_qqzz , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_qqzz , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_ggzz_pdf_2D ("bkg_ggzz", "", bkg_ggzz_pdf  ,Conditional(plpdf_alternativeMELA_ggzz , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ggzz , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_zjets_pdf_2D("bkg_zjets","", bkg_zjets_pdf ,Conditional(plpdf_alternativeMELA_zjets, RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_zjets, RooArgSet(CMS_zz4l_melaLD))); 
           
            w.import(bkg_qqzz_pdf_2D); 
            w.import(bkg_ggzz_pdf_2D); 
            w.import(bkg_zjets_pdf_2D); 
            w.import(sig_ggH_pdf_2D); 
            w.import(sig_VBF_pdf_2D);
            w.import(sig_WHi_pdf_2D); 
            w.import(sig_ZHi_pdf_2D);
            w.import(sig_ttH_pdf_2D);
            w.import(sig_ggHA_pdf_2D); 
            w.import(sig_VBFA_pdf_2D);
            w.import(sig_WHiA_pdf_2D); 
            w.import(sig_ZHiA_pdf_2D);
            w.import(sig_ttHA_pdf_2D);
        }
        else {
            RooProdPdf sig_ggH_pdf_2D  ("ggH",      "", signalCB_ggH  ,Conditional(plpdf_alternativeMELA_ggH  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ggH  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBF_pdf_2D  ("qqH",      "", signalCB_VBF  ,Conditional(plpdf_alternativeMELA_VBF  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_VBF  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_WHi_pdf_2D  ("WH" ,      "", signalCB_WHi  ,Conditional(plpdf_alternativeMELA_WHi  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_WHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ZHi_pdf_2D  ("ZH" ,      "", signalCB_ZHi  ,Conditional(plpdf_alternativeMELA_ZHi  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ZHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ttH_pdf_2D  ("ttH",      "", signalCB_ttH  ,Conditional(plpdf_alternativeMELA_ttH  , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ttH  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ggHA_pdf_2D ("ggH_ALT",  "", signalCB_ggHA ,Conditional(plpdf_alternativeMELA_ggHA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ggHA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBFA_pdf_2D ("qqH_ALT",  "", signalCB_VBFA ,Conditional(plpdf_alternativeMELA_VBFA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_VBFA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_WHiA_pdf_2D ("WH_ALT" ,  "", signalCB_WHiA ,Conditional(plpdf_alternativeMELA_WHiA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_WHiA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ZHiA_pdf_2D ("ZH_ALT" ,  "", signalCB_ZHiA ,Conditional(plpdf_alternativeMELA_ZHiA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ZHiA , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ttHA_pdf_2D ("ttH_ALT",  "", signalCB_ttHA ,Conditional(plpdf_alternativeMELA_ttHA , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ttHA , RooArgSet(CMS_zz4l_melaLD))); 

            RooProdPdf bkg_qqzz_pdf_2D ("bkg_qqzz", "", bkg_qqzz_pdf  ,Conditional(plpdf_alternativeMELA_qqzz , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_qqzz , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_ggzz_pdf_2D ("bkg_ggzz", "", bkg_ggzz_pdf  ,Conditional(plpdf_alternativeMELA_ggzz , RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_ggzz , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_zjets_pdf_2D("bkg_zjets","", bkg_zjets_pdf ,Conditional(plpdf_alternativeMELA_zjets, RooArgSet(CMS_zz4l_AlternativeMELA)),Conditional(plpdf_MELA_zjets, RooArgSet(CMS_zz4l_melaLD))); 
           
            w.import(bkg_qqzz_pdf_2D); 
            w.import(bkg_ggzz_pdf_2D); 
            w.import(bkg_zjets_pdf_2D); 
            w.import(sig_ggH_pdf_2D); 
            w.import(sig_VBF_pdf_2D);
            w.import(sig_WHi_pdf_2D); 
            w.import(sig_ZHi_pdf_2D);
            w.import(sig_ttH_pdf_2D);
            w.import(sig_ggHA_pdf_2D); 
            w.import(sig_VBFA_pdf_2D);
            w.import(sig_WHiA_pdf_2D); 
            w.import(sig_ZHiA_pdf_2D);
            w.import(sig_ttHA_pdf_2D);
        }

        
        w.writeToFile(workspace.c_str());

        delete histxsecbrggh;
        delete histxsecbrvbf;
        delete histxsecbrwhi;
        delete histxsecbrzhi;
        delete histxsecbrtth;

    }

    float getMassCut(float cardmass, bool low) {
        std::map<float, float> higgswidth;
        higgswidth[90.0] = 2.22E-03;
        higgswidth[95.0] = 2.35E-03;
        higgswidth[100.0] = 2.49E-03;
        higgswidth[105.0] = 2.65E-03;
        higgswidth[110.0] = 2.85E-03;
        higgswidth[110.5] = 2.87E-03;
        higgswidth[111.0] = 2.90E-03;
        higgswidth[111.5] = 2.92E-03;
        higgswidth[112.0] = 2.95E-03;
        higgswidth[112.5] = 2.98E-03;
        higgswidth[113.0] = 3.00E-03;
        higgswidth[113.5] = 3.03E-03;
        higgswidth[114.0] = 3.06E-03;
        higgswidth[114.5] = 3.09E-03;
        higgswidth[115.0] = 3.12E-03;
        higgswidth[115.5] = 3.16E-03;
        higgswidth[116.0] = 3.19E-03;
        higgswidth[116.5] = 3.23E-03;
        higgswidth[117.0] = 3.26E-03;
        higgswidth[117.5] = 3.30E-03;
        higgswidth[118.0] = 3.34E-03;
        higgswidth[118.5] = 3.37E-03;
        higgswidth[119.0] = 3.41E-03;
        higgswidth[119.5] = 3.46E-03;
        higgswidth[120.0] = 3.50E-03;
        higgswidth[120.5] = 3.55E-03;
        higgswidth[121.0] = 3.60E-03;
        higgswidth[121.5] = 3.65E-03;
        higgswidth[122.0] = 3.70E-03;
        higgswidth[122.5] = 3.76E-03;
        higgswidth[123.0] = 3.82E-03;
        higgswidth[123.5] = 3.88E-03;
        higgswidth[124.0] = 3.94E-03;
        higgswidth[124.5] = 4.00E-03;
        higgswidth[125.0] = 4.07E-03;
        higgswidth[125.5] = 4.14E-03;
        higgswidth[126.0] = 4.21E-03;
        higgswidth[126.5] = 4.29E-03;
        higgswidth[127.0] = 4.37E-03;
        higgswidth[127.5] = 4.45E-03;
        higgswidth[128.0] = 4.53E-03;
        higgswidth[128.5] = 4.62E-03;
        higgswidth[129.0] = 4.71E-03;
        higgswidth[129.5] = 4.81E-03;
        higgswidth[130.0] = 4.91E-03;
        higgswidth[130.5] = 5.01E-03;
        higgswidth[131.0] = 5.12E-03;
        higgswidth[131.5] = 5.23E-03;
        higgswidth[132.0] = 5.35E-03;
        higgswidth[132.5] = 5.47E-03;
        higgswidth[133.0] = 5.60E-03;
        higgswidth[133.5] = 5.74E-03;
        higgswidth[134.0] = 5.88E-03;
        higgswidth[134.5] = 6.03E-03;
        higgswidth[135.0] = 6.18E-03;
        higgswidth[135.5] = 6.34E-03;
        higgswidth[136.0] = 6.51E-03;
        higgswidth[136.5] = 6.69E-03;
        higgswidth[137.0] = 6.87E-03;
        higgswidth[137.5] = 7.06E-03;
        higgswidth[138.0] = 7.27E-03;
        higgswidth[138.5] = 7.48E-03;
        higgswidth[139.0] = 7.70E-03;
        higgswidth[139.5] = 7.93E-03;
        higgswidth[140.0] = 8.18E-03;
        higgswidth[141.0] = 8.70E-03;
        higgswidth[142.0] = 9.28E-03;
        higgswidth[143.0] = 9.93E-03;
        higgswidth[144.0] = 1.07E-02;
        higgswidth[145.0] = 1.15E-02;
        higgswidth[146.0] = 1.23E-02;
        higgswidth[147.0] = 1.34E-02;
        higgswidth[148.0] = 1.45E-02;
        higgswidth[149.0] = 1.58E-02;
        higgswidth[150.0] = 1.73E-02;
        higgswidth[151.0] = 1.91E-02;
        higgswidth[152.0] = 2.11E-02;
        higgswidth[153.0] = 2.36E-02;
        higgswidth[154.0] = 2.66E-02;
        higgswidth[155.0] = 3.03E-02;
        higgswidth[156.0] = 3.51E-02;
        higgswidth[157.0] = 4.14E-02;
        higgswidth[158.0] = 5.02E-02;
        higgswidth[159.0] = 6.32E-02;
        higgswidth[160.0] = 8.31E-02;
        higgswidth[162.0] = 1.47E-01;
        higgswidth[164.0] = 2.15E-01;
        higgswidth[165.0] = 2.46E-01;
        higgswidth[166.0] = 2.76E-01;
        higgswidth[168.0] = 3.30E-01;
        higgswidth[170.0] = 3.80E-01;
        higgswidth[172.0] = 4.29E-01;
        higgswidth[174.0] = 4.77E-01;
        higgswidth[175.0] = 5.01E-01;
        higgswidth[176.0] = 5.25E-01;
        higgswidth[178.0] = 5.75E-01;
        higgswidth[180.0] = 6.31E-01;
        higgswidth[182.0] = 7.00E-01;
        higgswidth[184.0] = 7.88E-01;
        higgswidth[185.0] = 8.32E-01;
        higgswidth[186.0] = 8.76E-01;
        higgswidth[188.0] = 9.60E-01;
        higgswidth[190.0] = 1.04E+00;
        higgswidth[192.0] = 1.12E+00;
        higgswidth[194.0] = 1.20E+00;
        higgswidth[195.0] = 1.24E+00;
        higgswidth[196.0] = 1.28E+00;
        higgswidth[198.0] = 1.35E+00;
        higgswidth[200.0] = 1.43E+00;
        higgswidth[202.0] = 1.51E+00;
        higgswidth[204.0] = 1.59E+00;
        higgswidth[206.0] = 1.68E+00;
        higgswidth[208.0] = 1.76E+00;
        higgswidth[210.0] = 1.85E+00;
        higgswidth[212.0] = 1.93E+00;
        higgswidth[214.0] = 2.02E+00;
        higgswidth[216.0] = 2.12E+00;
        higgswidth[218.0] = 2.21E+00;
        higgswidth[220.0] = 2.31E+00;
        higgswidth[222.0] = 2.40E+00;
        higgswidth[224.0] = 2.50E+00;
        higgswidth[226.0] = 2.61E+00;
        higgswidth[228.0] = 2.71E+00;
        higgswidth[230.0] = 2.82E+00;
        higgswidth[232.0] = 2.93E+00;
        higgswidth[234.0] = 3.04E+00;
        higgswidth[236.0] = 3.16E+00;
        higgswidth[238.0] = 3.27E+00;
        higgswidth[240.0] = 3.40E+00;
        higgswidth[242.0] = 3.52E+00;
        higgswidth[244.0] = 3.64E+00;
        higgswidth[246.0] = 3.77E+00;
        higgswidth[248.0] = 3.91E+00;
        higgswidth[250.0] = 4.04E+00;
        higgswidth[252.0] = 4.18E+00;
        higgswidth[254.0] = 4.32E+00;
        higgswidth[256.0] = 4.46E+00;
        higgswidth[258.0] = 4.61E+00;
        higgswidth[260.0] = 4.76E+00;
        higgswidth[262.0] = 4.91E+00;
        higgswidth[264.0] = 5.07E+00;
        higgswidth[266.0] = 5.23E+00;
        higgswidth[268.0] = 5.39E+00;
        higgswidth[270.0] = 5.55E+00;
        higgswidth[272.0] = 5.72E+00;
        higgswidth[274.0] = 5.89E+00;
        higgswidth[276.0] = 6.07E+00;
        higgswidth[278.0] = 6.25E+00;
        higgswidth[280.0] = 6.43E+00;
        higgswidth[282.0] = 6.61E+00;
        higgswidth[284.0] = 6.80E+00;
        higgswidth[286.0] = 6.99E+00;
        higgswidth[288.0] = 7.19E+00;
        higgswidth[290.0] = 7.39E+00;
        higgswidth[295.0] = 7.90E+00;
        higgswidth[300.0] = 8.43E+00;
        higgswidth[305.0] = 8.99E+00;
        higgswidth[310.0] = 9.57E+00;
        higgswidth[315.0] = 1.02E+01;
        higgswidth[320.0] = 1.08E+01;
        higgswidth[325.0] = 1.14E+01;
        higgswidth[330.0] = 1.21E+01;
        higgswidth[335.0] = 1.28E+01;
        higgswidth[340.0] = 1.35E+01;
        higgswidth[345.0] = 1.42E+01;
        higgswidth[350.0] = 1.52E+01;
        higgswidth[360.0] = 1.76E+01;
        higgswidth[370.0] = 2.02E+01;
        higgswidth[380.0] = 2.31E+01;
        higgswidth[390.0] = 2.61E+01;
        higgswidth[400.0] = 2.92E+01;
        higgswidth[410.0] = 3.25E+01;
        higgswidth[420.0] = 3.59E+01;
        higgswidth[430.0] = 3.94E+01;
        higgswidth[440.0] = 4.30E+01;
        higgswidth[450.0] = 4.68E+01;
        higgswidth[460.0] = 5.08E+01;
        higgswidth[470.0] = 5.49E+01;
        higgswidth[480.0] = 5.91E+01;
        higgswidth[490.0] = 6.35E+01;
        higgswidth[500.0] = 6.80E+01;
        higgswidth[510.0] = 7.27E+01;
        higgswidth[520.0] = 7.75E+01;
        higgswidth[530.0] = 8.25E+01;
        higgswidth[540.0] = 8.77E+01;
        higgswidth[550.0] = 9.30E+01;
        higgswidth[560.0] = 9.86E+01;
        higgswidth[570.0] = 1.04E+02;
        higgswidth[580.0] = 1.10E+02;
        higgswidth[590.0] = 1.16E+02;
        higgswidth[600.0] = 1.23E+02;
        higgswidth[610.0] = 1.29E+02;
        higgswidth[620.0] = 1.36E+02;
        higgswidth[630.0] = 1.43E+02;
        higgswidth[640.0] = 1.50E+02;
        higgswidth[650.0] = 1.58E+02;
        higgswidth[660.0] = 1.65E+02;
        higgswidth[670.0] = 1.73E+02;
        higgswidth[680.0] = 1.82E+02;
        higgswidth[690.0] = 1.90E+02;
        higgswidth[700.0] = 1.99E+02;
        higgswidth[710.0] = 2.08E+02;
        higgswidth[720.0] = 2.17E+02;
        higgswidth[730.0] = 2.27E+02;
        higgswidth[740.0] = 2.37E+02;
        higgswidth[750.0] = 2.47E+02;
        higgswidth[760.0] = 2.58E+02;
        higgswidth[770.0] = 2.69E+02;
        higgswidth[780.0] = 2.80E+02;
        higgswidth[790.0] = 2.92E+02;
        higgswidth[800.0] = 3.04E+02;
        higgswidth[810.0] = 3.17E+02;
        higgswidth[820.0] = 3.30E+02;
        higgswidth[830.0] = 3.43E+02;
        higgswidth[840.0] = 3.57E+02;
        higgswidth[850.0] = 3.71E+02;
        higgswidth[860.0] = 3.86E+02;
        higgswidth[870.0] = 4.01E+02;
        higgswidth[880.0] = 4.16E+02;
        higgswidth[890.0] = 4.32E+02;
        higgswidth[900.0] = 4.49E+02;
        higgswidth[910.0] = 4.66E+02;
        higgswidth[920.0] = 4.84E+02;
        higgswidth[930.0] = 5.02E+02;
        higgswidth[940.0] = 5.21E+02;
        higgswidth[950.0] = 5.40E+02;
        higgswidth[960.0] = 5.60E+02;
        higgswidth[970.0] = 5.81E+02;
        higgswidth[980.0] = 6.02E+02;
        higgswidth[990.0] = 6.24E+02;
        higgswidth[1000.0] = 6.47E+02;


        double windowVal = max(higgswidth[cardmass], float(1.));
        double lowside = (cardmass >= 275) ? 180. : 100.;
        if (low) return std::max((cardmass - 20.*windowVal), lowside);
        else return std::min((cardmass + 15.*windowVal), 1600.);
    }


    void makeCards() {
        for (float i = 114.; i <= 160.; i += 1.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
        for (float i = 162.; i <= 290.; i += 2.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
        for (float i = 295.; i <= 350.; i += 5.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
        for (float i = 360.; i <= 400.; i += 10.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
        for (float i = 420.; i <= 1000.; i += 20.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
    }

    void makeCard(float i) {
        createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
        createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
        createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
    }


};

void makeSpinCPMeasurementDatacards(Int_t AlternativeSignalType = 1) {

    string alternativeMELAFilename = "";
    if (AlternativeSignalType == 1) {
      alternativeMELAFilename = "PSMelaVsMass2DShapes.root";
    } 
    if (AlternativeSignalType == 2) {
      alternativeMELAFilename = "SpinTwoMinimalMelaVsMass2DShapes.root";
    } 


    HiggsMassPointInfo hmpi7;
    hmpi7.lumi = 5.05;
    hmpi7.base_lumi = 5.05;
    hmpi7.z1min = 40.;
    hmpi7.z2min = 12.;
    hmpi7.massLowBkgFit = 100.;
    hmpi7.massHighBkgFit = 1600.;
    hmpi7.melacut = -1.0;
    hmpi7.do1D = true;
    hmpi7.do7TeV = true;
    hmpi7.doFFT = true;
    hmpi7.doMassError = false;
    hmpi7.alternativeSignalType = AlternativeSignalType;
    hmpi7.treeFolder = "/data1/sixie/ntuples/step3new/hcp/2011/";
    hmpi7.melafilename = "mela2DShapes.root";
//     hmpi7.melafilename = "PSMelaVsMass2DShapes.root";
    hmpi7.alternativeMELAfilename = alternativeMELAFilename;

    init(hmpi7.do7TeV);

    hmpi7.ymaker_data.fill(hmpi7.treeFolder+"hcp.root");
    hmpi7.do1D = false;
    hmpi7.makeCard(125.);

   

    HiggsMassPointInfo hmpi8;
    hmpi8.lumi = 12.2; 
    hmpi8.base_lumi = 12.2; 
    hmpi8.z1min = 40.;
    hmpi8.z2min = 12.;
    hmpi8.massLowBkgFit = 100.;
    hmpi8.massHighBkgFit = 1600.;
    hmpi8.melacut = -1.0;
    hmpi8.do1D = true;
    hmpi8.do7TeV = false;
    hmpi8.doFFT = true;
    hmpi8.doMassError = false;
    hmpi8.alternativeSignalType = AlternativeSignalType;
    hmpi8.treeFolder = "/data1/sixie/ntuples/step3new/hcp/2012/";
    hmpi8.melafilename = "mela2DShapes.root";
//     hmpi8.melafilename = "PSMelaVsMass2DShapes.root";
    hmpi8.alternativeMELAfilename = alternativeMELAFilename;

    init(hmpi8.do7TeV);


    hmpi8.ymaker_data.fill(hmpi8.treeFolder+"hcp.root");

    hmpi8.do1D = false;
    hmpi8.makeCard(125.);

}


