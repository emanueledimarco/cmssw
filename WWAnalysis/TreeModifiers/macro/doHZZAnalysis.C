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

float getZXSystematicsUp(bool is7, int ch) {

    /*
    if (is7) {
        if      (ch == 0) return 1.0 + 0.631941/0.968059;
        else if (ch == 1) return 1.0 + 1.67584/1.14416;
        else              return 1.0 + 1.37835/2.20165;
    }
    else {
        if      (ch == 0) return 1.0 + 1.29856/2.97144;
        else if (ch == 1) return 1.0 + 3.54288/6.61712;
        else              return 1.0 + 4.28711/9.42289;
    }
    */

    return 1.4;
}

float getZXSystematicsDown(bool is7, int ch) {

    /*
    if (is7) {
        if      (ch == 0) return 1.0 - 0.248059/0.968059;
        else if (ch == 1) return 1.0 - 0.884161/1.14416;
        else              return 1.0 - 1.89165/2.20165;
    }
    else {
        if      (ch == 0) return 1.0 - 2.16144/2.97144;
        else if (ch == 1) return 1.0 - 3.24712/6.61712;
        else              return 1.0 - 5.26289/9.42289;
    }
    */

    return 0.6;
}


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
    std::string treeFolder;
    std::string melafilename;
    
    DataYieldMaker ymaker_data;


    void createCard(float mass, float massLow, float massHigh, int ch) {

        if (mass > 180.) doFFT = true;

        std::string chstr;
        if (ch == 0) chstr = "4mu";
        if (ch == 1) chstr = "4e";
        if (ch == 2) chstr = "2e2mu";
        if (ch == 3) chstr = "2mu2e";

        stringstream mass_str_ss;
        mass_str_ss << mass;
        std::string mass_str = mass_str_ss.str();
        
        std::cout << "Creating datacard for " << mass_str << " GeV mass point and channel " << chstr << " ... " << std::endl;
        
        std::string card_name   = do1D ? (std::string("card_1D_m")+mass_str+"_"+(do7TeV ? "7TeV_" : "8TeV_")) : (std::string("card_2D_m")+mass_str+"_"+(do7TeV ? "7TeV_" : "8TeV_"));
        card_name += chstr;
        std::string workspace = card_name+"_workspace.root";
        
        float yield_dt  = ymaker_data.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
         
        RooWorkspace w("w", "");
        
        RooRealVar CMS_zz4l_melaLD ("CMS_zz4l_melaLD" , "MELA"            ,   0    , 1       , "");
        RooRealVar CMS_zz4l_mass_1D("CMS_zz4l_mass_1D", "M(4l)"           , massLow, massHigh, "GeV/c^{2}");
        RooRealVar CMS_zz4l_massErr("CMS_zz4l_massErr", "M(4l) rel. Error",   0.002, 0.2     , "");
        CMS_zz4l_massErr.setBins(30);
        
        if (doFFT) CMS_zz4l_mass_1D.setBins(10000);
        
        if (do1D) {
          if(!doMassError) {
            RooArgSet argset_obs(CMS_zz4l_mass_1D, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
            
            ymaker_data.getDataSet1D(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D);
            
            w.import(data_obs);
          } else {
            RooArgSet argset_obs(CMS_zz4l_mass_1D, CMS_zz4l_massErr, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
            
            ymaker_data.getDataSet1DEBE(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_massErr);
            
            w.import(data_obs);            
          }
        }
        
        else {
          if(!doMassError) {
            RooArgSet argset_obs(CMS_zz4l_mass_1D, CMS_zz4l_melaLD, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
            
            ymaker_data.getDataSet2D(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            
            w.import(data_obs);
          } else {
            RooArgSet argset_obs(CMS_zz4l_mass_1D, CMS_zz4l_massErr, CMS_zz4l_melaLD, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
            
            ymaker_data.getDataSet2DEBE(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_massErr, CMS_zz4l_melaLD);
            
            w.import(data_obs);
          }
        }
        
        
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

        // ebe
        float ebe_qqldm = 0.;
        float ebe_qqlds = 0.;
        float ebe_qqlnk = 0.;
        float ebe_qqgam = 0.;
        float ebe_qqgas = 0.;
        float ebe_qqf   = 0.;

        float ebe_ggldm = 0.;
        float ebe_gglds = 0.;
        float ebe_gglnk = 0.;
        float ebe_gggam = 0.;
        float ebe_gggas = 0.;
        float ebe_ggf   = 0.;

        float ebe_zxldm = 0.;
        float ebe_zxlds = 0.;
        float ebe_zxlnk = 0.;
        float ebe_zxgam = 0.;
        float ebe_zxgas = 0.;
        float ebe_zxf   = 0.;

        if (do7TeV) {

            if (ch == 0) {
                zxme  = 135.829;
                zxsi = 24.5744;


                qqa0  = 107.564;
                qqa1  = 7.30727;
                qqa2  = 112.302;
                qqa3  = 0.0437908;
                qqa4  = 185.487;
                qqa5  = 8.55048;
                qqa6  = 35.2851;
                qqa7  = 0.0965921;
                qqa8  = 59.2775;
                qqa9  = 0.0506361;
                qqa10 = 99.5362;
                qqa11 = -7.24557;
                qqa12 = 1377.43;
                qqa13 = 0.0860447;


                gga0  = 153.748;
                gga1  = 64.0604;
                gga2  = 111.833;
                gga3  = 0.0460111;
                gga4  = 183.815;
                gga5  = 8.78774;
                gga6  = 41.9141;
                gga7  = 0.443452;
                gga8  = 40.531;
                gga9  = -0.208836;

                
                ebe_qqldm = 0.00836569;
                ebe_qqlds = 0.00159688;
                ebe_qqlnk = 0.726805;
                ebe_qqf   = 0.29648;
                
                ebe_ggldm = 0.00777304;
                ebe_gglds = 0.00124103;
                ebe_gglnk = 0.750389;
                ebe_ggf   = 0.423512;
                
                ebe_zxldm = 0.00699555;
                ebe_zxlds = 0.000605924;
                ebe_zxlnk = 1.45644;
                ebe_zxf   = 0.999984;
            }

            else if (ch == 1) {
                zxme  = 129.593;
                zxsi = 25.1211;


                qqa0  = 112.053;
                qqa1  = 25.6835;
                qqa2  = 120.663;
                qqa3  = 0.0396448;
                qqa4  = 184.122;
                qqa5  = 10.9143;
                qqa6  = 34.9035;
                qqa7  = 0.113666;
                qqa8  = 58.6916;
                qqa9  = 0.0588379;
                qqa10 = 95.1808;
                qqa11 = -7.05076;
                qqa12 = 1788.69;
                qqa13 = 0.0772151;


                gga0  = 164.332;
                gga1  = 58.5182;
                gga2  = 105.699;
                gga3  = 0.0521393;
                gga4  = 182.742;
                gga5  = 10.6714;
                gga6  = 42.681;
                gga7  = 0.44492;
                gga8  = 44.8908;
                gga9  = -0.198437;

                ebe_qqldm = 0.00931315;
                ebe_qqlds = 0.00140546;
                ebe_qqgam = 0.0190073;
                ebe_qqgas = 0.00504234;
                ebe_qqf   = 0.506115;
                
                ebe_ggldm = 0.00980677;
                ebe_gglds = 0.001633;
                ebe_gggam = 0.019325;
                ebe_gggas = 0.00449822;
                ebe_ggf   = 0.676257;
                
                ebe_zxldm = 0.0149304;
                ebe_zxlds = 0.0037406;
                ebe_zxgam = 0.0333999;
                ebe_zxgas = 0.0112965;
                ebe_zxf   = 0.662721;
            }

            else if (ch == 2) {
                zxme  = 142.496;
                zxsi = 21.878;


                qqa0  = 108.924;
                qqa1  = 15.594;
                qqa2  = 117.297;
                qqa3  = 0.0402424;
                qqa4  = 184.601;
                qqa5  = 9.84846;
                qqa6  = 35.3433;
                qqa7  = 0.109418;
                qqa8  = 60.1712;
                qqa9  = 0.0591363;
                qqa10 = 95.0982;
                qqa11 = -6.36596;
                qqa12 = 1174.68;
                qqa13 = 0.0913671;


                gga0  = 162.829;
                gga1  = 54.3068;
                gga2  = 109.33;
                gga3  = 0.0521449;
                gga4  = 183.862;
                gga5  = 10.2186;
                gga6  = 42.7888;
                gga7  = 0.443823;
                gga8  = 44.8107;
                gga9  = -0.195945;


                ebe_qqldm = 0.0075814;
                ebe_qqlds = 0.00112148;
                ebe_qqgam = 0.0163012;
                ebe_qqgas = 0.00463604;
                ebe_qqf   = 0.595006;
                
                ebe_ggldm = 0.00766245;
                ebe_gglds = 0.00101953;
                ebe_gggas = 0.015324;
                ebe_gggas = 0.00511548;
                ebe_ggf   = 0.62782;
                
                ebe_zxldm = 0.0103013;
                ebe_zxlds = 0.0013906;
                ebe_zxgam = 0.0256454;
                ebe_zxgas = 0.00820503;
                ebe_zxf   = 0.643808;
            }

            else {
                zxme  = 142.496;
                zxsi = 21.878;


                qqa0  = 108.924;
                qqa1  = 15.594;
                qqa2  = 117.297;
                qqa3  = 0.0402424;
                qqa4  = 184.601;
                qqa5  = 9.84846;
                qqa6  = 35.3433;
                qqa7  = 0.109418;
                qqa8  = 60.1712;
                qqa9  = 0.0591363;
                qqa10 = 95.0982;
                qqa11 = -6.36596;
                qqa12 = 1174.68;
                qqa13 = 0.0913671;


                gga0  = 162.829;
                gga1  = 54.3068;
                gga2  = 109.33;
                gga3  = 0.0521449;
                gga4  = 183.862;
                gga5  = 10.2186;
                gga6  = 42.7888;
                gga7  = 0.443823;
                gga8  = 44.8107;
                gga9  = -0.195945;


                ebe_qqldm = 0.0075814;
                ebe_qqlds = 0.00112148;
                ebe_qqgam = 0.0163012;
                ebe_qqgas = 0.00463604;
                ebe_qqf   = 0.595006;

                ebe_ggldm = 0.00766245;
                ebe_gglds = 0.00101953;
                ebe_gggas = 0.015324;
                ebe_gggas = 0.00511548;
                ebe_ggf   = 0.62782;

                ebe_zxldm = 0.0103013;
                ebe_zxlds = 0.0013906;
                ebe_zxgam = 0.0256454;
                ebe_zxgas = 0.00820503;
                ebe_zxf   = 0.643808;
            }

        }

        else {
            if (ch == 0) {
                zxme  = 141.507;
                zxsi  = 21.6726;

                qqa0  = 107.25;
                qqa1  = 8.57192;
                qqa2  = 120.212;
                qqa3  = 0.0454114;
                qqa4  = 185.411;
                qqa5  = 8.65822;
                qqa6  = 40.5369;
                qqa7  = 0.0991484;
                qqa8  = 48.5337;
                qqa9  = 0.0450332;
                qqa10 = 100.644;
                qqa11 = -8.71111;
                qqa12 = 1866.37;
                qqa13 = 0.0656402;


                gga0  = 109.444;
                gga1  = 71.2501;
                gga2  = 146.592;
                gga3  = 0.048661;
                gga4  = 185.434;
                gga5  = 9.69742;
                gga6  = 44.6195;
                gga7  = 0.493789;
                gga8  = 40.4904;
                gga9  = -0.155732;


                ebe_qqldm = 0.00850474;
                ebe_qqlds = 0.00158463;
                ebe_qqlnk = 1.34272;
                ebe_qqf   = 0.502923;

                ebe_ggldm = 0.00817017;
                ebe_gglds = 0.00113282;
                ebe_gglnk = 0.723619;
                ebe_ggf   = 0.590413;

                ebe_zxldm = 0.00854642;
                ebe_zxlds = 0.00132492;
                ebe_zxlnk = 0.821483;
                ebe_zxf   = 0.99839;
            }

            else if (ch == 1) {
                zxme  = 133.203;
                zxsi  = 21.9732;

                qqa0  = 113.21;
                qqa1  = 18.7447;
                qqa2  = 124.458;
                qqa3  = 0.0322816;
                qqa4  = 184.45;
                qqa5  = 10.4299;
                qqa6  = 34.5919;
                qqa7  = 0.0883106;
                qqa8  = 58.1311;
                qqa9  = 0.0612426;
                qqa10 = 97.8627;
                qqa11 = -7.36006;
                qqa12 = 1594.26;
                qqa13 = 0.0762025;


                gga0  = 133.555;
                gga1  = 52.092;
                gga2  = 143.361;
                gga3  = 0.0396186;
                gga4  = 184.579;
                gga5  = 10.3012;
                gga6  = 41.1703;
                gga7  = 0.496282;
                gga8  = 34.9484;
                gga9  = -0.213565;


                ebe_qqldm = 0.0106893;
                ebe_qqlds = 0.00170811;
                ebe_qqgam = 0.0197709;
                ebe_qqgas = 0.00449331;
                ebe_qqf   = 0.596535;
                
                ebe_ggldm = 0.0102151;
                ebe_gglds = 0.00145494;
                ebe_gggam = 0.0191593;
                ebe_gggas = 0.00480876;
                ebe_ggf   = 0.615308;
                
                ebe_zxldm = 0.0149392;
                ebe_zxlds = 0.0038365;
                ebe_zxgam = 0.0499944;
                ebe_zxgas = 0.0391288;
                ebe_zxf   = 0.353372;
            }
        
            else if (ch == 2) {
                zxme  = 127.551;
                zxsi  = 21.4145;

                qqa0  = 106.153;
                qqa1  = 8.19048;
                qqa2  = 124.026;
                qqa3  = 0.0447146;
                qqa4  = 185.25;
                qqa5  = 8.84078;
                qqa6  = 38.6016;
                qqa7  = 0.0870263;
                qqa8  = 53.7872;
                qqa9  = 0.0448448;
                qqa10 = 96.0439;
                qqa11 = -5.57062;
                qqa12 = 1467.78;
                qqa13 = 0.166461;


                gga0  = 137.284;
                gga1  = 43.5852;
                gga2  = 150.161;
                gga3  = 0.0455927;
                gga4  = 185.97;
                gga5  = 10.8979;
                gga6  = 44.84;
                gga7  = 0.500988;
                gga8  = 41.6484;
                gga9  = -0.165554;


                ebe_qqldm = 0.00797046;
                ebe_qqlds = 0.000902227;
                ebe_qqgam = 0.0181578;
                ebe_qqgas = 0.00320777;
                ebe_qqf   = 0.698307;

                ebe_ggldm = 0.00779737;
                ebe_gglds = 0.00102945;
                ebe_gggam = 0.0170503;
                ebe_gggas = 0.00371655;
                ebe_ggf   = 0.664479;
                
                ebe_zxldm = 0.00728841;
                ebe_zxlds = 0.000587849;
                ebe_zxgam = 0.0183173;
                ebe_zxgas = 0.00255046;
                ebe_zxf   = 0.459353;
            }

            else {
                zxme  = 126.093;
                zxsi  = 20.0526;

                qqa0  = 116.147;
                qqa1  = 17.7373;
                qqa2  = 125.913;
                qqa3  = 0.0348915;
                qqa4  = 185.392;
                qqa5  = 10.4563;
                qqa6  = 33.8311;
                qqa7  = 0.130637;
                qqa8  = 59.8135;
                qqa9  = 0.0710823;
                qqa10 = 95.4298;
                qqa11 = -6.95096;
                qqa12 = 817.716;
                qqa13 = 0.104338;


                gga0  = 170.057;
                gga1  = 70.965;
                gga2  = 124.605;
                gga3  = 0.0437131;
                gga4  = 185.146;
                gga5  = 9.79001;
                gga6  = 43.8616;
                gga7  = 0.494783;
                gga8  = 44.36;
                gga9  = -0.184552;


                ebe_qqldm = 0.00971269;
                ebe_qqlds = 0.00169118;
                ebe_qqgam = 0.0157131;
                ebe_qqgas = 0.00372416;
                ebe_qqf   = 0.700673;

                ebe_ggldm = 0.00959196;
                ebe_gglds = 0.00129222;
                ebe_gggam = 0.0150871;
                ebe_gggas = 0.00436361;
                ebe_ggf   = 0.481612;

                ebe_zxldm = 0.0149989;
                ebe_zxlds = 0.0039411;
                ebe_zxgam = 0.0499894;
                ebe_zxgas = 0.0271768;
                ebe_zxf   = 0.461463;
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
        
        RooRealVar sig_mean_err_4mu (("sig_mean_err_4mu" +tevstr).c_str()          , "", 0., -10., 10.);
        RooRealVar sig_mean_err_4e  (("sig_mean_err_4e"  +tevstr).c_str()          , "", 0., -10., 10.);
        RooRealVar sig_sigma_err_4mu(("sig_sigma_err_4mu"+tevstr).c_str()          , "", 0., -10., 10.);
        RooRealVar sig_sigma_err_4e (("sig_sigma_err_4e" +tevstr).c_str()          , "", 0., -10., 10.);
        RooRealVar sig_gamma_err    ("sig_gamma_err"                               , "", 0., -10., 10.);

        RooRealVar qqzz_ebe_LdM(("bkg_qqzz_"+chstr+tevstr+"_ebe_LdM").c_str(), "EBE Landau Mean", ebe_qqldm);
        RooRealVar ggzz_ebe_LdM(("bkg_ggzz_"+chstr+tevstr+"_ebe_LdM").c_str(), "EBE Landau Mean", ebe_ggldm);
        RooRealVar zjets_ebe_LdM(("bkg_zjets_"+chstr+tevstr+"_ebe_LdM").c_str(), "EBE Landau Mean", ebe_zxldm);

        RooRealVar qqzz_ebe_LdS(("bkg_qqzz_"+chstr+tevstr+"_ebe_LdS").c_str(), "EBE Landau Sigma", ebe_qqlds);
        RooRealVar ggzz_ebe_LdS(("bkg_ggzz_"+chstr+tevstr+"_ebe_LdS").c_str(), "EBE Landau Sigma", ebe_gglds);
        RooRealVar zjets_ebe_LdS(("bkg_zjets_"+chstr+tevstr+"_ebe_LdS").c_str(), "EBE Landau Sigma", ebe_zxlds);

        RooRealVar qqzz_ebe_LnK(("bkg_qqzz_"+chstr+tevstr+"_ebe_LnK").c_str(), "EBE LogNormal Kappa", ebe_qqlnk);
        RooRealVar ggzz_ebe_LnK(("bkg_ggzz_"+chstr+tevstr+"_ebe_LnK").c_str(), "EBE LogNormal Kappa", ebe_gglnk);
        RooRealVar zjets_ebe_LnK(("bkg_zjets_"+chstr+tevstr+"_ebe_LnK").c_str(), "EBE LogNormal Kappa", ebe_zxlnk);

        RooRealVar qqzz_ebe_GaM(("bkg_qqzz_"+chstr+tevstr+"_ebe_GaM").c_str(), "EBE Gaussian Mean", ebe_qqgam);
        RooRealVar ggzz_ebe_GaM(("bkg_ggzz_"+chstr+tevstr+"_ebe_GaM").c_str(), "EBE Gaussian Mean", ebe_gggam);
        RooRealVar zjets_ebe_GaM(("bkg_zjets_"+chstr+tevstr+"_ebe_GaM").c_str(), "EBE Gaussian Mean", ebe_zxgam);

        RooRealVar qqzz_ebe_GaS(("bkg_qqzz_"+chstr+tevstr+"_ebe_GaS").c_str(), "EBE Gaussian Sigma", ebe_qqgas);
        RooRealVar ggzz_ebe_GaS(("bkg_ggzz_"+chstr+tevstr+"_ebe_GaS").c_str(), "EBE Gaussian Sigma", ebe_gggas);
        RooRealVar zjets_ebe_GaS(("bkg_zjets_"+chstr+tevstr+"_ebe_GaS").c_str(), "EBE Gaussian Sigma", ebe_zxgas);

        RooRealVar qqzz_ebe_frac(("bkg_qqzz_"+chstr+tevstr+"_ebe_frac").c_str(), "EBE Landau Fraction", ebe_qqf);
        RooRealVar ggzz_ebe_frac(("bkg_ggzz_"+chstr+tevstr+"_ebe_frac").c_str(), "EBE Landau Fraction", ebe_ggf);
        RooRealVar zjets_ebe_frac(("bkg_zjets_"+chstr+tevstr+"_ebe_frac").c_str(), "EBE Landau Fraction", ebe_zxf);

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
                gghxsecbry[(int)(i-110.)] = getXsecggHByChannel(i, ch, true);
                vbfxsecbry[(int)(i-110.)] = getXsecVBFByChannel(i, ch, true);
                whixsecbry[(int)(i-110.)] = getXsecWHiByChannel(i, ch, true);
                zhixsecbry[(int)(i-110.)] = getXsecZHiByChannel(i, ch, true);
                tthxsecbry[(int)(i-110.)] = getXsecttHByChannel(i, ch, true);
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
                gghxsecbry[((int)(i-160.))/2] = getXsecggHByChannel(i, ch, true);
                vbfxsecbry[((int)(i-160.))/2] = getXsecVBFByChannel(i, ch, true);
                whixsecbry[((int)(i-160.))/2] = getXsecWHiByChannel(i, ch, true);
                zhixsecbry[((int)(i-160.))/2] = getXsecZHiByChannel(i, ch, true);
                tthxsecbry[((int)(i-160.))/2] = getXsecttHByChannel(i, ch, true);
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
                gghxsecbry[((int)(i-290.))/5] = getXsecggHByChannel(i, ch, true);
                vbfxsecbry[((int)(i-290.))/5] = getXsecVBFByChannel(i, ch, true);
                whixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecWHiByChannel(i, ch, true) : 0.0;
                zhixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecZHiByChannel(i, ch, true) : 0.0;
                tthxsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecttHByChannel(i, ch, true) : 0.0;
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
                gghxsecbry[((int)(i-350.))/10] = getXsecggHByChannel(i, ch, true);
                vbfxsecbry[((int)(i-350.))/10] = getXsecVBFByChannel(i, ch, true);
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
                gghxsecbry[((int)(i-400.))/20] = getXsecggHByChannel(i, ch, true);
                vbfxsecbry[((int)(i-400.))/20] = getXsecVBFByChannel(i, ch, true);
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

        RooHistFunc ggh_xsecbr(("xsecbr_ggH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrggh);
        RooHistFunc vbf_xsecbr(("xsecbr_VBF_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrvbf);
        RooHistFunc whi_xsecbr(("xsecbr_WHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrwhi);
        RooHistFunc zhi_xsecbr(("xsecbr_ZHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrzhi);
        RooHistFunc tth_xsecbr(("xsecbr_ttH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrtth);

        RooArgList* sig_mean_err_al;
        if      (ch == 0) sig_mean_err_al = new RooArgList(masshiggs, sig_mean_err_4mu);
        else if (ch == 1) sig_mean_err_al = new RooArgList(masshiggs, sig_mean_err_4e);
        else              sig_mean_err_al = new RooArgList(masshiggs, sig_mean_err_4mu, sig_mean_err_4e);

        RooArgList* sig_sigma_err_al;
        if      (ch == 0) sig_sigma_err_al = new RooArgList(masshiggs, sig_sigma_err_4mu);
        else if (ch == 1) sig_sigma_err_al = new RooArgList(masshiggs, sig_sigma_err_4e);
        else              sig_sigma_err_al = new RooArgList(masshiggs, sig_sigma_err_4mu, sig_sigma_err_4e);

        RooFormulaVar ggh_mean_CB  (("sig_ggh_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch, do7TeV, doFFT).c_str()           , *sig_mean_err_al);
        RooFormulaVar ggh_sigma_CB (("sig_ggh_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str()                  , *sig_sigma_err_al);
        RooFormulaVar ggh_alphaL   (("sig_ggh_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str()                            , RooArgList(masshiggs));
        RooFormulaVar ggh_alphaR   (("sig_ggh_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str()                            , RooArgList(masshiggs));
        RooFormulaVar ggh_nL       (("sig_ggh_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()                                , RooArgList(masshiggs));
        RooFormulaVar ggh_nR       (("sig_ggh_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()                                , RooArgList(masshiggs));
        RooFormulaVar ggh_gamma_BW (("sig_ggh_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str()                             , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar ggh_sigma_EBE(("sig_ggh_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                                                          , RooArgList(CMS_zz4l_massErr, ggh_sigma_CB));
        RooFormulaVar ggh_sigma_EBM(("sig_ggh_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar ggh_norm     ("ggH_norm"    , ("@0*@1*@2*"+lumistr).c_str(), RooArgList(cs_scale_z2_ggh, ggh_xsecbr, yield_var_ggh));

        RooFormulaVar vbf_mean_CB  (("sig_vbf_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , *sig_mean_err_al);
        RooFormulaVar vbf_sigma_CB (("sig_vbf_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , *sig_sigma_err_al);
        RooFormulaVar vbf_alphaL   (("sig_vbf_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar vbf_alphaR   (("sig_vbf_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar vbf_nL       (("sig_vbf_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar vbf_nR       (("sig_vbf_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar vbf_gamma_BW (("sig_vbf_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar vbf_sigma_EBE(("sig_vbf_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, vbf_sigma_CB));
        RooFormulaVar vbf_sigma_EBM(("sig_vbf_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar vbf_norm     ("qqH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_vbf, vbf_xsecbr, yield_var_vbf));

        RooFormulaVar whi_mean_CB  (("sig_whi_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch, do7TeV, doFFT).c_str() , *sig_mean_err_al);
        RooFormulaVar whi_sigma_CB (("sig_whi_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , *sig_sigma_err_al);
        RooFormulaVar whi_alphaL   (("sig_whi_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar whi_alphaR   (("sig_whi_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar whi_nL       (("sig_whi_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar whi_nR       (("sig_whi_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar whi_gamma_BW (("sig_whi_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar whi_sigma_EBE(("sig_whi_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, whi_sigma_CB));
        RooFormulaVar whi_sigma_EBM(("sig_whi_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar whi_norm     ("WH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_whi, whi_xsecbr, yield_var_whi));

        RooFormulaVar zhi_mean_CB  (("sig_zhi_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , *sig_mean_err_al);
        RooFormulaVar zhi_sigma_CB (("sig_zhi_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , *sig_sigma_err_al);
        RooFormulaVar zhi_alphaL   (("sig_zhi_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar zhi_alphaR   (("sig_zhi_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar zhi_nL       (("sig_zhi_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar zhi_nR       (("sig_zhi_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar zhi_gamma_BW (("sig_zhi_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar zhi_sigma_EBE(("sig_zhi_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, zhi_sigma_CB));
        RooFormulaVar zhi_sigma_EBM(("sig_zhi_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar zhi_norm     ("ZH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_zhi, zhi_xsecbr, yield_var_zhi));

        RooFormulaVar tth_mean_CB  (("sig_tth_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch, do7TeV, doFFT).c_str()  , *sig_mean_err_al);
        RooFormulaVar tth_sigma_CB (("sig_tth_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch, do7TeV).c_str() , *sig_sigma_err_al);
        RooFormulaVar tth_alphaL   (("sig_tth_"+chstr+tevstr+"_alphaL"   ).c_str(), getSignalCBAlphaLString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar tth_alphaR   (("sig_tth_"+chstr+tevstr+"_alphaR"   ).c_str(), getSignalCBAlphaRString(mass, ch, do7TeV).c_str(), RooArgList(masshiggs));
        RooFormulaVar tth_nL       (("sig_tth_"+chstr+tevstr+"_nL"       ).c_str(), getSignalCBNLString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar tth_nR       (("sig_tth_"+chstr+tevstr+"_nR"       ).c_str(), getSignalCBNRString(mass, ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar tth_gamma_BW (("sig_tth_"+chstr+tevstr+"_gamma_BW" ).c_str(), getSignalBWGammaString(mass, ch, do7TeV).c_str() , RooArgList(masshiggs, sig_gamma_err));
        RooFormulaVar tth_sigma_EBE(("sig_tth_"+chstr+tevstr+"_sigma_EBE").c_str(), "@0*(1+@1)"                              , RooArgList(CMS_zz4l_massErr, tth_sigma_CB));
        RooFormulaVar tth_sigma_EBM(("sig_tth_"+chstr+tevstr+"_sigma_EBM").c_str(), "(2.4e-11*@0*@0*@0*@0 + 1.1e-07*@0*@0*@0 - 3.3e-05*@0*@0 + 0.013*@0)", RooArgList(masshiggs));
        RooFormulaVar tth_norm     ("ttH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_tth, tth_xsecbr, yield_var_tth));

        RooFormulaVar *CMS_zz4l_absMassErr;
        if      (ch == 0) CMS_zz4l_absMassErr = new RooFormulaVar("CMS_zz4l_absMassErr","@0*@1*(1+@2)", RooArgList(masshiggs, CMS_zz4l_massErr, sig_sigma_err_4mu));
        else if (ch == 1) CMS_zz4l_absMassErr = new RooFormulaVar("CMS_zz4l_absMassErr","@0*@1*(1+@2)", RooArgList(masshiggs, CMS_zz4l_massErr, sig_sigma_err_4e));
        else              CMS_zz4l_absMassErr = new RooFormulaVar("CMS_zz4l_absMassErr", "@0*@1*TMath::Sqrt((1+@2)*(1+@3))", RooArgList(masshiggs, CMS_zz4l_massErr, sig_sigma_err_4mu, sig_sigma_err_4e));

        RooFormulaVar ggh_ebe_LdM(("sig_ggh_"+chstr+tevstr+"_ebe_LdM").c_str(), getSignalEBELandauMeanString(ch, do7TeV).c_str()     , RooArgList(masshiggs));
        RooFormulaVar ggh_ebe_LdS(("sig_ggh_"+chstr+tevstr+"_ebe_LdS").c_str(), getSignalEBELandauSigmaString(ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar ggh_ebe_LnK(("sig_ggh_"+chstr+tevstr+"_ebe_LnK").c_str(), getSignalEBELogNKString(ch, do7TeV).c_str()          , RooArgList(masshiggs));
        RooFormulaVar ggh_ebe_GaM(("sig_ggh_"+chstr+tevstr+"_ebe_GaM").c_str(), getSignalEBEGaussianMeanString(ch, do7TeV).c_str()   , RooArgList(masshiggs));
        RooFormulaVar ggh_ebe_GaS(("sig_ggh_"+chstr+tevstr+"_ebe_GaS").c_str(), getSignalEBEGaussianSigmaString(ch, do7TeV).c_str()  , RooArgList(masshiggs));
        RooFormulaVar ggh_ebe_frac(("sig_ggh_"+chstr+tevstr+"_ebe_frac").c_str(), getSignalEBELandauFracString(ch, do7TeV).c_str()   , RooArgList(masshiggs));
    
        RooFormulaVar vbf_ebe_LdM(("sig_vbf_"+chstr+tevstr+"_ebe_LdM").c_str(), getSignalEBELandauMeanString(ch, do7TeV).c_str()     , RooArgList(masshiggs));
        RooFormulaVar vbf_ebe_LdS(("sig_vbf_"+chstr+tevstr+"_ebe_LdS").c_str(), getSignalEBELandauSigmaString(ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar vbf_ebe_LnK(("sig_vbf_"+chstr+tevstr+"_ebe_LnK").c_str(), getSignalEBELogNKString(ch, do7TeV).c_str()          , RooArgList(masshiggs));
        RooFormulaVar vbf_ebe_GaM(("sig_vbf_"+chstr+tevstr+"_ebe_GaM").c_str(), getSignalEBEGaussianMeanString(ch, do7TeV).c_str()   , RooArgList(masshiggs));
        RooFormulaVar vbf_ebe_GaS(("sig_vbf_"+chstr+tevstr+"_ebe_GaS").c_str(), getSignalEBEGaussianSigmaString(ch, do7TeV).c_str()  , RooArgList(masshiggs));
        RooFormulaVar vbf_ebe_frac(("sig_vbf_"+chstr+tevstr+"_ebe_frac").c_str(), getSignalEBELandauFracString(ch, do7TeV).c_str()   , RooArgList(masshiggs));

        RooFormulaVar whi_ebe_LdM(("sig_whi_"+chstr+tevstr+"_ebe_LdM").c_str(), getSignalEBELandauMeanString(ch, do7TeV).c_str()     , RooArgList(masshiggs));
        RooFormulaVar whi_ebe_LdS(("sig_whi_"+chstr+tevstr+"_ebe_LdS").c_str(), getSignalEBELandauSigmaString(ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar whi_ebe_LnK(("sig_whi_"+chstr+tevstr+"_ebe_LnK").c_str(), getSignalEBELogNKString(ch, do7TeV).c_str()          , RooArgList(masshiggs));
        RooFormulaVar whi_ebe_GaM(("sig_whi_"+chstr+tevstr+"_ebe_GaM").c_str(), getSignalEBEGaussianMeanString(ch, do7TeV).c_str()   , RooArgList(masshiggs));
        RooFormulaVar whi_ebe_GaS(("sig_whi_"+chstr+tevstr+"_ebe_GaS").c_str(), getSignalEBEGaussianSigmaString(ch, do7TeV).c_str()  , RooArgList(masshiggs));
        RooFormulaVar whi_ebe_frac(("sig_whi_"+chstr+tevstr+"_ebe_frac").c_str(), getSignalEBELandauFracString(ch, do7TeV).c_str()   , RooArgList(masshiggs));

        RooFormulaVar zhi_ebe_LdM(("sig_zhi_"+chstr+tevstr+"_ebe_LdM").c_str(), getSignalEBELandauMeanString(ch, do7TeV).c_str()     , RooArgList(masshiggs));
        RooFormulaVar zhi_ebe_LdS(("sig_zhi_"+chstr+tevstr+"_ebe_LdS").c_str(), getSignalEBELandauSigmaString(ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar zhi_ebe_LnK(("sig_zhi_"+chstr+tevstr+"_ebe_LnK").c_str(), getSignalEBELogNKString(ch, do7TeV).c_str()          , RooArgList(masshiggs));
        RooFormulaVar zhi_ebe_GaM(("sig_zhi_"+chstr+tevstr+"_ebe_GaM").c_str(), getSignalEBEGaussianMeanString(ch, do7TeV).c_str()   , RooArgList(masshiggs));
        RooFormulaVar zhi_ebe_GaS(("sig_zhi_"+chstr+tevstr+"_ebe_GaS").c_str(), getSignalEBEGaussianSigmaString(ch, do7TeV).c_str()  , RooArgList(masshiggs));
        RooFormulaVar zhi_ebe_frac(("sig_zhi_"+chstr+tevstr+"_ebe_frac").c_str(), getSignalEBELandauFracString(ch, do7TeV).c_str()   , RooArgList(masshiggs));

        RooFormulaVar tth_ebe_LdM(("sig_tth_"+chstr+tevstr+"_ebe_LdM").c_str(), getSignalEBELandauMeanString(ch, do7TeV).c_str()     , RooArgList(masshiggs));
        RooFormulaVar tth_ebe_LdS(("sig_tth_"+chstr+tevstr+"_ebe_LdS").c_str(), getSignalEBELandauSigmaString(ch, do7TeV).c_str()    , RooArgList(masshiggs));
        RooFormulaVar tth_ebe_LnK(("sig_tth_"+chstr+tevstr+"_ebe_LnK").c_str(), getSignalEBELogNKString(ch, do7TeV).c_str()          , RooArgList(masshiggs));
        RooFormulaVar tth_ebe_GaM(("sig_tth_"+chstr+tevstr+"_ebe_GaM").c_str(), getSignalEBEGaussianMeanString(ch, do7TeV).c_str()   , RooArgList(masshiggs));
        RooFormulaVar tth_ebe_GaS(("sig_tth_"+chstr+tevstr+"_ebe_GaS").c_str(), getSignalEBEGaussianSigmaString(ch, do7TeV).c_str()  , RooArgList(masshiggs));
        RooFormulaVar tth_ebe_frac(("sig_tth_"+chstr+tevstr+"_ebe_frac").c_str(), getSignalEBELandauFracString(ch, do7TeV).c_str()   , RooArgList(masshiggs));

        RooRealVar ggh_scale_BW (("sig_ggh_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar vbf_scale_BW (("sig_vbf_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar whi_scale_BW (("sig_whi_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar zhi_scale_BW (("sig_zhi_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);
        RooRealVar tth_scale_BW (("sig_tth_"+chstr+tevstr+"_scale_BW" ).c_str(), "BW Gamma Scale", 1.0);

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
      
        qqzz_ebe_LdM .setConstant(kTRUE);
        qqzz_ebe_LdS .setConstant(kTRUE);
        qqzz_ebe_LnK .setConstant(kTRUE);
        qqzz_ebe_GaM .setConstant(kTRUE);
        qqzz_ebe_GaS .setConstant(kTRUE);
        qqzz_ebe_frac.setConstant(kTRUE);

        ggzz_ebe_LdM .setConstant(kTRUE);
        ggzz_ebe_LdS .setConstant(kTRUE);
        ggzz_ebe_LnK .setConstant(kTRUE);
        ggzz_ebe_GaM .setConstant(kTRUE);
        ggzz_ebe_GaS .setConstant(kTRUE);
        ggzz_ebe_frac.setConstant(kTRUE);

        zjets_ebe_LdM .setConstant(kTRUE);
        zjets_ebe_LdS .setConstant(kTRUE);
        zjets_ebe_LnK .setConstant(kTRUE);
        zjets_ebe_GaM .setConstant(kTRUE);
        zjets_ebe_GaS .setConstant(kTRUE);
        zjets_ebe_frac.setConstant(kTRUE);

        masshiggs        .setConstant(kTRUE);
        sig_mean_err_4mu .setConstant(kTRUE);
        sig_mean_err_4e  .setConstant(kTRUE);
        sig_sigma_err_4mu.setConstant(kTRUE);
        sig_sigma_err_4e .setConstant(kTRUE);
        sig_gamma_err    .setConstant(kTRUE);

        ggh_scale_BW.setConstant(kTRUE);
        vbf_scale_BW.setConstant(kTRUE);
        whi_scale_BW.setConstant(kTRUE);
        zhi_scale_BW.setConstant(kTRUE);
        tth_scale_BW.setConstant(kTRUE);
        
        ////////////////// Define the PDFs /////////////////////////////////
        
        const char* bkg_qqzz_pdf_name  = (do1D && !doMassError) ? "bkg_qqzz"  : "bkg_qqzz_1D" ;
        const char* bkg_ggzz_pdf_name  = (do1D && !doMassError) ? "bkg_ggzz"  : "bkg_ggzz_1D" ;
        const char* bkg_zjets_pdf_name = (do1D && !doMassError) ? "bkg_zjets" : "bkg_zjets_1D";
        const char* sig_ggH_pdf_name   = (do1D && !doMassError) ? "ggH"       : "ggH_1D"  ;
        const char* sig_VBF_pdf_name   = (do1D && !doMassError) ? "qqH"       : "qqH_1D"  ;
        const char* sig_WHi_pdf_name   = (do1D && !doMassError) ? "WH"        : "WH_1D"  ;
        const char* sig_ZHi_pdf_name   = (do1D && !doMassError) ? "ZH"        : "ZH_1D"  ;
        const char* sig_ttH_pdf_name   = (do1D && !doMassError) ? "ttH"       : "ttH_1D"  ;
       
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
      

        RooLandau  bkgLD_qqzz_EBE("bkgLD_qqzz_EBE", "", CMS_zz4l_massErr, qqzz_ebe_LdM, qqzz_ebe_LdS);
        RooLandau  bkgLD_ggzz_EBE("bkgLD_ggzz_EBE", "", CMS_zz4l_massErr, ggzz_ebe_LdM, ggzz_ebe_LdS);
        RooLandau  bkgLD_zjets_EBE("bkgLD_zjets_EBE", "", CMS_zz4l_massErr, zjets_ebe_LdM, zjets_ebe_LdS);

        RooLognormal  bkgLN_qqzz_EBE("bkgLN_qqzz_EBE", "", CMS_zz4l_massErr, qqzz_ebe_LdM, qqzz_ebe_LnK);
        RooLognormal  bkgLN_ggzz_EBE("bkgLN_ggzz_EBE", "", CMS_zz4l_massErr, ggzz_ebe_LdM, ggzz_ebe_LnK);
        RooLognormal  bkgLN_zjets_EBE("bkgLN_zjets_EBE", "", CMS_zz4l_massErr, zjets_ebe_LdM, zjets_ebe_LnK);

        RooGaussian  bkgGA_qqzz_EBE("bkgGA_qqzz_EBE", "", CMS_zz4l_massErr, qqzz_ebe_GaM, qqzz_ebe_GaS);
        RooGaussian  bkgGA_ggzz_EBE("bkgGA_ggzz_EBE", "", CMS_zz4l_massErr, ggzz_ebe_GaM, ggzz_ebe_GaS);
        RooGaussian  bkgGA_zjets_EBE("bkgGA_zjets_EBE", "", CMS_zz4l_massErr, zjets_ebe_GaM, zjets_ebe_GaS);        
        
        RooAddPdf *bkg_qqzz_EBE, *bkg_ggzz_EBE, *bkg_zjets_EBE;
        if(ch==0) {
          bkg_qqzz_EBE = new RooAddPdf("bkg_qqzz_EBE", "", bkgLD_qqzz_EBE, bkgLN_qqzz_EBE, qqzz_ebe_frac);
          bkg_ggzz_EBE = new RooAddPdf("bkg_ggzz_EBE", "", bkgLD_ggzz_EBE, bkgLN_ggzz_EBE, ggzz_ebe_frac);
          bkg_zjets_EBE = new RooAddPdf("bkg_zjets_EBE", "", bkgLD_zjets_EBE, bkgLN_zjets_EBE, zjets_ebe_frac);
        } else {
          bkg_qqzz_EBE = new RooAddPdf("bkg_qqzz_EBE", "", bkgLD_qqzz_EBE, bkgGA_qqzz_EBE, qqzz_ebe_frac);
          bkg_ggzz_EBE = new RooAddPdf("bkg_ggzz_EBE", "", bkgLD_ggzz_EBE, bkgGA_ggzz_EBE, ggzz_ebe_frac);
          bkg_zjets_EBE = new RooAddPdf("bkg_zjets_EBE", "", bkgLD_zjets_EBE, bkgGA_zjets_EBE, zjets_ebe_frac);
        }

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

        RooDoubleCB  signalCB_ggH   (doFFT ? "signalCB_ggH" : sig_ggH_pdf_name, "", CMS_zz4l_mass_1D, ggh_mean_CB,ggh_sigma_CB,ggh_alphaL,ggh_nL,ggh_alphaR,ggh_nR);
        RooDoubleCB  signalCB_VBF   (doFFT ? "signalCB_VBF" : sig_VBF_pdf_name, "", CMS_zz4l_mass_1D, vbf_mean_CB,vbf_sigma_CB,vbf_alphaL,vbf_nL,vbf_alphaR,vbf_nR);
        RooDoubleCB  signalCB_WHi   (doFFT ? "signalCB_WHi" : sig_WHi_pdf_name, "", CMS_zz4l_mass_1D, whi_mean_CB,whi_sigma_CB,whi_alphaL,whi_nL,whi_alphaR,whi_nR);
        RooDoubleCB  signalCB_ZHi   (doFFT ? "signalCB_ZHi" : sig_ZHi_pdf_name, "", CMS_zz4l_mass_1D, zhi_mean_CB,zhi_sigma_CB,zhi_alphaL,zhi_nL,zhi_alphaR,zhi_nR);
        RooDoubleCB  signalCB_ttH   (doFFT ? "signalCB_ttH" : sig_ttH_pdf_name, "", CMS_zz4l_mass_1D, tth_mean_CB,tth_sigma_CB,tth_alphaL,tth_nL,tth_alphaR,tth_nR);

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


        RooDoubleCB  signalCB_ggH_EBE("signalCB_ggH_EBE", "", CMS_zz4l_mass_1D, ggh_mean_CB,*CMS_zz4l_absMassErr,ggh_alphaL,ggh_alphaR,ggh_nL,ggh_nR);
        RooDoubleCB  signalCB_VBF_EBE("signalCB_VBF_EBE", "", CMS_zz4l_mass_1D, vbf_mean_CB,*CMS_zz4l_absMassErr,vbf_alphaL,vbf_alphaR,vbf_nL,vbf_nR);
        RooDoubleCB  signalCB_WHi_EBE("signalCB_WHi_EBE", "", CMS_zz4l_mass_1D, whi_mean_CB,*CMS_zz4l_absMassErr,whi_alphaL,whi_alphaR,whi_nL,whi_nR);
        RooDoubleCB  signalCB_ZHi_EBE("signalCB_ZHi_EBE", "", CMS_zz4l_mass_1D, zhi_mean_CB,*CMS_zz4l_absMassErr,zhi_alphaL,zhi_alphaR,zhi_nL,zhi_nR);
        RooDoubleCB  signalCB_ttH_EBE("signalCB_ttH_EBE", "", CMS_zz4l_mass_1D, tth_mean_CB,*CMS_zz4l_absMassErr,tth_alphaL,tth_alphaR,tth_nL,tth_nR);

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


        RooLandau  sigLD_ggH_EBE("sigLD_ggH_EBE", "", CMS_zz4l_massErr, ggh_ebe_LdM, ggh_ebe_LdS);
        RooLandau  sigLD_VBF_EBE("sigLD_VBF_EBE", "", CMS_zz4l_massErr, vbf_ebe_LdM, vbf_ebe_LdS);
        RooLandau  sigLD_WHi_EBE("sigLD_WHi_EBE", "", CMS_zz4l_massErr, whi_ebe_LdM, whi_ebe_LdS);
        RooLandau  sigLD_ZHi_EBE("sigLD_ZHi_EBE", "", CMS_zz4l_massErr, zhi_ebe_LdM, zhi_ebe_LdS);
        RooLandau  sigLD_ttH_EBE("sigLD_ttH_EBE", "", CMS_zz4l_massErr, tth_ebe_LdM, tth_ebe_LdS);

        RooLognormal  sigLN_ggH_EBE("sigLN_ggH_EBE", "", CMS_zz4l_massErr, ggh_ebe_LdM, ggh_ebe_LnK);
        RooLognormal  sigLN_VBF_EBE("sigLN_VBF_EBE", "", CMS_zz4l_massErr, vbf_ebe_LdM, vbf_ebe_LnK);
        RooLognormal  sigLN_WHi_EBE("sigLN_WHi_EBE", "", CMS_zz4l_massErr, whi_ebe_LdM, whi_ebe_LnK);
        RooLognormal  sigLN_ZHi_EBE("sigLN_ZHi_EBE", "", CMS_zz4l_massErr, zhi_ebe_LdM, zhi_ebe_LnK);
        RooLognormal  sigLN_ttH_EBE("sigLN_ttH_EBE", "", CMS_zz4l_massErr, tth_ebe_LdM, tth_ebe_LnK);

        RooGaussian  sigGA_ggH_EBE("sigGA_ggH_EBE", "", CMS_zz4l_massErr, ggh_ebe_GaM, ggh_ebe_GaS);
        RooGaussian  sigGA_VBF_EBE("sigGA_VBF_EBE", "", CMS_zz4l_massErr, vbf_ebe_GaM, vbf_ebe_GaS);
        RooGaussian  sigGA_WHi_EBE("sigGA_WHi_EBE", "", CMS_zz4l_massErr, whi_ebe_GaM, whi_ebe_GaS);
        RooGaussian  sigGA_ZHi_EBE("sigGA_ZHi_EBE", "", CMS_zz4l_massErr, zhi_ebe_GaM, zhi_ebe_GaS);
        RooGaussian  sigGA_ttH_EBE("sigGA_ttH_EBE", "", CMS_zz4l_massErr, tth_ebe_GaM, tth_ebe_GaS);

        RooAddPdf *sig_ggH_EBE, *sig_VBF_EBE, *sig_WHi_EBE, *sig_ZHi_EBE, *sig_ttH_EBE;
        if(ch==0) {
          sig_ggH_EBE = new RooAddPdf("sig_ggH_EBE", "", sigLD_ggH_EBE, sigLN_ggH_EBE, ggh_ebe_frac);
          sig_VBF_EBE = new RooAddPdf("sig_VBF_EBE", "", sigLD_VBF_EBE, sigLN_VBF_EBE, vbf_ebe_frac);
          sig_WHi_EBE = new RooAddPdf("sig_WHi_EBE", "", sigLD_WHi_EBE, sigLN_WHi_EBE, whi_ebe_frac);
          sig_ZHi_EBE = new RooAddPdf("sig_ZHi_EBE", "", sigLD_ZHi_EBE, sigLN_ZHi_EBE, zhi_ebe_frac);
          sig_ttH_EBE = new RooAddPdf("sig_ttH_EBE", "", sigLD_ttH_EBE, sigLN_ttH_EBE, tth_ebe_frac);
        } else {
          sig_ggH_EBE = new RooAddPdf("sig_ggH_EBE", "", sigLD_ggH_EBE, sigGA_ggH_EBE, ggh_ebe_frac);
          sig_VBF_EBE = new RooAddPdf("sig_VBF_EBE", "", sigLD_VBF_EBE, sigGA_VBF_EBE, vbf_ebe_frac);
          sig_WHi_EBE = new RooAddPdf("sig_WHi_EBE", "", sigLD_WHi_EBE, sigGA_WHi_EBE, whi_ebe_frac);
          sig_ZHi_EBE = new RooAddPdf("sig_ZHi_EBE", "", sigLD_ZHi_EBE, sigGA_ZHi_EBE, zhi_ebe_frac);
          sig_ttH_EBE = new RooAddPdf("sig_ttH_EBE", "", sigLD_ttH_EBE, sigGA_ttH_EBE, tth_ebe_frac);
        }


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
            qqzz_fullyield = do7TeV ? 20.7 : 91.6;
            ggzz_fullyield = do7TeV ? 1.29 : 5.00;
            zjet_fullyield = do7TeV ? 1.03 : 3.1;
        }

        else if (ch == 1) {
            qqzz_fullyield = do7TeV ? 13.5 : 59.5;
            ggzz_fullyield = do7TeV ? 0.92 : 3.59;
            zjet_fullyield = do7TeV ? 0.22 : 6.4;
        }

        else if (ch == 2){
            qqzz_fullyield = do7TeV ? 32.5 : 70.0;
            ggzz_fullyield = do7TeV ? 2.19 : 4.99;
            zjet_fullyield = do7TeV ? 0.50 : 1.5;
        }

        else {
            qqzz_fullyield = do7TeV ? 32.5 : 75.6;
            ggzz_fullyield = do7TeV ? 2.19 : 5.54;
            zjet_fullyield = do7TeV ? 0.50 : 5.5;
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
        if (mass<=300) card += createCardTemplate(do7TeV, ch, do1D, workspace.c_str(), true);
        else card += createCardTemplateNoVH(do7TeV, ch, do1D, workspace.c_str(), mass<400. ? true : false);

        std::string binname;
        if (ch == 0) binname = "hzz4mu";
        if (ch == 1) binname = "hzz4e";
        if (ch == 2) binname = "hzz2e2mu";
        if (ch == 3) binname = "hzz2mu2e";

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
        card = findAndReplace(card, "ZX_SYST"        , getZXSystematicsDown(do7TeV, ch), getZXSystematicsUp(do7TeV, ch));
        card = findAndReplace(card, "SIG_GGH_YIELD"  , (getXsecggHByChannel(mass, ch, true) + getXsecVBFByChannel(mass, ch, true) + getXsecZHiByChannel(mass, ch, true) + getXsecWHiByChannel(mass, ch, true) + getXsecttHByChannel(mass, ch, true))/getXsecggHByChannel(mass, ch, true));
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
        if (do1D) {
          if (!doMassError) { 
            if (doFFT) {
              w.import(sig_ggH_pdf);
              w.import(sig_VBF_pdf);
              w.import(sig_WHi_pdf);
              w.import(sig_ZHi_pdf);
              w.import(sig_ttH_pdf);
              w.import(bkg_qqzz_pdf);
              w.import(bkg_ggzz_pdf);
              w.import(bkg_zjets_pdf);
            }        
            else {
              w.import(signalCB_ggH);
              w.import(signalCB_VBF);
              w.import(signalCB_WHi);
              w.import(signalCB_ZHi);
              w.import(signalCB_ttH);
              w.import(bkg_qqzz_pdf);
              w.import(bkg_ggzz_pdf);
              w.import(bkg_zjets_pdf);
            }       
          }
          else {
            if (doFFT) {
              RooProdPdf sig_ggH_pdf_m4l_merr("ggH",  "", sig_ggH_pdf_EBE  ,Conditional(*sig_ggH_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_VBF_pdf_m4l_merr("qqH",  "", sig_VBF_pdf_EBE  ,Conditional(*sig_VBF_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_WHi_pdf_m4l_merr("WH" ,  "", sig_WHi_pdf_EBE  ,Conditional(*sig_WHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ZHi_pdf_m4l_merr("ZH" ,  "", sig_ZHi_pdf_EBE  ,Conditional(*sig_ZHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ttH_pdf_m4l_merr("ttH",  "", sig_ttH_pdf_EBE  ,Conditional(*sig_ttH_EBE, RooArgSet(CMS_zz4l_massErr)));
              
              RooProdPdf bkg_qqzz_pdf_m4l_merr("bkg_qqzz", "", bkg_qqzz_pdf, Conditional(*bkg_qqzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_ggzz_pdf_m4l_merr("bkg_ggzz", "", bkg_ggzz_pdf, Conditional(*bkg_ggzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_zjets_pdf_m4l_merr("bkg_zjets", "", bkg_zjets_pdf, Conditional(*bkg_zjets_EBE, RooArgSet(CMS_zz4l_massErr)));
              
              w.import(bkg_qqzz_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_ggzz_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_zjets_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ggH_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_VBF_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_WHi_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ZHi_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ttH_pdf_m4l_merr,RooFit::RecycleConflictNodes());
            }
            
            else {
              RooProdPdf sig_ggH_pdf_m4l_merr("ggH",  "", signalCB_ggH_EBE  ,Conditional(*sig_ggH_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_VBF_pdf_m4l_merr("qqH",  "", signalCB_VBF_EBE  ,Conditional(*sig_VBF_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_WHi_pdf_m4l_merr("WH" ,  "", signalCB_WHi_EBE  ,Conditional(*sig_WHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ZHi_pdf_m4l_merr("ZH" ,  "", signalCB_ZHi_EBE  ,Conditional(*sig_ZHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ttH_pdf_m4l_merr("ttH",  "", signalCB_ttH_EBE  ,Conditional(*sig_ttH_EBE, RooArgSet(CMS_zz4l_massErr)));                
              
              RooProdPdf bkg_qqzz_pdf_m4l_merr("bkg_qqzz", "", bkg_qqzz_pdf, Conditional(*bkg_qqzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_ggzz_pdf_m4l_merr("bkg_ggzz", "", bkg_ggzz_pdf, Conditional(*bkg_ggzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_zjets_pdf_m4l_merr("bkg_zjets", "", bkg_zjets_pdf, Conditional(*bkg_zjets_EBE, RooArgSet(CMS_zz4l_massErr)));
              
              w.import(bkg_qqzz_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_ggzz_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_zjets_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ggH_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_VBF_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_WHi_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ZHi_pdf_m4l_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ttH_pdf_m4l_merr,RooFit::RecycleConflictNodes());
            }
          }
        }

        else {
          TFile melafile(melafilename.c_str());
          RooArgList v2dList(CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
          RooArgSet  v2dSet (CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
          
          std::string chstrsmall;
          if (ch == 0) chstrsmall = "mm";
          if (ch == 1) chstrsmall = "ee";
          if (ch == 2) chstrsmall = "em";
          if (ch == 3) chstrsmall = "me";
          
          TH2F* melashape_qz = (TH2F*)(melafile.Get(("hist2D_bkg_"+chstrsmall).c_str()));
          TH2F* melashape_gz = (TH2F*)(melafile.Get(("hist2D_b2g_"+chstrsmall).c_str()));
          TH2F* melashape_si = (TH2F*)(melafile.Get(("hist2D_sig_"+chstrsmall).c_str()));
          TH2F* melashape_zu = (TH2F*)(melafile.Get(("hist2D_bkg_"+chstrsmall+"_up").c_str()));
          TH2F* melashape_zd = (TH2F*)(melafile.Get(("hist2D_bkg_"+chstrsmall+"_dn").c_str()));
          
          RooDataHist rhist_qqzz (("rhist_qqzz_" +chstr+tevstr).c_str(), "", v2dList, melashape_qz);
          RooDataHist rhist_ggzz (("rhist_ggzz_" +chstr+tevstr).c_str(), "", v2dList, melashape_gz);
          RooDataHist rhist_zjets(("rhist_zjets_"+chstr+tevstr).c_str(), "", v2dList, melashape_qz);
          RooDataHist rhist_ggH  (("rhist_ggH_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
          RooDataHist rhist_VBF  (("rhist_VBF_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
          RooDataHist rhist_WHi  (("rhist_WHi_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
          RooDataHist rhist_ZHi  (("rhist_ZHi_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
          RooDataHist rhist_ttH  (("rhist_ttH_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
          RooDataHist rhist_zjup (("rhist_zjup_" +chstr+tevstr).c_str(), "", v2dList, melashape_zu);
          RooDataHist rhist_zjdn (("rhist_zjdn_" +chstr+tevstr).c_str(), "", v2dList, melashape_zd);
          RooHistPdf rpdf_qqzz (("bkg_qqzz_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_qqzz);
          RooHistPdf rpdf_ggzz (("bkg_ggzz_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_ggzz);
          RooHistPdf rpdf_zjets(("bkg_zjets_mela2D_pdf_"+chstr+tevstr).c_str(), "", v2dSet , rhist_zjets);
          RooHistPdf rpdf_ggH  (("sig_ggH_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_ggH);
          RooHistPdf rpdf_VBF  (("sig_VBF_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_VBF);
          RooHistPdf rpdf_WHi  (("sig_WHi_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_WHi);
          RooHistPdf rpdf_ZHi  (("sig_ZHi_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_ZHi);
          RooHistPdf rpdf_ttH  (("sig_ttH_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_ttH);
          
          RooHistPdf rpdf_zjup (("bkg_zjup_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_zjup);
          RooHistPdf rpdf_zjdn (("bkg_zjdn_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_zjdn);
          
          RooArgList list_zjets;
          list_zjets.add(rpdf_zjets);
          list_zjets.add(rpdf_zjup);
          list_zjets.add(rpdf_zjdn);
          
          RooRealVar CMS_zz4l_bkgMELA("CMS_zz4l_bkgMELA" ,"" ,0,-10,10); 
          
          FastVerticalInterpHistPdf2D plpdf_qqzz (("bkg_qqzz_FVIHP_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_qqzz) ,RooArgList()                ,1.0,1);
          FastVerticalInterpHistPdf2D plpdf_ggzz (("bkg_ggzz_FVIHP_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggzz) ,RooArgList()                ,1.0,1);
          FastVerticalInterpHistPdf2D plpdf_zjets(("bkg_zjets_FVIHP_"+chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,list_zjets            ,RooArgList(CMS_zz4l_bkgMELA),1.0,1);
          FastVerticalInterpHistPdf2D plpdf_ggH  (("sig_ggH_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggH)  ,RooArgList()                ,1.0,1);
          FastVerticalInterpHistPdf2D plpdf_VBF  (("sig_VBF_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_VBF)  ,RooArgList()                ,1.0,1);
          FastVerticalInterpHistPdf2D plpdf_WHi  (("sig_WHi_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_WHi)  ,RooArgList()                ,1.0,1);
          FastVerticalInterpHistPdf2D plpdf_ZHi  (("sig_ZHi_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ZHi)  ,RooArgList()                ,1.0,1);
          FastVerticalInterpHistPdf2D plpdf_ttH  (("sig_ttH_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ttH)  ,RooArgList()                ,1.0,1);
          
          if (doFFT) { 
            if (doMassError) {
              RooProdPdf sig_ggH_pdf_m4l_merr("ggH_m4l_merr",  "", sig_ggH_pdf_EBE  ,Conditional(*sig_ggH_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_VBF_pdf_m4l_merr("VBF_m4l_merr",  "", sig_VBF_pdf_EBE  ,Conditional(*sig_VBF_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_WHi_pdf_m4l_merr("WHi_m4l_merr" ,  "", sig_WHi_pdf_EBE  ,Conditional(*sig_WHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ZHi_pdf_m4l_merr("ZHi_m4l_merr" ,  "", sig_ZHi_pdf_EBE  ,Conditional(*sig_ZHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ttH_pdf_m4l_merr("ttH_m4l_merr",  "", sig_ttH_pdf_EBE  ,Conditional(*sig_ttH_EBE, RooArgSet(CMS_zz4l_massErr)));
              
              RooProdPdf bkg_qqzz_pdf_m4l_merr("bkg_qqzz_m4l_merr", "", bkg_qqzz_pdf, Conditional(*bkg_qqzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_ggzz_pdf_m4l_merr("bkg_ggzz_m4l_merr", "", bkg_ggzz_pdf, Conditional(*bkg_ggzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_zjets_pdf_m4l_merr("bkg_zjets_m4l_merr", "", bkg_zjets_pdf, Conditional(*bkg_zjets_EBE, RooArgSet(CMS_zz4l_massErr)));
              
              
              RooProdPdf sig_ggH_pdf_2D_merr("ggH",  "", sig_ggH_pdf_m4l_merr, Conditional(plpdf_ggH, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_VBF_pdf_2D_merr("qqH",  "", sig_VBF_pdf_m4l_merr, Conditional(plpdf_VBF, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_WHi_pdf_2D_merr("WH" ,  "", sig_WHi_pdf_m4l_merr, Conditional(plpdf_WHi, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_ZHi_pdf_2D_merr("ZH" ,  "", sig_ZHi_pdf_m4l_merr, Conditional(plpdf_ZHi, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_ttH_pdf_2D_merr("ttH",  "", sig_ttH_pdf_m4l_merr, Conditional(plpdf_ttH, RooArgSet(CMS_zz4l_melaLD)));
              
              RooProdPdf bkg_qqzz_pdf_2D_merr ("bkg_qqzz" , "", bkg_qqzz_pdf_m4l_merr  ,Conditional(plpdf_qqzz, RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf bkg_ggzz_pdf_2D_merr ("bkg_ggzz" , "", bkg_ggzz_pdf_m4l_merr  ,Conditional(plpdf_ggzz, RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf bkg_zjets_pdf_2D_merr("bkg_zjets", "", bkg_zjets_pdf_m4l_merr ,Conditional(plpdf_zjets,RooArgSet(CMS_zz4l_melaLD))); 
                
              w.import(bkg_qqzz_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_ggzz_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_zjets_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ggH_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_VBF_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_WHi_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ZHi_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ttH_pdf_2D_merr,RooFit::RecycleConflictNodes());
            }

            else {
              RooProdPdf sig_ggH_pdf_2D  ("ggH"    ,  "", sig_ggH_pdf   ,Conditional(plpdf_ggH  , RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf sig_VBF_pdf_2D  ("qqH"    ,  "", sig_VBF_pdf   ,Conditional(plpdf_VBF  , RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf sig_WHi_pdf_2D  ("WH"     ,  "", sig_WHi_pdf   ,Conditional(plpdf_WHi  , RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf sig_ZHi_pdf_2D  ("ZH"     ,  "", sig_ZHi_pdf   ,Conditional(plpdf_ZHi  , RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf sig_ttH_pdf_2D  ("ttH"    ,  "", sig_ttH_pdf   ,Conditional(plpdf_ttH  , RooArgSet(CMS_zz4l_melaLD))); 
              
              RooProdPdf bkg_qqzz_pdf_2D ("bkg_qqzz", "", bkg_qqzz_pdf  ,Conditional(plpdf_qqzz , RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf bkg_ggzz_pdf_2D ("bkg_ggzz", "", bkg_ggzz_pdf  ,Conditional(plpdf_ggzz , RooArgSet(CMS_zz4l_melaLD))); 
              RooProdPdf bkg_zjets_pdf_2D("bkg_zjets","", bkg_zjets_pdf ,Conditional(plpdf_zjets, RooArgSet(CMS_zz4l_melaLD))); 
              
              w.import(bkg_qqzz_pdf_2D); 
              w.import(bkg_ggzz_pdf_2D); 
              w.import(bkg_zjets_pdf_2D); 
              w.import(sig_ggH_pdf_2D); 
              w.import(sig_VBF_pdf_2D);
              w.import(sig_WHi_pdf_2D); 
              w.import(sig_ZHi_pdf_2D);
              w.import(sig_ttH_pdf_2D);
            }
            
          }
          else {
            if (doMassError) {
              RooProdPdf sig_ggH_pdf_m4l_merr("ggH_m4l_merr",  "", signalCB_ggH_EBE  ,Conditional(*sig_ggH_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_VBF_pdf_m4l_merr("VBF_m4l_merr",  "", signalCB_VBF_EBE  ,Conditional(*sig_VBF_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_WHi_pdf_m4l_merr("WHi_m4l_merr",  "", signalCB_WHi_EBE  ,Conditional(*sig_WHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ZHi_pdf_m4l_merr("ZHi_m4l_merr",  "", signalCB_ZHi_EBE  ,Conditional(*sig_ZHi_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf sig_ttH_pdf_m4l_merr("ttH_m4l_merr",  "", signalCB_ttH_EBE  ,Conditional(*sig_ttH_EBE, RooArgSet(CMS_zz4l_massErr)));
              
              RooProdPdf bkg_qqzz_pdf_m4l_merr("bkg_qqzz_m4l_merr", "", bkg_qqzz_pdf, Conditional(*bkg_qqzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_ggzz_pdf_m4l_merr("bkg_ggzz_m4l_merr", "", bkg_ggzz_pdf, Conditional(*bkg_ggzz_EBE, RooArgSet(CMS_zz4l_massErr)));
              RooProdPdf bkg_zjets_pdf_m4l_merr("bkg_zjets_m4l_merr", "", bkg_zjets_pdf,Conditional(*bkg_zjets_EBE, RooArgSet(CMS_zz4l_massErr)));
              
              RooProdPdf sig_ggH_pdf_2D_merr("ggH",  "", sig_ggH_pdf_m4l_merr, Conditional(plpdf_ggH, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_VBF_pdf_2D_merr("qqH",  "", sig_VBF_pdf_m4l_merr, Conditional(plpdf_VBF, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_WHi_pdf_2D_merr("WH" ,  "", sig_WHi_pdf_m4l_merr, Conditional(plpdf_WHi, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_ZHi_pdf_2D_merr("ZH" ,  "", sig_ZHi_pdf_m4l_merr, Conditional(plpdf_ZHi, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf sig_ttH_pdf_2D_merr("ttH",  "", sig_ttH_pdf_m4l_merr, Conditional(plpdf_ttH, RooArgSet(CMS_zz4l_melaLD)));
              
              RooProdPdf bkg_qqzz_pdf_2D_merr ("bkg_qqzz" , "", bkg_qqzz_pdf_m4l_merr  ,Conditional(plpdf_qqzz, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf bkg_ggzz_pdf_2D_merr ("bkg_ggzz" , "", bkg_ggzz_pdf_m4l_merr  ,Conditional(plpdf_ggzz, RooArgSet(CMS_zz4l_melaLD)));
              RooProdPdf bkg_zjets_pdf_2D_merr("bkg_zjets", "", bkg_zjets_pdf_m4l_merr  ,Conditional(plpdf_zjets,RooArgSet(CMS_zz4l_melaLD)));
              
              w.import(bkg_qqzz_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_ggzz_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(bkg_zjets_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ggH_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_VBF_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_WHi_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ZHi_pdf_2D_merr,RooFit::RecycleConflictNodes());
              w.import(sig_ttH_pdf_2D_merr,RooFit::RecycleConflictNodes());
              
            }
            else {
                    RooProdPdf sig_ggH_pdf_2D  ("ggH",      "", signalCB_ggH  ,Conditional(plpdf_ggH  , RooArgSet(CMS_zz4l_melaLD))); 
                    RooProdPdf sig_VBF_pdf_2D  ("qqH",      "", signalCB_VBF  ,Conditional(plpdf_VBF  , RooArgSet(CMS_zz4l_melaLD))); 
                    RooProdPdf sig_WHi_pdf_2D  ("WH" ,      "", signalCB_WHi  ,Conditional(plpdf_WHi  , RooArgSet(CMS_zz4l_melaLD))); 
                    RooProdPdf sig_ZHi_pdf_2D  ("ZH" ,      "", signalCB_ZHi  ,Conditional(plpdf_ZHi  , RooArgSet(CMS_zz4l_melaLD))); 
                    RooProdPdf sig_ttH_pdf_2D  ("ttH",      "", signalCB_ttH  ,Conditional(plpdf_ttH  , RooArgSet(CMS_zz4l_melaLD))); 

                    RooProdPdf bkg_qqzz_pdf_2D ("bkg_qqzz", "", bkg_qqzz_pdf  ,Conditional(plpdf_qqzz , RooArgSet(CMS_zz4l_melaLD))); 
                    RooProdPdf bkg_ggzz_pdf_2D ("bkg_ggzz", "", bkg_ggzz_pdf  ,Conditional(plpdf_ggzz , RooArgSet(CMS_zz4l_melaLD))); 
                    RooProdPdf bkg_zjets_pdf_2D("bkg_zjets","", bkg_zjets_pdf ,Conditional(plpdf_zjets, RooArgSet(CMS_zz4l_melaLD))); 
           
                    w.import(bkg_qqzz_pdf_2D); 
                    w.import(bkg_ggzz_pdf_2D); 
                    w.import(bkg_zjets_pdf_2D); 
                    w.import(sig_ggH_pdf_2D); 
                    w.import(sig_VBF_pdf_2D);
                    w.import(sig_WHi_pdf_2D); 
                    w.import(sig_ZHi_pdf_2D);
                    w.import(sig_ttH_pdf_2D);
                }
            }
        }
        
        w.writeToFile(workspace.c_str());

        delete histxsecbrggh;
        delete histxsecbrvbf;
        delete histxsecbrwhi;
        delete histxsecbrzhi;
        delete histxsecbrtth;

        delete sig_mean_err_al;
        delete sig_sigma_err_al;
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
        for (float i = 110.; i <= 160.; i += 1.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 3);
        }
        for (float i = 162.; i <= 290.; i += 2.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 3);
        }
        for (float i = 295.; i <= 350.; i += 5.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 3);
        }
        for (float i = 360.; i <= 400.; i += 10.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 3);
        }
        for (float i = 420.; i <= 1000.; i += 20.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 3);
        }
    }

    void makeCard(float i) {
        createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
        createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
        createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        createCard(i, getMassCut(i, true), getMassCut(i, false), 3);
    }


};

void doHZZAnalysis() {

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
    hmpi7.doFFT = false;
    hmpi7.doMassError = false;
    if(hmpi7.doMassError) hmpi7.doFFT = false;
    hmpi7.treeFolder = "/home/avartak/CMS/Higgs/Moriond/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";
    hmpi7.melafilename = "mela2DShapes.root";

    init(hmpi7.do7TeV);

    hmpi7.ymaker_data.fill(hmpi7.treeFolder+"data.root");

    hmpi7.makeCards();
    hmpi7.do1D = false;
    hmpi7.makeCards();

    //hmpi7.makeCard(126.);


    HiggsMassPointInfo hmpi8;
    hmpi8.lumi = 19.6; 
    hmpi8.base_lumi = 19.6; 
    hmpi8.z1min = 40.;
    hmpi8.z2min = 12.;
    hmpi8.massLowBkgFit = 100.;
    hmpi8.massHighBkgFit = 1600.;
    hmpi8.melacut = -1.0;
    hmpi8.do1D = true;
    hmpi8.do7TeV = false;
    hmpi8.doFFT = false;
    hmpi8.doMassError = false;
    if(hmpi8.doMassError) hmpi8.doFFT = false;
    hmpi8.treeFolder = "/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/";
    hmpi8.melafilename = "mela2DShapes.root";

    init(hmpi8.do7TeV);


    hmpi8.ymaker_data.fill(hmpi8.treeFolder+"data.root");

    hmpi8.makeCards();
    hmpi8.do1D = false;
    hmpi8.makeCards();

    //hmpi8.makeCard(126.);

}


