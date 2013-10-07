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
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
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
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>

#include "scales2.h"
#include "FakeRateCalculator.h"
#include "CardTemplate.h"
#include "FitMaker.h"
#include "YieldMaker.h"
#include "findAndReplace.h"

using namespace std;
using namespace RooFit;

double wavg(double x, double y, double ex, double ey) {
  return (x/ex/ex + y/ey/ey)/(1/ex/ex + 1/ey/ey);
}

pair<double,double> werror(double x, double y, double ex, double ey) {
  double max = TMath::Max(x+ex,y+ey);
  double min = TMath::Min(x-ex,y-ey);

  double val = wavg(x,y,ex,ey);
  double eu=fabs(max-val);
  double ed=fabs(val-min);
  
  return std::make_pair<double,double>(ed,eu);
}

struct QQZZ4muCardInfo {

  float mean_lan, sigma_lan, k_lnN, f1;

    QQZZ4muCardInfo() {
      mean_lan  = 0.0;
      sigma_lan = 0.0;
      k_lnN     = 0.0;
      f1        = 0.0;
    }


    QQZZ4muCardInfo(ZZEbE4muFitMaker fm) {
      mean_lan  = fm.getVarLanMean();
      sigma_lan = fm.getVarLanSigma();
      k_lnN     = fm.getVarLnNK();
      f1        = fm.getVarF1();
    }


    std::string printInfo() {
        
        stringstream info;
        info << "\t\t\t\tmLan = " << mean_lan  << ";" << std::endl;
        info << "\t\t\t\tsLan = " << sigma_lan << ";" << std::endl;
        info << "\t\t\t\tkLnN = " << k_lnN     << ";" << std::endl;
        info << "\t\t\t\tf1   = " << f1        << ";" << std::endl;
        info << "\n\n";
        return info.str();
    }
};

struct QQZZ4eAnd2e2muCardInfo {

  float mean_lan, sigma_lan, mean_gaus, sigma_gaus, f1;

    QQZZ4eAnd2e2muCardInfo() {
      mean_lan   = 0.0;
      sigma_lan  = 0.0;
      mean_gaus  = 0.0;
      sigma_gaus = 0.0;
      f1         = 0.0;
    }


    QQZZ4eAnd2e2muCardInfo(ZZEbE4eAnd2e2muFitMaker fm) {
      mean_lan   = fm.getVarLanMean();
      sigma_lan  = fm.getVarLanSigma();
      mean_gaus  = fm.getVarGauMean();
      sigma_gaus = fm.getVarGauSigma();
      f1         = fm.getVarF1();
    }


    std::string printInfo() {
        
        stringstream info;
        info << "\t\t\t\tmLan = " << mean_lan   << ";" << std::endl;
        info << "\t\t\t\tsLan = " << sigma_lan  << ";" << std::endl;
        info << "\t\t\t\tmLnN = " << mean_gaus  << ";" << std::endl;
        info << "\t\t\t\tkLnN = " << sigma_gaus << ";" << std::endl;
        info << "\t\t\t\tf1   = " << f1         << ";" << std::endl;
        info << "\n\n";
        return info.str();
    }
};

struct HiggsMassPointInfo {

    float z1min;
    float z2min;
    float massLow;
    float massHigh;
    float massErrLow;
    float massErrHigh;
    float melacut;
    bool  do7TeV;
    int   njetscut;
    bool  njetveto;
    float ptcut;
    float ptmax;

    ZXYieldMaker   ymaker_zxss;
    ZXYieldMaker   ymaker_zxos;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;

  void analyze(int ch, ofstream& file, double fitValues[5], double fitErrors[5], bool issignal=false, int mass=115) {

        float yield_qq  = ymaker_qqzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
        float yield_gg = 0;
        float yield_zjss = 0;
        if(!issignal) {
          yield_gg  = ymaker_ggzz.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
          yield_zjss  = ymaker_zxss.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
          //        float yield_zjos  = ymaker_zxos.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
        }

        float yielde_qq  = ymaker_qqzz.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut);
        float yielde_gg  = 0;
        float yielde_zjss = 0;
        if(!issignal) {
          yielde_gg  = ymaker_ggzz.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut);
          yielde_zjss  = ymaker_zxss.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut);
          //        float yielde_zjos  = ymaker_zxos.getYieldStatError(ch, z1min, z2min, massLow, massHigh, melacut);
        }
        
        if (ch == 2) {
          yield_qq += ymaker_qqzz.getYield(3, z1min, z2min, massLow, massHigh, melacut);
          if(!issignal) {
            yield_gg += ymaker_ggzz.getYield(3, z1min, z2min, massLow, massHigh, melacut);
            yield_zjss += ymaker_zxss.getYield(3, z1min, z2min, massLow, massHigh, melacut);
            //            yield_zjos += ymaker_zxos.getYield(3, z1min, z2min, massLow, massHigh, melacut);
          }
          
          yielde_qq = sqrt(yielde_qq + pow(ymaker_qqzz.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut),2));
          if(!issignal) {
            yielde_gg = sqrt(yielde_gg + pow(ymaker_ggzz.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut),2));
            yielde_zjss = sqrt(yielde_zjss + pow(ymaker_zxss.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut),2));
            //            yielde_zjos = sqrt(yielde_zjos + pow(ymaker_zxos.getYieldStatError(3, z1min, z2min, massLow, massHigh, melacut),2));
          }
        }
        
        std::string chstr;
        
        if(ch==0) {
          ZZEbE4muFitMaker  fitmaker_qz ("bkg_qqzz" , massLow, massHigh);
          ZZEbE4muFitMaker  fitmaker_gz ("bkg_ggzz" , massLow, massHigh);
          ZZEbE4muFitMaker  fitmaker_zxss ("bkg_zjets_ss", 120, 130);

          fitmaker_qz.add(ymaker_qqzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
          if(!issignal) {
            fitmaker_gz.add(ymaker_ggzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
            fitmaker_zxss.add(ymaker_zxss.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
          }

          chstr = "4mu";
          
          std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
          stringstream massstr;
          if(issignal) massstr << mass << "_";

          fitmaker_qz.fit();
          fitmaker_qz.print("cat0/qqZZ_"+chstr+massstr.str()+tevstr, 200);
          fitmaker_qz.print("cat0/qqZZ_"+chstr+massstr.str()+"_zoom"+tevstr, 100, true);
          QQZZ4muCardInfo qzci(fitmaker_qz);
          if(!issignal) {
            fitmaker_gz.fit();
            fitmaker_gz.print("cat0/ggZZ_"+chstr+tevstr, 50, true);
            fitmaker_gz.print("cat0/ggZZ_"+chstr+"_zoom"+tevstr, 50, true);
            
            fitmaker_zxss.fit();
            fitmaker_zxss.print("cat0/zjets_ss_"+chstr+tevstr, 25);
            fitmaker_zxss.print("cat0/zjets_ss_"+chstr+"_zoom"+tevstr, 25);
          } 
          file << "\n\n\n***********************************************************";
          file << "Channel = " << chstr << "********************************************************************\n\n";
          file << qzci.printInfo();
          if(!issignal) {
            QQZZ4muCardInfo gzci(fitmaker_gz);          
            QQZZ4muCardInfo zci(fitmaker_zxss);
            file << gzci.printInfo();
            file << zci.printInfo() ;
          } else {
            fitValues[0] = fitmaker_qz.getVarLanMean();
            fitValues[1] = fitmaker_qz.getVarLanSigma();
            fitValues[2] = fitmaker_qz.getVarLanMean(); // the same as the ladau for 4mu
            fitValues[3] = fitmaker_qz.getVarLnNK();
            fitValues[4] = fitmaker_qz.getVarF1();

            fitErrors[0] = fitmaker_qz.getErrLanMean();
            fitErrors[1] = fitmaker_qz.getErrLanSigma();
            fitErrors[2] = fitmaker_qz.getErrLanMean(); // the same as the ladau for 4mu
            fitErrors[3] = fitmaker_qz.getErrLnNK();
            fitErrors[4] = fitmaker_qz.getErrF1();
          }
        } else {
          if (ch == 1) chstr = "4e";
          if (ch == 2) chstr = "2e2mu";
          
          ZZEbE4eAnd2e2muFitMaker  fitmaker_qz ("bkg_qqzz" , massLow, massHigh);
          ZZEbE4eAnd2e2muFitMaker  fitmaker_gz ("bkg_ggzz" , massLow, massHigh);
          ZZEbE4eAnd2e2muFitMaker  fitmaker_zxss ("bkg_zjets_ss", 120, 130);
          
          //           float f1val = (ch == 1) ? 0.50 : 0.65;
          //           fitmaker_gz.setVarF1(f1val);
          //           fitmaker_qz.setVarF1(f1val);
          
          fitmaker_qz.add(ymaker_qqzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
          if(!issignal) {
            fitmaker_gz.add(ymaker_ggzz.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
            fitmaker_zxss.add(ymaker_zxss.getFitDataSet(ch, z1min, z2min, massLow, massHigh, melacut));
          }

          if (ch == 2) {
            fitmaker_qz.add(ymaker_qqzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut));
            if(!issignal) {
              fitmaker_gz.add(ymaker_ggzz.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut));
              fitmaker_zxss.add(ymaker_zxss.getFitDataSet(3, z1min, z2min, massLow, massHigh, melacut));
            }
          } 
          
          std::string tevstr = do7TeV ? "_7TeV" : "_8TeV";           
          stringstream massstr;
          if(issignal) massstr << "_mH" << mass;

          fitmaker_qz.fit();
          fitmaker_qz.print("cat0/qqZZ_"+chstr+massstr.str()+tevstr, 200);
          fitmaker_qz.print("cat0/qqZZ_"+chstr+massstr.str()+"_zoom"+tevstr, 100, true);
          QQZZ4eAnd2e2muCardInfo qzci(fitmaker_qz);
          if(!issignal) {
            fitmaker_gz.fit();
            fitmaker_gz.print("cat0/ggZZ_"+chstr+tevstr, 50, true);
            fitmaker_gz.print("cat0/ggZZ_"+chstr+"_zoom"+tevstr, 50, true);
            
            fitmaker_zxss.fit();
            fitmaker_zxss.print("cat0/zjets_ss_"+chstr+tevstr, 25);
            fitmaker_zxss.print("cat0/zjets_ss_"+chstr+"_zoom"+tevstr, 25, true);
          }

          file << "\n\n\n***********************************************************";
          file << "Channel = " << chstr << "********************************************************************\n\n";
          file << qzci.printInfo();
          if(!issignal) {
            QQZZ4eAnd2e2muCardInfo gzci(fitmaker_gz);          
            QQZZ4eAnd2e2muCardInfo zci(fitmaker_zxss);
            file << gzci.printInfo();
            file << zci.printInfo() ;
          } else {
            fitValues[0] = fitmaker_qz.getVarLanMean();
            fitValues[1] = fitmaker_qz.getVarLanSigma();
            fitValues[2] = fitmaker_qz.getVarGauMean();
            fitValues[3] = fitmaker_qz.getVarGauSigma();
            fitValues[4] = fitmaker_qz.getVarF1();

            fitErrors[0] = fitmaker_qz.getErrLanMean();
            fitErrors[1] = fitmaker_qz.getErrLanSigma();
            fitErrors[2] = fitmaker_qz.getErrGauMean();
            fitErrors[3] = fitmaker_qz.getErrGauSigma();
            fitErrors[4] = fitmaker_qz.getErrF1();
          }
        }
      
    } 


};

void getParamsBkg() {

    float lumi7     = 5.05;
    float lumi8     = 19.6;

    std::string treeFolder7 = "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_42X_S1_V13_S2_V03/";
    std::string treeFolder8 = "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_53X_S1_V13_S2_V03/";

    enum zxalgo { ss=0, os=1};
    float Z1min[2] = { 40,   81.186 };
    float Z1max[2] = { 120, 101.186 };

    double fitValues[5];
    double fitErrors[5];

    /*
    HiggsMassPointInfo hmpi7;
    hmpi7.z1min    = 40.;
    hmpi7.z2min    = 12.;
    hmpi7.massLow  = 100.;
    hmpi7.massHigh = 160;
    hmpi7.melacut  = -1.;
    hmpi7.do7TeV   = true;
    hmpi7.njetscut = -1;
    hmpi7.ptcut    = 0;
    hmpi7.ptmax    = 1000000;
    hmpi7.njetveto = true;

    init(hmpi7.do7TeV);

    FakeRateCalculator FR_7TeV_SS(treeFolder7+"DATA/7TeV/data2011.root", hmpi7.do7TeV, Z1min[ss], Z1max[ss], 0.0, 0.0, true);
    FakeRateCalculator FR_7TeV_OS(treeFolder7+"DATA/7TeV/data2011.root", hmpi7.do7TeV, Z1min[os], Z1max[os], 0.0, 0.0, true);
    
    hmpi7.ymaker_zxss.fill(treeFolder7+"DATA/7TeV/data2011.root"          , 1.0, FR_7TeV_SS, true);
    //    hmpi7.ymaker_zxos.fill(treeFolder7+"DATA/7TeV/data2011.root"          , 1.0, FR_7TeV_OS, false);

    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id102.root" , getBkgXsec(102)*lumi7/evt_7TeV(102), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id103.root" , getBkgXsec(103)*lumi7/evt_7TeV(103), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id104.root" , getBkgXsec(104)*lumi7/evt_7TeV(104), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id105.root" , getBkgXsec(105)*lumi7/evt_7TeV(105), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id106.root" , getBkgXsec(106)*lumi7/evt_7TeV(106), 0.0, false);
    hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id107.root" , getBkgXsec(107)*lumi7/evt_7TeV(107), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"MC/7TeV/hzzTree_id101.root" , getBkgXsec(101)*lumi7/evt_7TeV(101), 0.0, false);
    hmpi7.ymaker_ggzz.fill(treeFolder7+"MC/7TeV/hzzTree_id100.root" , getBkgXsec(100)*lumi7/evt_7TeV(100), 0.0, false);

    ofstream file7;
    file7.open("cat0/info7TeV.txt");

    hmpi7.analyze(0, file7, fitValues, fitErrors);
    hmpi7.analyze(1, file7, fitValues, fitErrors);
    hmpi7.analyze(2, file7, fitValues, fitErrors);
    
    file7.close();
    */

    HiggsMassPointInfo hmpi8;
    hmpi8.z1min    = 40.;
    hmpi8.z2min    = 12.;
    hmpi8.massLow  = 100;
    hmpi8.massHigh = 160;
    hmpi8.do7TeV   = false;
    hmpi8.melacut  = -1.;
    hmpi8.njetscut = -1;
    hmpi8.ptcut    = 0;
    hmpi8.ptmax    = 1000000.;
    hmpi8.njetveto = true;


    setmc(false);    
    init(hmpi8.do7TeV);

    FakeRateCalculator FR_8TeV_SS(treeFolder8+"DATA/8TeV/data2012.root", hmpi8.do7TeV, Z1min[ss], Z1max[ss], 0.0, 0.0, true);
    FakeRateCalculator FR_8TeV_OS(treeFolder8+"DATA/8TeV/data2012.root", hmpi8.do7TeV, Z1min[os], Z1max[os], 0.0, 0.0, true);
    
    hmpi8.ymaker_zxss.fill(treeFolder8+"DATA/8TeV/data2012.root"          , 1.0, FR_8TeV_SS, true);
    hmpi8.ymaker_zxos.fill(treeFolder8+"DATA/8TeV/data2012.root"          , 1.0, FR_8TeV_OS, false);

    setmc(true);    
    init(hmpi8.do7TeV);

    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id102.root" , getBkgXsec(102)*lumi8/evt_8TeV(102), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id103.root" , getBkgXsec(103)*lumi8/evt_8TeV(103), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id104.root" , getBkgXsec(104)*lumi8/evt_8TeV(104), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id105.root" , getBkgXsec(105)*lumi8/evt_8TeV(105), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id106.root" , getBkgXsec(106)*lumi8/evt_8TeV(106), 0.0, false);
    hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id107.root" , getBkgXsec(107)*lumi8/evt_8TeV(107), 0.0, false);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"MC/8TeV/hzzTree_id101.root" , getBkgXsec(101)*lumi8/evt_8TeV(101), 0.0, false);
    hmpi8.ymaker_ggzz.fill(treeFolder8+"MC/8TeV/hzzTree_id100.root" , getBkgXsec(100)*lumi8/evt_8TeV(100), 0.0, false);

    ofstream file8;
    file8.open("cat0/info8TeV.txt");

    hmpi8.analyze(0, file8, fitValues, fitErrors);
    hmpi8.analyze(1, file8, fitValues, fitErrors);
    hmpi8.analyze(2, file8, fitValues, fitErrors);
    
    file8.close();

}

void getParamsSig() {

  float masses[20] =  {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150};
  float massesE[20] = {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
  float masses7TeV[5] =  {115,120,130,140,150};
  float massesE7TeV[5] = {  0,  0,  0,  0,  0};
  double fitValues[5];
  double fitErrors[5];
  float ldmVal[20],ldsVal[20],lnmVal[20],lnkVal[20],f1Val[20];
  float ldmErr[20],ldsErr[20],lnmErr[20],lnkErr[20],f1Err[20];

  float lumi7     = 5.05;
  float lumi8     = 19.6;

  std::string treeFolder7 = "root://pcmssd12//data/hzz4l/step2/HZZ4L_42X_S1_V13_S2_V03/";
  std::string treeFolder8 = "root://pcmssd12//data/hzz4l/step2/HZZ4L_53X_S1_V13_S2_V03/";

  enum zxalgo { ss=0, os=1};
  float Z1min[2] = { 40,   81.186 };
  float Z1max[2] = { 120, 101.186 };

  ofstream file8;
  file8.open("cat0/info8TeV_sig.txt");
  
  for(int channel=0; channel<3; ++channel) {
    for(int i=0; i<20; ++i) {
      cout << "PROCESSING MASS " << masses[i] << endl;
      stringstream idgg, idvbf, idwztt;
      idgg << "1" << masses[i];
      idvbf << "2" << masses[i];
      idwztt << "3" << masses[i];
        
      HiggsMassPointInfo hmpi8;
      hmpi8.z1min    = 40.;
      hmpi8.z2min    = 12.;
      hmpi8.massLow  = 100;
      hmpi8.massHigh = 160;
      hmpi8.do7TeV   = false;
      hmpi8.melacut  = -1.;
      hmpi8.njetscut = -1;
      hmpi8.ptcut    = 0;
      hmpi8.ptmax    = 1000000.;
      hmpi8.njetveto = true;
        
      setmc(true);    
      init(hmpi8.do7TeV);
        
      hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id"+idgg.str()+".root" , getXsecggH(masses[i])*lumi8/evt_8TeV(1000+masses[i]), 0.0, false);
      hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id"+idvbf.str()+".root" , getXsecVBF(masses[i])*lumi8/evt_8TeV(2000+masses[i]), 0.0, false);
      hmpi8.ymaker_qqzz.fill(treeFolder8+"MC/8TeV/hzzTree_id"+idwztt.str()+".root" , (getXsecWHi(masses[i])+getXsecZHi(masses[i])+getXsecttH(masses[i]))*lumi8/evt_8TeV(3000+masses[i]), 0.0, false);
        
      file8 << " ==> Mass MH = " << masses[i] << endl;
      hmpi8.analyze(channel, file8, fitValues, fitErrors, true, masses[i]);
        
      file8 << " ======================\n\n" << endl;
  
      cout << "landau m value,error: " << fitValues[0] << " , " << fitErrors[0] << endl;
      ldmVal[i]=fitValues[0]; ldmErr[i]=fitErrors[0];
        
      cout << "landau s value,error: " << fitValues[1] << " , " << fitErrors[1] << endl;
      ldsVal[i]=fitValues[1]; ldsErr[i]=fitErrors[1];
        
      if(channel==0) {
        cout << "logN m value,error: " << fitValues[2] << " , " << fitErrors[2] << endl;
        lnmVal[i]=fitValues[2]; lnmErr[i]=fitErrors[2];
          
        cout << "logN k value,error: " << fitValues[3] << " , " << fitErrors[3] << endl;
        lnkVal[i]=fitValues[3]; lnkErr[i]=fitErrors[3];
      } else {
        cout << "gauss m value,error: " << fitValues[2] << " , " << fitErrors[2] << endl;
        lnmVal[i]=fitValues[2]; lnmErr[i]=fitErrors[2];
          
        cout << "gauss sigma value,error: " << fitValues[3] << " , " << fitErrors[3] << endl;
        lnkVal[i]=fitValues[3]; lnkErr[i]=fitErrors[3];        
      }
        
      cout << "landau f value,error: " << fitValues[4] << " , " << fitErrors[4] << endl;
      f1Val[i]=fitValues[4]; f1Err[i]=fitErrors[4];

    }
    
    stringstream namefile;
    namefile << "parameters_masserr_channel" << channel << ".root";
    TFile *resultfile = TFile::Open(namefile.str().c_str(),"RECREATE");
      
    TGraphErrors *gLdM = new TGraphErrors(20,masses,ldmVal,massesE,ldmErr);
    TGraphErrors *gLdS = new TGraphErrors(20,masses,ldsVal,massesE,ldsErr);
    TGraphErrors *gLnM = new TGraphErrors(20,masses,lnmVal,massesE,lnmErr);
    TGraphErrors *gLnK = new TGraphErrors(20,masses,lnkVal,massesE,lnkErr);
    TGraphErrors *gF1 = new TGraphErrors(20,masses,f1Val,massesE,f1Err);
      
    gLdM->SetName("gLdM");  gLdM->SetMarkerStyle(20);   gLdM->SetMarkerSize(1);
    gLdS->SetName("gLdS");  gLdS->SetMarkerStyle(20);   gLdS->SetMarkerSize(1);
    gLnM->SetName("gLnM");  gLnM->SetMarkerStyle(20);   gLnM->SetMarkerSize(1);
    gLnK->SetName("gLnK");  gLnK->SetMarkerStyle(20);   gLnK->SetMarkerSize(1);
    gF1->SetName("gF1");    gF1->SetMarkerStyle(20);    gF1->SetMarkerSize(1);
      
    gLdM->SetTitle("");
    gLdM->GetXaxis()->SetTitle("mass (GeV)");
    gLdM->GetYaxis()->SetTitle("Landau m parameter");
      
    gLdS->SetTitle("");
    gLdS->GetXaxis()->SetTitle("mass (GeV)");
    gLdS->GetYaxis()->SetTitle("Landau #sigma parameter");
      
    if(channel==0) {
      gLnM->SetTitle("");
      gLnM->GetXaxis()->SetTitle("mass (GeV)");
      gLnM->GetYaxis()->SetTitle("LogN m parameter");
        
      gLnK->SetTitle("");
      gLnK->GetXaxis()->SetTitle("mass (GeV)");
      gLnK->GetYaxis()->SetTitle("LogN K parameter");
    } else {
      gLnM->SetTitle("");
      gLnM->GetXaxis()->SetTitle("mass (GeV)");
      gLnM->GetYaxis()->SetTitle("Gaussian m parameter");
        
      gLnK->SetTitle("");
      gLnK->GetXaxis()->SetTitle("mass (GeV)");
      gLnK->GetYaxis()->SetTitle("Gaussian #sigma parameter");
    }
      
    gF1->SetTitle("");
    gF1->GetXaxis()->SetTitle("mass (GeV)");
    gF1->GetYaxis()->SetTitle("f_{1} parameter");
      
    resultfile->cd();
    gLdM->Fit("pol2"); gLdM->Draw("Ap"); gLdM->Write();
    gLdS->Fit("pol2"); gLdS->Draw("Ap"); gLdS->Write();
    gLnM->Fit("pol2"); gLnM->Draw("Ap"); gLnM->Write();
    gLnK->Fit("pol2"); gLnK->Draw("Ap"); gLnK->Write();
    gF1->Fit("pol2");  gF1->Draw("Ap");  gF1->Write();
    resultfile->Close();
  }

  /*
  ofstream file7;
  file7.open("cat0/info7TeV_sig.txt");
  
  for(int channel=0; channel<3; ++channel) {
    for(int i=0; i<5; ++i) {
      cout << "PROCESSING MASS " << masses7TeV[i] << endl;
      stringstream idgg, idvbf, idwztt;
      idgg << "1" << masses7TeV[i];
      idvbf << "2" << masses7TeV[i];
      idwztt << "3" << masses7TeV[i];
        
      HiggsMassPointInfo hmpi7;
      hmpi7.z1min    = 40.;
      hmpi7.z2min    = 12.;
      hmpi7.massLow  = 100;
      hmpi7.massHigh = 160;
      hmpi7.do7TeV   = true;
      hmpi7.melacut  = -1.;
      hmpi7.njetscut = -1;
      hmpi7.ptcut    = 0;
      hmpi7.ptmax    = 1000000.;
      hmpi7.njetveto = true;
        
      setmc(true);    
      init(hmpi7.do7TeV);
        
      hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id"+idgg.str()+".root" , getXsecggH(masses7TeV[i])*lumi7/evt_7TeV(1000+masses7TeV[i]), 0.0, false);
      hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id"+idvbf.str()+".root" , getXsecVBF(masses7TeV[i])*lumi7/evt_7TeV(2000+masses7TeV[i]), 0.0, false);
      hmpi7.ymaker_qqzz.fill(treeFolder7+"MC/7TeV/hzzTree_id"+idwztt.str()+".root" , (getXsecWHi(masses7TeV[i])+getXsecZHi(masses7TeV[i])+getXsecttH(masses7TeV[i]))*lumi7/evt_7TeV(3000+masses7TeV[i]), 0.0, false);
        
      file7 << " ==> Mass MH = " << masses7TeV[i] << endl;
      hmpi7.analyze(channel, file7, fitValues, fitErrors, true, masses7TeV[i]);
        
      file7 << " ======================\n\n" << endl;
  
      cout << "landau m value,error: " << fitValues[0] << " , " << fitErrors[0] << endl;
      ldmVal[i]=fitValues[0]; ldmErr[i]=fitErrors[0];
        
      cout << "landau s value,error: " << fitValues[1] << " , " << fitErrors[1] << endl;
      ldsVal[i]=fitValues[1]; ldsErr[i]=fitErrors[1];
        
      if(channel==0) {
        cout << "logN m value,error: " << fitValues[2] << " , " << fitErrors[2] << endl;
        lnmVal[i]=fitValues[2]; lnmErr[i]=fitErrors[2];
          
        cout << "logN k value,error: " << fitValues[3] << " , " << fitErrors[3] << endl;
        lnkVal[i]=fitValues[3]; lnkErr[i]=fitErrors[3];
      } else {
        cout << "gauss m value,error: " << fitValues[2] << " , " << fitErrors[2] << endl;
        lnmVal[i]=fitValues[2]; lnmErr[i]=fitErrors[2];
          
        cout << "gauss sigma value,error: " << fitValues[3] << " , " << fitErrors[3] << endl;
        lnkVal[i]=fitValues[3]; lnkErr[i]=fitErrors[3];        
      }
        
      cout << "landau f value,error: " << fitValues[4] << " , " << fitErrors[4] << endl;
      f1Val[i]=fitValues[4]; f1Err[i]=fitErrors[4];

    }
    
    stringstream namefile;
    namefile << "parameters_masserr_channel" << channel << ".root";
    TFile *resultfile = TFile::Open(namefile.str().c_str(),"RECREATE");
      
    TGraphErrors *gLdM = new TGraphErrors(5,masses7TeV,ldmVal,massesE7TeV,ldmErr);
    TGraphErrors *gLdS = new TGraphErrors(5,masses7TeV,ldsVal,massesE7TeV,ldsErr);
    TGraphErrors *gLnM = new TGraphErrors(5,masses7TeV,lnmVal,massesE7TeV,lnmErr);
    TGraphErrors *gLnK = new TGraphErrors(5,masses7TeV,lnkVal,massesE7TeV,lnkErr);
    TGraphErrors *gF1 = new TGraphErrors(5,masses7TeV,f1Val,massesE7TeV,f1Err);
      
    gLdM->SetName("gLdM");  gLdM->SetMarkerStyle(20);   gLdM->SetMarkerSize(1);
    gLdS->SetName("gLdS");  gLdS->SetMarkerStyle(20);   gLdS->SetMarkerSize(1);
    gLnM->SetName("gLnM");  gLnM->SetMarkerStyle(20);   gLnM->SetMarkerSize(1);
    gLnK->SetName("gLnK");  gLnK->SetMarkerStyle(20);   gLnK->SetMarkerSize(1);
    gF1->SetName("gF1");    gF1->SetMarkerStyle(20);    gF1->SetMarkerSize(1);
      
    gLdM->SetTitle("");
    gLdM->GetXaxis()->SetTitle("mass (GeV)");
    gLdM->GetYaxis()->SetTitle("Landau m parameter");
      
    gLdS->SetTitle("");
    gLdS->GetXaxis()->SetTitle("mass (GeV)");
    gLdS->GetYaxis()->SetTitle("Landau #sigma parameter");
      
    if(channel==0) {
      gLnM->SetTitle("");
      gLnM->GetXaxis()->SetTitle("mass (GeV)");
      gLnM->GetYaxis()->SetTitle("LogN m parameter");
        
      gLnK->SetTitle("");
      gLnK->GetXaxis()->SetTitle("mass (GeV)");
      gLnK->GetYaxis()->SetTitle("LogN K parameter");
    } else {
      gLnM->SetTitle("");
      gLnM->GetXaxis()->SetTitle("mass (GeV)");
      gLnM->GetYaxis()->SetTitle("Gaussian m parameter");
        
      gLnK->SetTitle("");
      gLnK->GetXaxis()->SetTitle("mass (GeV)");
      gLnK->GetYaxis()->SetTitle("Gaussian #sigma parameter");
    }
      
    gF1->SetTitle("");
    gF1->GetXaxis()->SetTitle("mass (GeV)");
    gF1->GetYaxis()->SetTitle("f_{1} parameter");
      
    resultfile->cd();
    gLdM->Fit("pol2"); gLdM->Draw("Ap"); gLdM->Write();
    gLdS->Fit("pol2"); gLdS->Draw("Ap"); gLdS->Write();
    gLnM->Fit("pol2"); gLnM->Draw("Ap"); gLnM->Write();
    gLnK->Fit("pol2"); gLnK->Draw("Ap"); gLnK->Write();
    gF1->Fit("pol2");  gF1->Draw("Ap");  gF1->Write();
    resultfile->Close();
  }
  */

}

