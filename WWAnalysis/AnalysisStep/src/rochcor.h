#ifndef WWAnalysis_AnalysisStep_RochCor_h
#define WWAnalysis_AnalysisStep_RochCor_h

////  VERSION V4.1 DOWNLOADED FROM
////     http://www-cdf.fnal.gov/~jyhan/cms_momscl/rochcor_v4_new.tgz
////  ON 12/12/2012 (although not at 12:12, what a pity...)

#include <iostream>
#include <TChain.h>
#include <TClonesArray.h>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>


class rochcor {
 public:
  rochcor();
  rochcor(int seed);
  ~rochcor();
  
  void momcor_mc(TLorentzVector&, float, float, int, float&);
  void momcor_data(TLorentzVector&, float, float, int, float&);
  
  float zptcor(float);
  int etabin(float);
  int phibin(float);
  
 private:
  
  TRandom3 eran;
  TRandom3 sran;
  
  
  //  static float netabin[9] = {-2.4,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.4};
  static const float netabin[9];
  
  static const double pi;
  static const float genm_smr;
  static const float genm;
  
  static const float recmA;
  static const float drecmA;
  static const float mgsclA_stat;
  static const float mgsclA_syst;
  static const float dgsclA_stat;
  static const float dgsclA_syst;
  static const float recmB;
  static const float drecmB;
  static const float mgsclB_stat;
  static const float mgsclB_syst;
  static const float dgsclB_stat;
  static const float dgsclB_syst;
  
  //iteration2 after FSR : after Z Pt correction
  static const float deltaA;
  static const float deltaA_stat;
  static const float deltaA_syst;
  
  static const float sfA;
  static const float sfA_stat;
  static const float sfA_syst;
  
  static const float deltaB;
  static const float deltaB_stat;
  static const float deltaB_syst;
  
  static const float sfB;
  static const float sfB_stat;
  static const float sfB_syst;

  static const float apar;
  static const float bpar;
  static const float cpar;
  static const float d0par;
  static const float e0par;
  static const float d1par;
  static const float e1par;
  static const float d2par;
  static const float e2par;
 
  //---------------------------------------------------------------------------------------------
  
  static const float dcor_bfA[8][8];  
  static const float dcor_maA[8][8];
  static const float mcor_bfA[8][8];
  static const float mcor_maA[8][8];
  static const float dcor_bfAer[8][8];  
  static const float dcor_maAer[8][8];
  static const float mcor_bfAer[8][8];
  static const float mcor_maAer[8][8];

  static const float dcor_bfB[8][8];  
  static const float dcor_maB[8][8];
  static const float mcor_bfB[8][8];
  static const float mcor_maB[8][8];
  static const float dcor_bfBer[8][8];  
  static const float dcor_maBer[8][8];
  static const float mcor_bfBer[8][8];
  static const float mcor_maBer[8][8];

  //=======================================================================================================
  
  static const float dmavgA[8][8];  
  static const float dpavgA[8][8];  
  static const float mmavgA[8][8];  
  static const float mpavgA[8][8];

  static const float dmavgB[8][8];  
  static const float dpavgB[8][8];  
  static const float mmavgB[8][8];  
  static const float mpavgB[8][8];
  
  //===============================================================================================
  //parameters for Z pt correction
  static const int nptbins=84;
  static const float ptlow[85];    
  
  static const float zptscl[84];
  static const float zptscler[84];

  float mptsys_mc_dm[8][8];
  float mptsys_mc_da[8][8];
  float mptsys_da_dm[8][8];
  float mptsys_da_da[8][8];

  float gscler_mc_dev;
  float gscler_da_dev;
};

#endif  
