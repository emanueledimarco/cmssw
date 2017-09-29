#ifndef RecoEgamma_EgammaTools_EnergyScaleCorrectionDB_class_h
#define RecoEgamma_EgammaTools_EnergyScaleCorrectionDB_class_h
/// Read and get energy scale and smearings from conditions DB
// *  author Emanuele Di Marco

#include <TString.h>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <TChain.h>
#include <TRandom3.h>
#include <string>
#include <bitset> 

#include "CondFormats/EgammaObjects/interface/EgmCorrectorParameters.h"

class EnergyScaleCorrectionDB_class {
	
public:  

  class VariableValues {
  public:
    friend class EnergyScaleCorrectionDB_class;
    VariableValues();
    void setEta  (float fEta) {mEta=fEta;}
    void setR9   (float fR9)  {mR9=fR9;}
    void setET   (float fET)  {mET=fET;}
    void setGain (int fGain)  {mGain=fGain;}

    float mEta;
    float mR9;
    float mET;
    int mGain;
  };

  enum VarTypes {kEta, kR9, kET, kGain};

  bool doScale, doSmearings;
  
 public:
  // constructor from DB
 EnergyScaleCorrectionDB_class(): doScale(false), doSmearings(false), smearingType_(ECALELF) {}
  EnergyScaleCorrectionDB_class(){}; ///< dummy constructor needed in ElectronEnergyCalibratorRun2
  ~EnergyScaleCorrectionDB_class(void) {}
  
  // setup from conditions DB
  void setTokens(const edm::ParameterSet &iConfig);
  void setEventSetup(const edm::EventSetup &es);
  void initCorrector();

  //------------------------------ scales
  float ScaleCorrection(unsigned int runNumber, bool isEBEle, double R9Ele, double etaSCEle,
			double EtEle, unsigned int gainSeed=12) const; ///< method to get energy scale corrections
  
  float ScaleCorrectionUncertainty(unsigned int runNumber, bool isEBEle,
				   double R9Ele, double etaSCEle, double EtEle, unsigned int gainSeed, 
				   std::bitset<scAll> uncBitMask=scAll) const;///< method to get scale correction uncertainties: it is:
  /** 
   * bit 0 = stat
   * bit 1 = syst
   * but 2 = gain
   */  
  float getSmearingSigma(int runNumber, bool isEBEle, float R9Ele, float etaSCEle, float EtEle, unsigned int gainSeed, paramSmear_t par, float nSigma = 0.) const;
  float getSmearingSigma(int runNumber, bool isEBEle, float R9Ele, float etaSCEle, float EtEle, unsigned int gainSeed, float nSigma_rho, float nSigma_phi) const;
  float getSmearingRho(int runNumber, bool isEBEle, float R9Ele, float etaSCEle, float EtEle, unsigned int gainSeed) const; ///< public for sigmaE estimate  
  
 private:

  std::vector<float> fillVector(const std::vector<VarTypes>& fVarTypes, const VariableValues&) const;
  std::vector<VarTypes> mapping(const std::vector<std::string>& fNames) const;
  std::vector<float> getParameters(EgmCorrectorParameters *parameters, EnergyScaleCorrectionDB_class::VariableValues& iValues) const;

  const EgmCorrectorParameters *scaleCorrections_;
  const EgmCorrectorParameters *resolutionCorrections_;
  
  std::string energyScaleKey_;
  std::string energyResolutionKey_;

  std::vector<VarTypes> mScaleParTypes,mScaleBinTypes;
  std::vector<VarTypes> mResolutionParTypes,mResolutionBinTypes;
};


#endif
