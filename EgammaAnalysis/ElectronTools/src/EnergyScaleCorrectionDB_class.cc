#include <cassert>
#include <stdlib.h>
#include <float.h> 
#include <iomanip>
#include <sstream>

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CondFormats/DataRecord/interface/EgmCorrectorParametersRcd.h"
#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrectionDB_class.h"

void EnergyScaleCorrectionDB_class::setTokens(const edm::ParameterSet &iConfig) {
  energyScaleKey_ = iConfig.getParameter<std::string>("electronScaleKey");
  energyResolutionKey_ = iConfig.getParameter<std::string>("electronResolutionKey");
}

void EnergyScaleCorrectionDB_class::setEventSetup(const edm::EventSetup &es) {
  edm::ESHandle<EgmCorrectorParameters> energyScaleH;
  edm::ESHandle<EgmCorrectorParameters> energyResolutionH;

  es.get<EgmCorrectorParametersRcd>.get(energyScaleKey_, energyScaleH);
  es.get<EgmCorrectorParametersRcd>.get(energyResolutionKey_, energyResolutionH);

  scaleCorrections_ = energyScaleH.product();
  resolutionCorrections_ = energyResolutionH.product();  
}

void EnergyScaleCorrectionDB_class::initCorrectors() {
  mScaleBinTypes.push_back(mapping(scaleCorrections_->parameters().definitions().binVar()));
  mScaleParTypes.push_back(mapping(scaleCorrections_->parameters().definitions().parVar()));
  mResolutionBinTypes.push_back(mapping(resolutionCorrections_->parameters().definitions().binVar()));
  mResolutionParTypes.push_back(mapping(resolutionCorrections_->parameters().definitions().parVar()));
}


std::vector<EnergyScaleCorrectionDB_class::VarTypes> EnergyScaleCorrectionDB_class::mapping(const std::vector<std::string>& fNames) const {
  std::vector<int> result;
  for(unsigned i=0;i<fNames.size();i++) {
    std::string ss = fNames[i];
    if (ss=="eta")
      result.push_back(kEta);
    else if(ss=="r9")
      result.push_back(kPt);
    else if(ss=="ET")
      result.push_back(kET);
    else if(ss=="gain")
      result.push_back(kGain);
    else {
      std::stringstream sserr;
      sserr<<"unknown parameter name: "<<ss;
      handleError("FactorizedJetCorrectorCalculator",sserr.str());
    }
  }
  return result;
}

std::vector<float> EnergyScaleCorrectionDB_class::fillVector(const std::vector<VarTypes>& fVarTypes, const VariableValues& iValues) {
  std::vector<float> result;
  for(unsigned i=0;i<fVarTypes.size();i++) {
    if (fVarTypes[i] == kEta) {
      result.push_back(iValues.mEta);
    } else if (fVarTypes[i] == kR9) {
      result.push_back(iValues.mR9);
    } else if (fVarTypes[i] == kET) {
      result.push_back(iValues.mET);
    } else if (fVarTypes[i] == kGain) {
      result.push_back(iValues.mGain);
    }   else {
      std::stringstream sserr;
      sserr<<"unknown parameter "<<fVarTypes[i];
    }
  }
  return result;
}

std::vector<float> EnergyScaleCorrectionDB_class::getParameters(EgmCorrectorParameters *parameters, EnergyScaleCorrectionDB_class::VariableValues& iValues) const 
{
  std::vector<float> results;
  std::vector<float> vx = fillVector(mBinTypes[i],iValues);
  std::vector<float> vy = fillVector(mParTypes[i],iValues);
  int bin = mParameters.binIndex(vX);
  if (bin<0) return results;
  const std::vector<float> par = parameters->record(fBin).parameters();
  return par;
}

float EnergyScaleCorrectionDB_class::ScaleCorrection(double R9Ele, double etaSCEle, double EtEle, unsigned int gainSeed) const
{
  double correction = 1;
  if(doScale==false) return correction;
  VariableValues iValues;
  iValues.setR9(R9Ele);
  iValues.setEta(etaSCEle);
  iValues.setET(EtEle);
  iValues.setGain(gainSeed);
  std::vector<float> parameters = getParameters(scaleCorrections_,iValues);
  return parameters[0];
}

float EnergyScaleCorrectionDB_class::ScaleCorrectionUncertainty(unsigned int runNumber, bool isEBEle,
							      double R9Ele, double etaSCEle, double EtEle, unsigned int gainSeed, std::bitset<scAll> uncBitMask) const
{
  double totUncertainty(0);
  VariableValues iValues;
  iValues.setR9(R9Ele);
  iValues.setEta(etaSCEle);
  iValues.setET(EtEle);
  iValues.setGain(gainSeed);
  std::vector<float> parameters = getParameters(scaleCorrections_,iValues);

  for (unsigned int p=1;p<parameters.size();++i) {
    if(uncBitMask.test(p-1)) totUncertainty += parameters[p]*parameters[p];
  }
  return sqrt(totUncertainty);
}

float EnergyScaleCorrectionDB_class::getSmearingSigma(float R9Ele, float etaSCEle, float EtEle, unsigned int gainSeed, paramSmear_t par, float nSigma) const
{
  if (par == kRho) return getSmearingSigma(runNumber, isEBEle, R9Ele, etaSCEle, EtEle, gainSeed, nSigma, 0.);
  if (par == kPhi) return getSmearingSigma(runNumber, isEBEle, R9Ele, etaSCEle, EtEle, gainSeed, 0., nSigma);
  return getSmearingSigma(runNumber, isEBEle, R9Ele, etaSCEle, EtEle, gainSeed, 0., 0.);
}

float EnergyScaleCorrectionDB_class::getSmearingSigma(float R9Ele, float etaSCEle, float EtEle, unsigned int gainSeed, float nSigma_rho, float nSigma_phi) const
{  
  double correction = 1;
  if(doScale==false) return correction;
  VariableValues iValues;
  iValues.setR9(R9Ele);
  iValues.setEta(etaSCEle);
  iValues.setET(EtEle);
  iValues.setGain(gainSeed);
  std::vector<float> parameters = getParameters(resolutionCorrections_,iValues);

  double rho = parameters[2] + parameters[3] * nSigma_rho;
  double phi = parameters[4] + parameters[5] * nSigma_phi;

  double constTerm =  rho * sin(phi);
  double alpha =  rho *  parameters[0] * cos( phi);

  return sqrt(constTerm * constTerm + alpha * alpha / EtEle); 
}

float EnergyScaleCorrectionDB_class::getSmearingRho(float R9Ele, float etaSCEle, float EtEle, unsigned int gainSeed) const
{
  double correction = 1;
  if(doScale==false) return correction;
  VariableValues iValues;
  iValues.setR9(R9Ele);
  iValues.setEta(etaSCEle);
  iValues.setET(EtEle);
  iValues.setGain(gainSeed);
  std::vector<float> parameters = getParameters(resolutionCorrections_,iValues);
  return parameters[2];
}

