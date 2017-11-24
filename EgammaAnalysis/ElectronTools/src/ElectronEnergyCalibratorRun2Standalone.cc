#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2Standalone.h"
#include <CLHEP/Random/RandGaussQ.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"

ElectronEnergyCalibratorRun2Standalone::ElectronEnergyCalibratorRun2Standalone(bool isMC, 
							   bool synchronization, 
							   std::string correctionFile
							   ) :
  isMC_(isMC), synchronization_(synchronization),
  rng_(0),
  doEpCombination_(false),
  _correctionRetriever(correctionFile) // here is opening the files and reading the corrections
{
  if(isMC_) {
    _correctionRetriever.doScale = false; 
    _correctionRetriever.doSmearings = true;
  } else {
    _correctionRetriever.doScale = true; 
    _correctionRetriever.doSmearings = false;
  }
}

ElectronEnergyCalibratorRun2Standalone::~ElectronEnergyCalibratorRun2Standalone()
{}

void ElectronEnergyCalibratorRun2Standalone::initPrivateRng(TRandom *rnd) 
{
  rng_ = rnd;   
}

void ElectronEnergyCalibratorRun2Standalone::calibrate(SimpleElectron &electron) const 
{
  assert(isMC_ == electron.isMC());
  float smear = 0.0, scale = 1.0;
  float aeta = std::abs(electron.getEta()); //, r9 = electron.getR9();
  float et = electron.getNewEnergy()/cosh(aeta);
  
  scale = _correctionRetriever.ScaleCorrection(electron.getRunNumber(), electron.isEB(), electron.getR9(), aeta, et);
  smear = _correctionRetriever.getSmearingSigma(electron.getRunNumber(), electron.isEB(), electron.getR9(), aeta, et, 0., 0.); 

  double newEcalEnergy, newEcalEnergyError;
  if (isMC_) {
    double corr = 1.0 + smear * gauss();
    newEcalEnergy      = electron.getNewEnergy() * corr;
    newEcalEnergyError = std::hypot(electron.getNewEnergyError() * corr, smear * newEcalEnergy);
  } else {
    newEcalEnergy      = electron.getNewEnergy() * scale;
    newEcalEnergyError = std::hypot(electron.getNewEnergyError() * scale, smear * newEcalEnergy);
  }
  electron.setNewEnergy(newEcalEnergy); 
  electron.setNewEnergyError(newEcalEnergyError);
}

double ElectronEnergyCalibratorRun2Standalone::gauss() const 
{
  if (synchronization_) return 1.0;
  if (rng_) {
    return rng_->Gaus();
  } else {
    std::cout << "ERROR! TRandom not initialized! Returning correction 0" << std::endl;
    return 0.;
  }
}

