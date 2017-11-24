#ifndef ElectronEnergyCalibratorRun2Standalone_h
#define ElectronEnergyCalibratorRun2Standalone_h

#include <TRandom.h>
#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.hh"
#include "EgammaAnalysis/ElectronTools/interface/SimpleElectronStandalone.h"

#include <vector>

class ElectronEnergyCalibratorRun2Standalone {
 public:
  // dummy constructor for persistence
  ElectronEnergyCalibratorRun2Standalone() {}
  
  // constructor w/o combinator: do not combine E-p: intented for residual corrections
  ElectronEnergyCalibratorRun2Standalone(bool isMC, bool synchronization, std::string); 
  ~ElectronEnergyCalibratorRun2Standalone() ;
  
  /// Initialize with a random number generator (if not done, it will use the CMSSW service)
  /// Caller code owns the TRandom.
  void initPrivateRng(TRandom *rnd) ;
  
  /// Correct this electron. 
  void calibrate(SimpleElectron &electron) const ;
  
 protected:    
  // whatever data will be needed
  bool isMC_, synchronization_;
  TRandom *rng_;
  
  /// Return a number distributed as a unit gaussian, drawn from the private RNG if initPrivateRng was called, 
  /// or from the CMSSW RandomNumberGenerator service
  /// If synchronization is set to true, it returns a fixed number (1.0)
  double gauss() const ;
  bool doEpCombination_;
  EnergyScaleCorrection_class _correctionRetriever;
};

#endif
