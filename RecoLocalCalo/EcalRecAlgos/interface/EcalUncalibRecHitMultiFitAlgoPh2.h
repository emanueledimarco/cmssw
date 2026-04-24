#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgoPh2_h
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgoPh2_h

/** \class EcalUncalibRecHitMultiFitAlgoPh2
  *  Amplitude reconstucted from Phase 2 digis by the multi-template fit
  */

#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/PiecewiseCubicSpline.h"

class EcalUncalibRecHitMultiFitAlgoPh2 {
public:
  using SampleVector = typename EigenMatrixTypes<P>::SampleVector;
  using FullSampleVector = typename EigenMatrixTypes<P>::FullSampleVector;
  using BXVector = typename EigenMatrixTypes<P>::BXVector;
  using SampleGainVector = typename EigenMatrixTypes<P>::SampleGainVector;
  using SampleMatrix = typename EigenMatrixTypes<P>::SampleMatrix;
  using FullSampleMatrix = typename EigenMatrixTypes<P>::FullSampleMatrix;
  using SampleMatrixGainArray = typename EigenMatrixTypes<P>::SampleMatrixGainArray;

  EcalUncalibRecHitMultiFitAlgoPh2();
  ~EcalUncalibRecHitMultiFitAlgoPh2(){};
  EcalUncalibratedRecHit makeRecHit(const EcalDataFrame_Ph2 &dataFrame,
                                    const EcalLiteDTUPedestalsMap::Item *aped,
                                    const EcalCATIAGainRatio *aGain,
                                    const SampleMatrixGainArray &noisecors,
                                    const FullSampleVector &fullpulse,
                                    const FullSampleMatrix &fullpulsecov,
                                    const BXVector &activeBX
                                    const PiecewiseCubicSpline &spline);
  void disableErrorCalculation() { computeErrors_ = false; }
  void setDoPrefit(const bool b) { doPrefit_ = b; }
  void setPrefitMaxChiSq(const double x) { prefitMaxChiSq_ = x; }
  void setDynamicPedestals(const bool b) { dynamicPedestals_ = b; }
  void setMitigateBadSamples(const bool b) { mitigateBadSamples_ = b; }
  void setSelectiveBadSampleCriteria(const bool b) { selectiveBadSampleCriteria_ = b; }
  void setAddPedestalUncertainty(const double x) { addPedestalUncertainty_ = x; }
  void setSimplifiedNoiseModelForGainSwitch(const bool b) { simplifiedNoiseModelForGainSwitch_ = b; }
  void setGainSwitchUseMaxSample(const bool b) { gainSwitchUseMaxSample_ = b; }

private:
  PulseChiSqSNNLS<P> pulsefunc_;
  PulseChiSqSNNLS<P> pulsefuncSingle_;
  bool computeErrors_;
  bool doPrefit_;
  double prefitMaxChiSq_;
  bool dynamicPedestals_;
  bool mitigateBadSamples_;
  bool selectiveBadSampleCriteria_;
  double addPedestalUncertainty_;
  bool simplifiedNoiseModelForGainSwitch_;
  bool gainSwitchUseMaxSample_;
  BXVector singlebx_;
};



