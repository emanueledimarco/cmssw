#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitRecWeightsAlgo_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitRecWeightsAlgo_HH

/** \class EcalUncalibRecHitRecWeightsAlgo
  *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
  *  using a weights method
  *
  *  \author R. Bruneliere - A. Zabi
  *  
  *  The chi2 computation with matrix is replaced by the chi2express which is  moved outside the weight algo
  *  (need to clean up the interface in next iteration so that we do not pass-by useless arrays)
  *
  */

#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalShapeBase.h"

#include <vector>

template<class C> class EcalUncalibRecHitRecWeightsAlgo 
{
 private:
  std::vector<double> pileupDataFrame;
  bool subtract_pu;
  bool dyn_pedestal;

 public:
  // constructor
  EcalUncalibRecHitRecWeightsAlgo<C>() { 
    subtract_pu = false; 
    dyn_pedestal = true;
  };

  // destructor
  virtual ~EcalUncalibRecHitRecWeightsAlgo<C>() { };

  /// set the pileup dataframe to subtract
  void setPileupDataFrame(std::vector<double> &dataFrame) { 
    pileupDataFrame = dataFrame; 
    subtract_pu = true;
  }

  void SetDynamicPedestal(bool dyn_pede) { dyn_pedestal = dyn_pede; }

  /// Compute parameters
   virtual EcalUncalibratedRecHit makeRecHit(
					      const C& dataFrame 
					      , const double* pedestals
					      , const double* pedestalsRMS
					      , const double* gainRatios 
					      , const EcalWeightSet::EcalWeightMatrix** weights
					      , const EcalShapeBase & testbeamPulseShape
    ) {
    double amplitude_(-1.),  pedestal_(-1.), jitter_(-1.), chi2_(-1.);
    uint32_t flag = 0;

    float pedestalFromDB=0;

    // Get time samples
    ROOT::Math::SVector<double,C::MAXSAMPLES> frame, framePedSub;
    int gainId0 = 1;
    int iGainSwitch = 0;
    bool isSaturated = 0;
    for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
      int gainId = dataFrame.sample(iSample).gainId();
      //Handling saturation (treating saturated gainId as maximum gain)
      if ( gainId == 0 ) 
	{ 
	  gainId = 3;
	  //isSaturated = 1;
	  // in pileup run May2012 samples 7,8,9,10 have gainid ==0
          // fix it like this: it won't hurt for the future SA20120512
	  if(iSample==4 || iSample ==5 || iSample==6) isSaturated = 1;
	}

      //      if (gainId != gainId0) iGainSwitch = 1;
      // same problem as above: mark saturation only when physically
      // expected to occur SA20120513
      if ( (gainId != gainId0) && (iSample==4 || iSample ==5 || iSample==6) ) iGainSwitch = 1;
      if (!iGainSwitch)
	frame(iSample) = double(dataFrame.sample(iSample).adc());
      else
	frame(iSample) = double(((double)(dataFrame.sample(iSample).adc()) - pedestals[gainId-1]) * gainRatios[gainId-1]);

      if(subtract_pu) frame(iSample) -= pileupDataFrame[iSample];

      float gainratio = 1.;
      // now calculate the amplitude ped-subtracted for the amplitude
      if (gainId==0 || gainId==3) {
        pedestalFromDB = pedestals[2];
        gainratio = gainRatios[2];
      } else if (gainId==1) {
        pedestalFromDB = pedestals[0];
        gainratio = 1.;
      } else if(gainId==2) {
        pedestalFromDB = pedestals[1];
        gainratio = gainRatios[1];
      }
      framePedSub(iSample) = double(((double)(dataFrame.sample(iSample).adc()) - pedestalFromDB) * gainratio);
      if (gainId == 0)  framePedSub(iSample) = double((4095. - pedestalFromDB) * gainratio);

      if(subtract_pu) framePedSub(iSample) -= pileupDataFrame[iSample];

    }

    /*    
    std::cout << "=====In the weights algo. ==== " << std::endl;
    std::cout << "PU SUB = " << subtract_pu << std::endl;
    for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++)
      std::cout << "\tisample = " << iSample << " frame(iSample) = " << frame(iSample) << "\tframePedSub(iSample) = " << framePedSub(iSample) << std::endl;
    std::cout << "=====Done weights algo. ==== " << std::endl;
    */

    // Compute amplitude
    double amplitudeDBPed=0.;
    for(int iSample = 3; iSample < C::MAXSAMPLES; iSample++) 
      amplitudeDBPed += (weights[iGainSwitch])->At(0,iSample) * framePedSub(iSample);


    // Compute other parameters
    ROOT::Math::SVector <double,3> param = (*(weights[iGainSwitch])) * frame;
    double amplitudeDynPed = param(EcalUncalibRecHitRecAbsAlgo<C>::iAmplitude);
    double pedestalDyn = param(EcalUncalibRecHitRecAbsAlgo<C>::iPedestal);
    if (amplitude_) jitter_ = -param(EcalUncalibRecHitRecAbsAlgo<C>::iTime) / amplitude_;
    else jitter_ = 0.;

    amplitude_ = dyn_pedestal ? amplitudeDynPed : amplitudeDBPed;
    pedestal_ = dyn_pedestal ? pedestalDyn : pedestalFromDB;

    //When saturated gain flag i
    if (isSaturated)
      {
        flag = EcalUncalibratedRecHit::kSaturated;
	amplitude_ = double((4095. - pedestals[2]) * gainRatios[2]);
      }
    return EcalUncalibratedRecHit( dataFrame.id(), amplitude_, pedestal_, jitter_, chi2_, flag);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};
#endif
