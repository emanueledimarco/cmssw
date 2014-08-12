#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitOutOfTimeSubtractionAlgo_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitOutOfTimeSubtractionAlgo_HH

/** \class EcalUncalibRecHitOutOfTimeSubtractionAlgo
 *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
 *  using a ratio method
 *
 *  $Id: EcalUncalibRecHitOutOfTimeSubtractionAlgo.h,v 1.50 2012/06/11 21:02:13 wmtan Exp $
 *  $Date: 2012/06/11 21:02:13 $
 *  $Revision: 1.50 $
 *  \author E. Di Marco
 */

#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"
#include <vector>

template < class C > class EcalUncalibRecHitOutOfTimeSubtractionAlgo {
 public:

  struct CalculatedExtraHit {
    double amplitudeMax;
    double timeMax;
    double amplitudeExtapolated;
  };

  virtual ~ EcalUncalibRecHitOutOfTimeSubtractionAlgo < C > () { };
  // function to be able to compute
  // amplitude and time of the OOT pileup
  void init( const C &dataFrame, std::vector<C> neighbors, 
             const double *pedestals, const double *gainRatios, 
             std::vector< double > &amplitudeFitParameters, std::vector< double > &subtractionLimit );
  double pulseShapeFunction(double t);
  void computeAmplitudeOOT(std::vector<double> &extraHitFrame);
  
 protected:
  
  DetId          theDetId_;
  std::vector < double > amplitudes_;
  double pedestal_;
  double alpha_, beta_, alphabeta_;
  int nSameTimeHits_;
  float minAmplitudeOutOfTime_, lastOOTSampleUsed_;

  CalculatedExtraHit calculatedExtrahit_;
};

template<class C> double EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::pulseShapeFunction(double t){
  if( alphabeta_ <= 0 ) return((double)0.);
  double t0 = calculatedExtrahit_.timeMax;
  //  std::cout << "alphabeta_ = " << alphabeta_ << "   t0 = " << t0 << "   t = " << t << std::endl;
  double dtsbeta,variable,puiss;
  double dt = t-t0 ;
  if(dt > -alphabeta_)  {
    // std::cout << "deltat = " << dt << std::endl;
    dtsbeta=dt/beta_ ;
    variable=1.+dt/alphabeta_ ;
    puiss=pow(variable,alpha_);
    //    std::cout << "the normalization at t0 = " << t0 << " is A_out = " << amplitudes_[calculatedExtrahit_.timeMax] << std::endl;
    int intTime = int(calculatedExtrahit_.timeMax);
    if(intTime<0) intTime = 0;
    if(intTime>9) intTime = 9;
    return amplitudes_[intTime]*puiss*exp(-dtsbeta) ;
  }
  return  0.0 ;
}

template <class C>
void EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::init( const C &dataFrame, std::vector<C> neighbors, 
                                                         const double *pedestals, const double *gainRatios, 
                                                         std::vector< double > &amplitudeFitParameters, std::vector< double > &subtractionLimits ) {
    
  // these are the parameters of the pulse shape function
  alphabeta_ = amplitudeFitParameters[0]*amplitudeFitParameters[1];
  alpha_ = amplitudeFitParameters[0];
  beta_ = amplitudeFitParameters[1];

  // these are the limits of the subtraction algorithm
  lastOOTSampleUsed_ = subtractionLimits[0];
  minAmplitudeOutOfTime_ = subtractionLimits[1];
  
  theDetId_ = DetId(dataFrame.id().rawId());  
  calculatedExtrahit_.timeMax = 5;
  calculatedExtrahit_.amplitudeMax = 0;
  amplitudes_.clear();
  amplitudes_.reserve(C::MAXSAMPLES);
	
  // obtain the max sample in the NxN matrix. 
  // pedestal is super-simple: 200 ADC (we are interested only in the time of the max here)
  int max5x5Ampli=-100;
  int max5x5Sample=5;
  int timeNeighbors[neighbors.size()];
  int index=0;
  for(typename std::vector<C>::const_iterator itn=neighbors.begin(); itn!=neighbors.end(); ++itn, ++index) {
    int thisHitMaxAmpli=-100;
    int thisHitMaxSample=5;
    for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
      int gainId = itn->sample(iSample).gainId(); 

      int  sampleAdc = 0;
      if ( gainId == 1 ) sampleAdc = itn->sample(iSample).adc() - 200;
      else if ( gainId == 2 ) sampleAdc = (itn->sample(iSample).adc() - 200) * 2 ;
      else sampleAdc = (dataFrame.sample(iSample).adc() - 200) * 12 ;
      
      // hit max
      if( sampleAdc > thisHitMaxAmpli ) {
        thisHitMaxAmpli = sampleAdc;
        thisHitMaxSample = iSample;
      }

      // global max
      if( sampleAdc > max5x5Ampli ) {
        max5x5Ampli  = sampleAdc;
        max5x5Sample = iSample;
      }
    }// loop on samples
    timeNeighbors[index] = (thisHitMaxAmpli > minAmplitudeOutOfTime_) ? thisHitMaxSample : -900;
  } // loop over neighbors

  nSameTimeHits_=0;
  for(int iHit=0; iHit<neighbors.size(); ++iHit)
    if(timeNeighbors[iHit]>-500 && timeNeighbors[iHit]==max5x5Sample) nSameTimeHits_++;

  if(neighbors.size()>0) {
    calculatedExtrahit_.amplitudeMax = max5x5Ampli;
    calculatedExtrahit_.timeMax = max5x5Sample;
  } else calculatedExtrahit_.amplitudeMax = -100.;

  // now compute the samples for the reconstructed rechit 
  // with pedestal subtraction scheme that is consistent with the local reco 
  for (int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
    
    const EcalMGPASample &sample = dataFrame.sample(iSample);

    double amplitude = 0.;
    int gainId = sample.gainId();

    double gainratio = 1.;
    float pedestal=0;
    if (gainId==0 || gainId==3) {
      pedestal = pedestals[2];
      gainratio = gainRatios[2];
    } else if (gainId==1) {
      pedestal = pedestals[0];
      gainratio = 1.;
    } else if(gainId==2) {
      pedestal = pedestals[1];
      gainratio = gainRatios[1];
    }

    amplitude = double(((double)(sample.adc()) - pedestal) * gainratio);
    
    if (gainId == 0)  amplitude = double((4095. - pedestal) * gainratio);

    amplitudes_.push_back(amplitude);
  }
}

template<class C>
void EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::computeAmplitudeOOT(std::vector<double> &extraHitFrame)
{

  extraHitFrame.resize(C::MAXSAMPLES);

  // if the max sample of the max of the neighbors is 
  // in the pre-samples, and the amplitude is above noise, 
  // --> there is a hope to fit for the out-of-time contributions
  // -> not, return 0
  if(calculatedExtrahit_.amplitudeMax == -100 ||
     calculatedExtrahit_.timeMax > lastOOTSampleUsed_ ||
     calculatedExtrahit_.amplitudeMax < minAmplitudeOutOfTime_ ) {
    for(int iSample = 3; iSample < C::MAXSAMPLES; iSample++) {
      extraHitFrame[iSample] = .0;
    }
  } else {    

    /*
    std::cout << "===> timeMax = " << calculatedExtrahit_.timeMax << "\tcalculatedExtrahit_.amplitudeMax = " << calculatedExtrahit_.amplitudeMax << std::endl;
    std::cout << "===> # hits with consistent time = " << nSameTimeHits_ << std::endl;
    std::cout << "listing the samples ... " << std::endl;
    for (int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
      std::cout << "\tisample " << iSample << " has amplitudes_[iSample] = " << amplitudes_[iSample] << std::endl;
    }
    std::cout << "Done listing samples " << std::endl;
    */
    
    for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
      extraHitFrame[iSample] = (iSample<3) ? 0.0 : std::max(.0, pulseShapeFunction(iSample));
      // std::cout << "\t\tSUBTRACTING TO SAMPLE " << iSample << " extraHitFrame[iSample] = " << extraHitFrame[iSample] << std::endl;
    }

  }
}


#endif
