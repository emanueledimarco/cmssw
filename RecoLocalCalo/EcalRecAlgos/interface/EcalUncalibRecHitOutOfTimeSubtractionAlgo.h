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
    double amplitudeMax, amplitudeSecondMax;
    double timeMax, timeSecondMax;
    bool puTagged;
  };

  // constructor
  EcalUncalibRecHitOutOfTimeSubtractionAlgo < C > () {
    farTail_ = false;
  }

  // destructor 
  virtual ~ EcalUncalibRecHitOutOfTimeSubtractionAlgo < C > () { };

  CalculatedExtraHit calculatedExtrahit_;

  // function to be able to compute
  // amplitude and time of the OOT pileup
  void init( const C &dataFrame, std::vector<C> neighbors, 
             const double *pedestals, const double *gainRatios, 
             std::vector< double > &amplitudeFitParameters, std::vector< double > &subtractionLimit );
  double pulseShapeFunction(double t);
  double farTailFunction(double t);
  void computeAmplitudeOOT(std::vector<double> &extraHitFrame);

 protected:
  
  DetId          theDetId_;
  std::vector < double > amplitudes_;
  double pedestal_;
  double alpha_, beta_, alphabeta_;
  int minNConsistentBX_, nConsistentBX_, nConsistentBXSecond_;
  float minSeedAmplitudeOutOfTime_, minAmplitudeOutOfTime_, lastOOTSampleUsed_;

  // parameters for the endcap in case the max sample is the 0th
  bool farTail_;
  double slope_, offset_;

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

template<class C> double EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::farTailFunction(double t) {
  return slope_ * t + offset_;
}

template <class C>
void EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::init( const C &dataFrame, std::vector<C> neighbors, 
                                                         const double *pedestals, const double *gainRatios,
                                                         std::vector< double > &amplitudeFitParameters, std::vector< double > &subtractionParameters ) {
  
  // these are the parameters of the pulse shape function
  alphabeta_ = amplitudeFitParameters[0]*amplitudeFitParameters[1];
  alpha_ = amplitudeFitParameters[0];
  beta_ = amplitudeFitParameters[1];

  // these are the parameters of the subtraction algorithm
  minSeedAmplitudeOutOfTime_ = subtractionParameters[0];
  minAmplitudeOutOfTime_ = subtractionParameters[1];
  minNConsistentBX_ = subtractionParameters[2];
  lastOOTSampleUsed_ = subtractionParameters[3];
  
  // only for EE, assume a far tail if maxSample = 0
  slope_ = .0;
  offset_ = .0;

  theDetId_ = DetId(dataFrame.id().rawId());  
  calculatedExtrahit_.timeMax = 5;
  calculatedExtrahit_.amplitudeMax = 0;
  calculatedExtrahit_.timeSecondMax = 5;
  calculatedExtrahit_.amplitudeSecondMax = 0;
  calculatedExtrahit_.puTagged = false;
  farTail_ = false;

  amplitudes_.clear();
  amplitudes_.reserve(C::MAXSAMPLES);
	
  // obtain the max sample in the NxN matrix. 
  // pedestal is super-simple: 200 ADC (we are interested only in the time of the max here)
  int max5x5Ampli(-100), secondMax5x5Ampli(-100);
  int max5x5Sample(5), secondMax5x5Sample(5);
  int timeNeighbors[neighbors.size()];
  int index=0;
  double slopeMax(0.0), slopeSecondMax(0.0);
  for(typename std::vector<C>::const_iterator itn=neighbors.begin(); itn!=neighbors.end(); ++itn, ++index) {
    int thisHitMaxAmpli=-100;
    int thisHitMaxSample=5;
    double hitA(0.0),hitB(0.0);
    // double itPede = (double(itn->sample(0).adc()) + double(itn->sample(1).adc()) + double(itn->sample(2).adc()))/3.;
    for(int iSample = 0; iSample <= lastOOTSampleUsed_; iSample++) {
      int gainId = itn->sample(iSample).gainId(); 
      int  sampleAdc = 0;
      if ( gainId == 1 ) sampleAdc = itn->sample(iSample).adc() - 200;
      else if ( gainId == 2 ) sampleAdc = (itn->sample(iSample).adc() - 200) * 2 ;
      else sampleAdc = (dataFrame.sample(iSample).adc() - 200) * 12 ;
      
      // slope of the far tail (this is valid only for 3 points and used only for EE)
      if(dataFrame.id().subdetId()==EcalEndcap) {
        hitA += iSample * sampleAdc;
        hitB += sampleAdc;
      }

      // hit max
      if( sampleAdc > thisHitMaxAmpli ) {
        thisHitMaxAmpli = sampleAdc;
        thisHitMaxSample = iSample;
      }
      
      // global max
      if( sampleAdc > max5x5Ampli ) {
        max5x5Ampli  = sampleAdc;
        max5x5Sample = iSample;
      } else if( sampleAdc > secondMax5x5Ampli ) {
        secondMax5x5Ampli = sampleAdc;
        secondMax5x5Sample = iSample;
      }

    }// loop on samples
    timeNeighbors[index] = (thisHitMaxAmpli > minAmplitudeOutOfTime_) ? thisHitMaxSample : -900;
    if(thisHitMaxAmpli == max5x5Ampli) slopeMax =  0.5 * (hitA-hitB);
    else if(thisHitMaxAmpli == secondMax5x5Ampli) slopeSecondMax = 0.5 * (hitA-hitB);
    //    std::cout << "\t\tXXX neighbor " << index << " has maxTime = " << thisHitMaxSample << " and maxAmpli = " << thisHitMaxAmpli << std::endl;
  } // loop over neighbors

  nConsistentBX_ = 0; 
  nConsistentBXSecond_ = 0;
  for(int iHit=0; iHit<neighbors.size(); ++iHit) {
    if(timeNeighbors[iHit]>-500 && timeNeighbors[iHit]==max5x5Sample) nConsistentBX_++;
    if(timeNeighbors[iHit]>-500 && timeNeighbors[iHit]==secondMax5x5Sample) nConsistentBXSecond_++;
  }

  if(neighbors.size()>0) {
    calculatedExtrahit_.amplitudeMax = max5x5Ampli;
    calculatedExtrahit_.timeMax = max5x5Sample;
    calculatedExtrahit_.amplitudeSecondMax = secondMax5x5Ampli;
    calculatedExtrahit_.timeSecondMax = secondMax5x5Sample;
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

  // for the far tail subtraction
  if(dataFrame.id().subdetId()==EcalEndcap) {
    if(calculatedExtrahit_.timeMax==0 && 
       (calculatedExtrahit_.amplitudeMax >= minSeedAmplitudeOutOfTime_ && nConsistentBX_ >= minNConsistentBX_) ) {
      farTail_ = true;
      slope_ = slopeMax;
      offset_ = amplitudes_[0];
      //      std::cout << "far tail from the max sample detected: slope = " << slope_ << std::endl;
    } else if(calculatedExtrahit_.timeSecondMax==0 && 
              (calculatedExtrahit_.amplitudeSecondMax >= minSeedAmplitudeOutOfTime_ && nConsistentBXSecond_ >= minNConsistentBX_)) {
      farTail_ = true;
      slope_ = slopeSecondMax;
      offset_ = amplitudes_[0];
      // std::cout << "far tail from the second max sample detected: slope = " << slope_ << std::endl;
    }
  }

}

template<class C>
void EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::computeAmplitudeOOT(std::vector<double> &extraHitFrame)
{

  extraHitFrame.resize(C::MAXSAMPLES);

  /*
  if(calculatedExtrahit_.amplitudeMax >= minSeedAmplitudeOutOfTime_ && nConsistentBX_ >= minNConsistentBX_ ) {
      std::cout << "===> timeMax = " << calculatedExtrahit_.timeMax << "\tcalculatedExtrahit_.amplitudeMax = " << calculatedExtrahit_.amplitudeMax << std::endl;
      std::cout << "===> # hits with consistent time = " << nConsistentBX_ << std::endl;    
      std::cout << "===> timeSecondMax = " << calculatedExtrahit_.timeSecondMax << "\tcalculatedExtrahit_.amplitudeSecondMax = " << calculatedExtrahit_.amplitudeSecondMax << std::endl;
      std::cout << "===> # hits with consistent time of the second max = " << nConsistentBXSecond_ << std::endl;    
  }

  std::cout << "listing the samples before any subtraction... " << std::endl;
  for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
    std::cout << "\tisample " << iSample << " has amplitudes_[iSample] = " << amplitudes_[iSample] << std::endl;
  }
  */

  for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
    extraHitFrame[iSample] = .0;

    // first subtract (for EE) an offset extrapolating a far OOT tail
    if(farTail_) {
      extraHitFrame[iSample] = std::max(0.0,farTailFunction(iSample));
      amplitudes_[iSample] -= extraHitFrame[iSample];
      calculatedExtrahit_.puTagged = true;
      //  std::cout << "\t\t\tFAR TAIL SUBTRACTING TO SAMPLE " << iSample << " extraHitFrame[iSample] = " << extraHitFrame[iSample] << std::endl;
    } 
    if(calculatedExtrahit_.amplitudeMax >= minSeedAmplitudeOutOfTime_ && nConsistentBX_ >= minNConsistentBX_ ) {
      // std::cout << "\tbefore pulse shape sub: isample " << iSample << " has amplitudes_[iSample] = " << amplitudes_[iSample] << std::endl;
      if(!farTail_ || calculatedExtrahit_.timeMax>0) extraHitFrame[iSample] += (iSample<3) ? 0.0 : std::max(0.0, pulseShapeFunction(iSample));
      calculatedExtrahit_.puTagged = true; 
    }
    
    //     if(calculatedExtrahit_.puTagged) std::cout << "\t\t\tFINAL SUBTRACTING TO SAMPLE " << iSample << " extraHitFrame[iSample] = " << extraHitFrame[iSample] 
    //                                               << "\t(farTail = " << farTail_ << ")" << std::endl;
  }
}


#endif
