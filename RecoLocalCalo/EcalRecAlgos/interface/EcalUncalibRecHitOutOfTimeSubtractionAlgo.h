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
  void init( const C &dataFrame, const EcalSampleMask &sampleMask, const double * pedestals, const double * pedestalRMSes, const double * gainRatios, std::vector<C> neighbors, double Nconst, double Cconst );
  double pulseShapeFunction(double t);
  void computeAmplitudeOOT( std::vector< double > &amplitudeFitParameters, std::vector< double > &subtractionLimits, double time );
  CalculatedExtraHit getCalculatedExtraHit() { return calculatedExtrahit_; };
  
 protected:
  
  EcalSampleMask sampleMask_;
  DetId          theDetId_;
  std::vector < double > amplitudes_;
  double pedestal_;
  double alpha_, beta_, alphabeta_;
  double maxSampleOutOfTime_, minAmplitudeOutOfTime_;

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
void EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::init( const C &dataFrame, const EcalSampleMask &sampleMask,
                                                         const double * pedestals, const double * pedestalRMSes, const double * gainRatios, 
                                                         std::vector<C> neighbors, double Nconst, double Cconst ) {
    
  theDetId_ = DetId(dataFrame.id().rawId());  
  calculatedExtrahit_.timeMax = 5;
  calculatedExtrahit_.amplitudeMax = 0;
  calculatedExtrahit_.amplitudeExtapolated = 0;
  amplitudes_.clear();
  amplitudes_.reserve(C::MAXSAMPLES);
	
  // obtain the max sample in the 3x3 matrix. 
  // pedestal is super-simple: 200 ADC (we are interested only in the time of the max here)
  double ampliAverageNeighbors = 0.0;
  double timeAverageNum = 0;
  double timeAverageDen = 0;
  for(typename std::vector<C>::const_iterator itn=neighbors.begin(); itn!=neighbors.end(); ++itn) {
    double maxSampleThisNeighbor = 5;
    double maxAmpThisNeighbor = -999.;
    for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
      int gainId = itn->sample(iSample).gainId(); 

      int  sampleAdc = 0;
      if ( gainId == 1 ) sampleAdc = itn->sample(iSample).adc() - 200;
      else if ( gainId == 2 ) sampleAdc = (itn->sample(iSample).adc() - 200) * 2 ;
      else sampleAdc = (dataFrame.sample(iSample).adc() - 200) * 12 ;

      if( sampleAdc > maxAmpThisNeighbor ) {
        maxAmpThisNeighbor  = sampleAdc;
        maxSampleThisNeighbor = iSample;
      }
    }// loop on samples
    double cterm         = Cconst;
    double sigmaped      = pedestalRMSes[0];  // approx
    if(maxAmpThisNeighbor<1) maxAmpThisNeighbor = 1;
    double nterm         = Nconst*sigmaped/maxAmpThisNeighbor;
    double sigmat        = std::sqrt( nterm*nterm  + cterm*cterm   )/25.0;
    timeAverageNum      += maxSampleThisNeighbor/sigmat/sigmat;
    timeAverageDen      += 1.0/sigmat/sigmat;
    ampliAverageNeighbors += maxAmpThisNeighbor;
  } // loop over neighbors
  
  if(neighbors.size()>0) {
    calculatedExtrahit_.amplitudeMax = ampliAverageNeighbors/double(neighbors.size());
    calculatedExtrahit_.timeMax = timeAverageNum/timeAverageDen;
  } else calculatedExtrahit_.amplitudeMax = -100.;

  // to obtain gain 12 pedestal:
  // -> if it's in gain 12, use first sample
  // --> average it with second sample if in gain 12 and 3-sigma-noise compatible (better LF noise cancellation)
  // -> else use pedestal from database
  pedestal_ = 0;
  int num_  = 0;
  if (dataFrame.sample(0).gainId() == 1 &&
      sampleMask_.useSample(0, theDetId_ ) 
      ) {
    pedestal_ += double (dataFrame.sample(0).adc());
    num_++;
  }
  if (num_!=0 &&
      dataFrame.sample(1).gainId() == 1 && 
      sampleMask_.useSample(1, theDetId_) &&
      fabs(dataFrame.sample(1).adc()-dataFrame.sample(0).adc())<3*pedestalRMSes[0]) {
    pedestal_ += double (dataFrame.sample(1).adc());
    num_++;
  }
  if (num_ != 0)
    pedestal_ /= num_;
  else
    pedestal_ = pedestals[0];
  
  // fill vector of amplitudes
  // ped-subtracted and gain-renormalized samples. It is VERY
  // IMPORTANT to have samples one clock apart which means to
  // have vector size equal to MAXSAMPLES
  double sample;
  int GainId;
  for (int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
    
    GainId = dataFrame.sample(iSample).gainId();
    
    // only use normally samples which are desired; if sample not to be used
    // inflate error so won't generate ratio considered for the measurement
    if (!sampleMask_.useSample(iSample, theDetId_ ) ) 
      sample      = 1e-9;
    else if (GainId == 1) 
      sample      = double (dataFrame.sample(iSample).adc() - pedestal_);
    else if (GainId == 2 || GainId == 3) 
      sample      = (double (dataFrame.sample(iSample).adc() - pedestals[GainId - 1])) *gainRatios[GainId - 1];
    else
      sample      = 1e-9;  // GainId=0 case falls here, from saturation

    amplitudes_.push_back(sample);
  }

}

template<class C>
void EcalUncalibRecHitOutOfTimeSubtractionAlgo<C>::computeAmplitudeOOT( std::vector< double > &amplitudeFitParameters, std::vector< double > &subtractionLimits, double time )
{

  // these are the parameters of the pulse shape function
  alphabeta_ = amplitudeFitParameters[0]*amplitudeFitParameters[1];
  alpha_ = amplitudeFitParameters[0];
  beta_ = amplitudeFitParameters[1];

  // these are the limits of the subtraction algorithm
  maxSampleOutOfTime_ = subtractionLimits[0];
  minAmplitudeOutOfTime_ = subtractionLimits[1];

  // if the max sample of the max of the neighbors is 
  // in the pre-samples, and the amplitude is above noise, 
  // --> there is a hope to fit for the out-of-time contributions
  // -> not, return 0
  double amplitudeExtapolated_= -1.0;
  if (calculatedExtrahit_.amplitudeMax == -100) amplitudeExtapolated_ = 0.0;
  else if (calculatedExtrahit_.timeMax > maxSampleOutOfTime_ ) amplitudeExtapolated_ = 0.0;
  else if (calculatedExtrahit_.amplitudeMax < minAmplitudeOutOfTime_) amplitudeExtapolated_ = 0.0;
  else {

    /*
    std::cout << "===> timeMax = " << calculatedExtrahit_.timeMax << "\tcalculatedExtrahit_.amplitudeMax = " << calculatedExtrahit_.amplitudeMax << std::endl;
    std::cout << "listing the samples ... " << std::endl;
    for (int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
      std::cout << "\tisample " << iSample << " has amplitudes_[iSample] = " << amplitudes_[iSample] << std::endl;
    }
    std::cout << "Done listing samples " << std::endl;
    */

    amplitudeExtapolated_ = pulseShapeFunction(time);
    //    std::cout << "==> AMPLITUDE EXTRAP = " << amplitudeExtapolated_ << std::endl;
  }

  calculatedExtrahit_.amplitudeExtapolated = amplitudeExtapolated_;
}


#endif
