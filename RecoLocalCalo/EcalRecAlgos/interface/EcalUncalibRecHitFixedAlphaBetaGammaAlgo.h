#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitFixedAlphaBetaGammaAlgo_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitFixedAlphaBetaGammaAlgo_HH

/** \class EcalUncalibRecHitFixedAlphaBetaGammaAlgo
  *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
  *  using an analytical function fit, with the pulse parameters alpha and beta fixed.
  *  It follows a fast fit algorithms devolped on test beam data by P. Jarry
  *  If the pedestal is not given, it is calculated from the first 2 pre-samples; 
  *  FIXME: conversion gainID (1,2,3) with gain factor (12,6,1) is hardcoded here !  
  *
  *  \author E. Di Marco
  */

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
//#include "CLHEP/Matrix/Matrix.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <string>

//#include "TROOT.h"
//#include "TMinuit.h"
//#include "TGraph.h"
//#include "TF1.h"
//#include "TMatrixD.h"
//#include "TVectorD.h"

template<class C> class EcalUncalibRecHitFixedAlphaBetaGammaAlgo : public EcalUncalibRecHitRecAbsAlgo<C>
{


 private:
  double MinAmpl_;
  bool dyn_pedestal;
  double fAlpha_;//parameter of the shape
  double fBeta_;//parameter of the shape
  double fGamma_;//parameter of the oot average shape
  double fAmp_max_;// peak amplitude 
  double fTim_max_ ;// time of the peak (in 25ns units)
  double fAmpOut_; // oot amplitude
  double alfabeta_; 
  double chi2P_;
  
  int fNb_iter_;
  int  fNum_samp_bef_max_ ;
  int fNum_samp_after_max_;
  
  float fSigma_ped; 
  double un_sur_sigma;
  //temporary solution for different alpha and beta
  float alpha_table_[36][1701];
  float beta_table_[36][1701];

  bool doFit_;

  double pulseShapeFunction(double t);
  float PerformAnalyticFit(double* samples, int max_sample);
  void InitFitParameters(double* samples, int max_sample);
  CLHEP::HepSymMatrix DM1_ ; CLHEP::HepVector temp_;
   public:

  EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>():fAlpha_(0.),fBeta_(0.),fGamma_(0.),fAmp_max_(-1.),fTim_max_(-1),fAmpOut_(-1.),alfabeta_(0),chi2P_(-1.),fNb_iter_(4),fNum_samp_bef_max_(5),fNum_samp_after_max_(3),fSigma_ped(1.1),DM1_(4),temp_(4){
    un_sur_sigma = 1./double(fSigma_ped) ;
    for (int i=0;i<36;i++){
      for(int j=0;j<1701;j++){
	alpha_table_[i][j] = 1.2 ;
	beta_table_[i][j]  =  1.7 ;
      }
    }
    doFit_ = false;
    MinAmpl_ = 16;
    dyn_pedestal = true;
  }
    EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>(int n_iter, int n_bef_max =5, int n_aft_max =3, float sigma_ped = 1.1):fAlpha_(0.),fBeta_(0.),fGamma_(0.),fAmp_max_(-1.),fTim_max_(-1),fAmpOut_(-1.),alfabeta_(0),chi2P_(-1.),DM1_(4),temp_(4){
    
   fNb_iter_ = n_iter ;
   fNum_samp_bef_max_   = n_bef_max ;
   fNum_samp_after_max_ = n_aft_max ;
   fSigma_ped = sigma_ped;
   un_sur_sigma = 1./double(fSigma_ped) ;
   for (int i=0;i<36;i++){
     for(int j=0;j<1701;j++){
       alpha_table_[i][j] = 1.2 ;
       beta_table_[i][j]  =  1.7 ;
     }
   }
   doFit_ = false;
   MinAmpl_ = 16;
   dyn_pedestal = true;
  };

  virtual ~EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>() { };
  virtual EcalUncalibratedRecHit makeRecHit(const C& dataFrame, const double* pedestals,
					    const double* gainRatios,
					    const EcalWeightSet::EcalWeightMatrix** weights, 
					    const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix); 
  void SetAlphaBeta( double alpha, double beta);
  void SetMinAmpl(double ampl);
  void SetPedestalRMS(double rms) { fSigma_ped = rms; }
  void SetDynamicPedestal(bool dyn_pede);
};



  /// Compute parameters
template<class C> EcalUncalibratedRecHit  EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>::makeRecHit(const C& dataFrame, const double* pedestals,
											     const double* gainRatios,
											     const EcalWeightSet::EcalWeightMatrix** weights, 
											     const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix){
  double chi2_(-1.);
  
  //  double Gain12Equivalent[4]={0,1,2,12};
  double frame[C::MAXSAMPLES];// will contain the ADC values
  double pedestal =0;     // carries pedestal for gain12 i.e. gainId==1
  
  int gainId0 = 1;        // expected gainId at the beginning of dataFrame 
  int iGainSwitch = 0;    // flags whether there's any gainId other than gainId0
  int GainId= 0;          // stores gainId at every sample
  double maxsample(-1);   // ADC value of maximal ped-subtracted sample
  int imax(-1);           // sample number of maximal ped-subtracted sample
  bool external_pede = false;
  bool isSaturated = 0;   // flag reporting whether gain0 has been found

  // Get time samples checking for Gain Switch and pedestals
  if(pedestals){
    external_pede = true;
    if(dyn_pedestal) { pedestal = (double(dataFrame.sample(0).adc()) + double(dataFrame.sample(1).adc()))/2.;}
    else{ pedestal  = pedestals[0];}
    for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
	//create frame in adc gain 12 equivalent
	GainId = dataFrame.sample(iSample).gainId();

	// FIX-ME: warning: the vector pedestal is supposed to have in the order G12, G6 and G1
        // if GainId is zero treat it as 3 temporarily to protect against undefined
	// frame will be set to ~max of gain1
	if ( GainId == 0 )
	  { 
	    GainId = 3;
	    isSaturated = 1;
	  }

	if (GainId != gainId0) iGainSwitch = 1;

	if(GainId==gainId0){frame[iSample] = double(dataFrame.sample(iSample).adc())-pedestal ;}
	else {frame[iSample] = (double(dataFrame.sample(iSample).adc())-pedestals[GainId-1])*gainRatios[GainId-1];}

	if( frame[iSample]>maxsample ) {
          maxsample = frame[iSample];
          imax = iSample;
	}
      }
  }
  else {// pedestal from pre-sample
    external_pede = false;
    pedestal = (double(dataFrame.sample(0).adc()) + double(dataFrame.sample(1).adc()))/2.;

    for(int iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
      //create frame in adc gain 12 equivalent
      GainId = dataFrame.sample(iSample).gainId();
      //no gain switch forseen if there is no external pedestal
      if ( GainId == 0 ) 
	{
	  GainId = 3;
	  isSaturated = 1;
	}

      frame[iSample] = double(dataFrame.sample(iSample).adc())-pedestal ;
      // if gain has switched but no pedestals are available, no much good you can do...
      if (GainId > gainId0) iGainSwitch = 1;
      if( frame[iSample]>maxsample ) {
	maxsample = frame[iSample];
	imax = iSample;
      }
    } 
  }

  if( (iGainSwitch==1 && external_pede==false) ||  // ... thus you return dummy rechit
      imax ==-1 ){                                 // protect against all frames being <-1
    return EcalUncalibratedRecHit( dataFrame.id(), -1., -100., -1. , -1.);
  }
  
  InitFitParameters(frame, imax);
  chi2_ = PerformAnalyticFit(frame,imax);
  uint32_t flags = 0;
  if (isSaturated) flags = EcalUncalibratedRecHit::kSaturated;

  /*
  std::cout << "separate fits\nA: " << fAmp_max_  
            <<", pedestal: "<<pedestal << ", tPeak " << fTim_max_ 
            <<", pileup ampli: " << fAmpOut_ << ", chi2 pileup:" << chi2P_ << std::endl;
  */
  
  EcalUncalibratedRecHit theHit( dataFrame.id(),fAmp_max_, pedestal, fTim_max_ - 5 , chi2_, flags );
  theHit.setOutOfTimeEnergy(fAmpOut_);
  theHit.setOutOfTimeChi2(chi2P_);
  return theHit;
}

template<class C> double EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>::pulseShapeFunction(double t){
  if( alfabeta_ <= 0 ) return((double)0.);
  double dtsbeta,variable,puiss;
  double dt = t-fTim_max_ ;
  double ootdecay=exp(-t/fGamma_);
  if(dt > -alfabeta_)  {
    dtsbeta=dt/fBeta_ ;
    variable=1.+dt/alfabeta_ ;
    puiss=pow(variable,fAlpha_);
    return fAmp_max_*puiss*exp(-dtsbeta)+fAmpOut_*ootdecay ;
  }
  return  fAmpOut_*ootdecay;
}

template<class C> void  EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>::InitFitParameters(double* samples, int max_sample){

  // in a first attempt just use the value of the maximum sample 
  fAmp_max_ = samples[max_sample];
  fTim_max_ = max_sample;
  fAmpOut_  = 1.; // in 2012 in EB ~40 MeV/channel
  fGamma_   = 1.; // almost flat in the first 10 samples

  // amplitude too low for fit to converge
  // timing set correctly is assumed
  // oot amplitude from the average of the pre-samples
  if(fAmp_max_ <  MinAmpl_){
    fAmp_max_      = samples[5];
    double sumA    = samples[5]+samples[4]+samples[6];
    if(sumA != 0) { fTim_max_ = 5+(samples[6]-samples[4])/sumA; }
    else{ fTim_max_ = -993; }//-999+6
    fAmpOut_       = (samples[0]+samples[1])/2.;
    fGamma_        = 100.; // flat
    doFit_ = false;
  }
  // if timing very badly off, that just use max sample
  else if(max_sample <1 || max_sample > 7)
    {    doFit_ = false;}
  else{
    //y=a*(x-xM)^2+b*(x-xM)+c
    doFit_ = true;
    float a = float(samples[max_sample-1]+samples[max_sample+1]-2*samples[max_sample])/2.;
    if(a==0){doFit_ =false; return;}
    float b = float(samples[max_sample+1]-samples[max_sample-1])/2.;
  
    fTim_max_ = max_sample - b/(2*a);
    fAmp_max_ =  samples[max_sample] - b*b/(4*a); 
    fAmpOut_  = (samples[0]+samples[1])/2.;
  }  
  
} 

template<class C> float EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>::PerformAnalyticFit(double* samples, int max_sample){
  
  //int fValue_tim_max = max_sample;  
  //! fit electronic function from simulation
  //! parameters fAlpha_ and fBeta_ are fixed and fit is providing the 5 following parameters
  //! the maximum amplitude ( fAmp_max_ ) 
  //! the time of the maximum  ( fTim_max_)
  //! the out-of-time pileup amplitude (fAmpOut_) 	
  //! the left tail parameter (fGamma_)

  double chi2=-1 , chi2P=-1 , db[4] ;
  

  //HepSymMatrix DM1(5) ; CLHEP::HepVector temp(5) ;
  
  int num_fit_min =(int)(max_sample - fNum_samp_bef_max_ ) ;
  int num_fit_max =(int)(max_sample + fNum_samp_after_max_) ;
  
  if (num_fit_min<0) num_fit_min = 0 ; 
  //if (num_fit_max>=fNsamples-1) num_fit_max = fNsamples-2 ;
  if (num_fit_max>= C::MAXSAMPLES ) {num_fit_max = C::MAXSAMPLES-1 ;}
  
  if(! doFit_ ) {
    LogDebug("EcalUncalibRecHitFixedAlphaBetaGammaAlgo")<<"No fit performed. The amplitude of sample 5 will be used";
    return -1;
  }
  
  double func,delta,deltaP ;
  double variation_func_max = 0. ;double variation_tim_max = 0. ;
  double variation_ampo = 0. ; double variation_gamma = 0. ;
  //!          Loop on iterations
  for (int iter=0 ; iter < fNb_iter_ ; iter ++) {
    //!          initialization inside iteration loop !
    chi2 = 0. ; //PROD.Zero() ;  DM1.Zero() ;
    chi2P = 0. ;

    for(int i1=0 ; i1<4 ; i1++) {
      temp_[i1]=0;
      for(int j1=i1 ; j1<4 ; j1++) { 
        DM1_.fast(j1+1,i1+1) = 0; }
    }

    fAmp_max_ += variation_func_max ;
    fTim_max_ += variation_tim_max ;
    fAmpOut_  += variation_ampo ;
    //! limit the largely negative oot aplitude
    if(fAmpOut_< -5*fSigma_ped) fAmpOut_ = -5*fSigma_ped;
    fGamma_   += variation_gamma ;
    //! limit the exponent decay (assume only < 0 OOT PU)
    if(fGamma_<0.5) fGamma_ = 0.5;
    //! limit gamma to 1000 (flat distribution already)
    if(fGamma_>1000.) fGamma_ = 1000.;
    
    //! Then we loop on samples to be fitted
    for( int i = num_fit_min ; i <= num_fit_max ; i++) {
      //if(i>fsamp_edge_fit && i<num_fit_min) continue ; // remove front edge samples
      //! calculate function to be fitted
      func = pulseShapeFunction( (double)i  ) ;
      // std::cout << "\t\tsample " << i << "\tfunc = " << func << std::endl; // debug
      //! then calculate derivatives of function to be fitted
      double dt =(double)i - fTim_max_ ;
      if(dt > -alfabeta_)  {      
	double dt_sur_beta = dt/fBeta_ ;
	double variable = (double)1. + dt/alfabeta_ ;
	double expo = exp(-dt_sur_beta) ;	 
	double puissance = pow(variable,fAlpha_) ;
	
	db[0]=un_sur_sigma*puissance*expo ;
	db[1]=fAmp_max_*db[0]*dt_sur_beta/(alfabeta_*variable) ;
      }
      else {
	db[0]=0. ; db[1]=0. ; 
      }
      if(i<3) {
        double t(i);
        double ootdecay = exp(-t/fGamma_);
        db[2]=ootdecay;
        db[3]=t/fGamma_/fGamma_*ootdecay;
      }
      else {
        db[2]=0. ; db[3]=0. ;
      }
      //! compute matrix elements DM1
      for(int i1=0 ; i1<4 ; i1++) {
	for(int j1=i1 ; j1<4 ; j1++) { 
	  //double & fast(int row, int col);
	  DM1_.fast(j1+1,i1+1) += db[i1]*db[j1]; }
      }
      //! compute delta
      delta = (samples[i]-func)*un_sur_sigma ;
      //! compute vector elements PROD
      for(int ii=0 ; ii<4 ;ii++) {temp_[ii] += delta*db[ii] ;}
      chi2 += delta * delta ;
      //! compute the chi2 for OOT (inverse definition, check consistency with flat pedestal)
      if(i<3) {
        deltaP = samples[i]*un_sur_sigma ;
        chi2P += deltaP * deltaP;
      }
    }//! end of loop on samples 
    
    // debug
    // std::cout << "Iter = " << iter << " Chi2 = " << chi2 << std::endl;
    // std::cout << "fGamma_ = " << fGamma_ << " fAmpOut_ = " << fAmpOut_ << std::endl;

    int fail=0 ;
    DM1_.invert(fail) ;
    if(fail != 0.) {
      //just a guess from the value of the parameters in the previous interaction;
      //printf("wH4PulseFitWithFunction =====> determinant error --> No Fit Provided !\n") ;
      InitFitParameters(samples,max_sample);
      return -101 ;
    }
    /*     for(int i1=0 ; i1<5 ; i1++) { */
    /*       for(int j1=0 ; j1<5 ; j1++) {  */
    /* 	//double & fast(int row, int col); */
    /* 	std::cout<<"inverted: "<<DM1[j1][i1]<<std::endl;;} */
    /*     } */
    /*     std::cout<<"vector temp: "<< temp[0]<<" "<<temp[1]<<" "<<temp[2]<<std::endl; */
    //! compute variations of parameters fAmp_max and fTim_max 
    CLHEP::HepVector PROD = DM1_*temp_ ;
    //    std::cout<<"vector PROD: "<< PROD[0]<<" "<<PROD[1]<<" "<<PROD[2]<<std::endl;

    // Probably the fastest way to protect against
    // +-inf value in the matrix DM1_ after inversion
    // (which is nevertheless flagged as successfull...)
    
    if ( std::isnan( PROD[0] ) ) {
      InitFitParameters(samples,max_sample);
      return -103 ;
    }

    variation_func_max = PROD[0] ;
    variation_tim_max = PROD[1] ;
    variation_ampo    = PROD[2] ;
    variation_gamma   = PROD[3] ;
    //chi2 = chi2/((double)nsamp_used - 3.) ;
  }//!end of loop on iterations       
  
  //!   protection again diverging/unstable fit 
  if( variation_func_max > 2000. || variation_func_max < -1000. ) {
    InitFitParameters(samples,max_sample);
    return -102;
  }
  
  //!      results of the fit are calculated
  fAmp_max_ += variation_func_max ;
  fTim_max_ += variation_tim_max ;
  fGamma_   += variation_gamma ;
  fAmpOut_  += variation_ampo ; // normalization of the full tail

  //! limit the largely negative oot aplitude
  if(fAmpOut_< -5*fSigma_ped) fAmpOut_ = -5*fSigma_ped;
  fGamma_   += variation_gamma ;
  //! limit the exponent decay (assume only < 0 OOT PU)
  if(fGamma_<0.5) fGamma_ = 0.5;
  //! limit gamma to 1000 (flat distribution already)
  if(fGamma_>1000.) fGamma_ = 1000.;

  // debug
  // std::cout << "hit fAmp_max_ = " << fAmp_max_ << "\tfTim_max_= " << fTim_max_ 
  //          << "\tfAmpOut_ = " << fAmpOut_ << "\tfGamma_ = " << fGamma_ << std::endl;
 
  // protection against unphysical results:
  // ampli mismatched to MaxSample, ampli largely negative, time off preselected range
  if( fAmp_max_>2*samples[max_sample]  ||  fAmp_max_<-100 ||  fTim_max_<0  ||  9<fTim_max_ ) {
    InitFitParameters(samples,max_sample);
    return -104;
  }
  
  chi2P_ = chi2P;
  //std::cout <<"chi2: "<<chi2<<" ampl: "<<fAmp_max_<<" time: "<<fTim_max_<<std::endl;
  return chi2;
}

template<class C> void EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>::SetAlphaBeta( double alpha, double beta){
  fAlpha_ = alpha;
  fBeta_=  beta;
  alfabeta_ = fAlpha_*fBeta_;
}

template<class C> void EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>::SetMinAmpl( double ampl){
  MinAmpl_ = ampl;
}
template<class C> void EcalUncalibRecHitFixedAlphaBetaGammaAlgo<C>::SetDynamicPedestal (bool p){
  dyn_pedestal = p;
}
#endif
