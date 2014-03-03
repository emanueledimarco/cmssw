#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitRecAnalFitAlgo_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitRecAnalFitAlgo_HH

/** \class EcalUncalibRecHitRecAnalFitAlgo
 *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
 *  using an analytical fit
 *
 *  $Id: EcalUncalibRecHitRecAnalFitAlgo.h,v 1.11 2009/03/27 18:07:38 ferriff Exp $
 *  $Date: 2009/03/27 18:07:38 $
 *  $Revision: 1.11 $
 *  \author E. Di Marco
 */

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/PulseShapeTemplate.hh"
#include <vector>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooPolynomial.h"
#include "RooPlot.h"

template<class C> class EcalUncalibRecHitRecAnalFitAlgo : public EcalUncalibRecHitRecAbsAlgo<C>
{


 private:
  double MinAmpl_;
  bool dyn_pedestal;
  bool doFit_;

 public:
  // constructor
  EcalUncalibRecHitRecAnalFitAlgo<C>(){
    doFit_ = false;
    MinAmpl_ = 16;
    dyn_pedestal = true;
  }
  // destructor
  virtual ~EcalUncalibRecHitRecAnalFitAlgo<C>() { };
  virtual EcalUncalibratedRecHit makeRecHit(const C& dataFrame, const double* pedestals,
					    const double* gainRatios,
					    const EcalWeightSet::EcalWeightMatrix** weights, 
					    const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix); 

  void SetMinAmpl(double ampl);
  void SetDynamicPedestal(bool dyn_pede);
};

/// Compute parameters
template<class C> EcalUncalibratedRecHit  EcalUncalibRecHitRecAnalFitAlgo<C>::makeRecHit(const C& dataFrame, const double* pedestals,
                                                                                         const double* gainRatios,
                                                                                         const EcalWeightSet::EcalWeightMatrix** weights, 
                                                                                         const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix) { 
  double amplitude_(-1.),  pedestal_(-1.), jitter_(-1.), chi2_(-1.);
  
  // Get time samples
  //HepMatrix frame(C::MAXSAMPLES, 1);
  double frame[C::MAXSAMPLES]; // will contain the ADC values
  double pedestal =0;     // carries pedestal for gain12 i.e. gainId==1 
  int gainId0 = 1;        // expected gainId at the beginning of dataFrame
  int iGainSwitch = 0;    // flags whether there's any gainId other than gainId0
  int GainId= 0;          // stores gainId at every sample 
  double maxsample(-1);   // ADC value of maximal ped-subtracted sample
  int imax(-1);           // sample number of maximal ped-subtracted sample   
  bool isSaturated = 0;   // flag reporting whether gain0 has been found
  bool external_pede = false;
  uint32_t flag = 0;
  TH1D histo("histo","",C::MAXSAMPLES,0,C::MAXSAMPLES);

  cout << "digis:" << endl;

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
      histo.SetBinContent(iSample+1,frame[iSample]);
      cout << "pedestal = " << pedestal 
           << "\tadc = " << dataFrame.sample(iSample).adc()
           << "\tframe[iSample] = " << frame[iSample] << endl;
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
      histo.SetBinContent(iSample+1,frame[iSample]);
    } 
  }
    
  if( (iGainSwitch==1 && external_pede==false) ||  // ... thus you return dummy rechit
      imax ==-1 ){                                 // protect against all frames being <-1
    return EcalUncalibratedRecHit( dataFrame.id(), -1., -100., -1. , -1.);
  }

  // Get the default ECAL templates
  TFile *fileTemplates = TFile::Open("RecoLocalCalo/EcalRecAlgos/data/EcalShapes.root");
  TH1D *shape = (TH1D*)fileTemplates->Get("EcalBarrelShape");

  RooRealVar time("time","time",0,10);
  RooDataHist theData("data","data",RooArgList(time),&histo); 
    
  //--- pedestal
/*   RooRealVar pedSlope("pedSlope","pedSlope",0); */
/*   RooPolynomial pedestalPoly("pedestal","pedestal",time,pedSlope); */
  
  // --- real pulse
  RooRealVar meanGauss("meanGauss","meanGauss", 0, -2.5, 2.5);
  RooRealVar sigmaGauss("sigmaGauss","sigmaGauss", 0.1, 0, 1);
  ECALShapeConvGaussian pulse("0T",time,&(*shape),&meanGauss,&sigmaGauss);

  // --- OOT -1 bx pulse
  RooRealVar meanGauss1m("meanGauss1m","meanGauss1m", -10);
  RooRealVar sigmaGauss1m("sigmaGauss1m","sigmaGauss1m", 1);
  ECALShapeConvGaussian pulse1m("1mT",time,&(*shape),&meanGauss1m,&sigmaGauss1m);
    
  RooRealVar N0("N0","N0",0,0,5e+6);
  RooExtendPdf extSigPulse("extSigPulse","extSigPulse",*(pulse.model),N0);
  
  RooRealVar N1m("N1m","N1m",0,0,10000);
  RooExtendPdf ext1mPulse("ext1mPulse","ext1mPulse",*(pulse1m.model),N1m);
    
/*   RooRealVar Nped("Nped","Nped",pedestal*C::MAXSAMPLES); */
/*   RooExtendPdf extPed("extPed","extPed",pedestalPoly,Nped); */
    
  RooAddPdf model("model","model",RooArgList(extSigPulse,ext1mPulse));
    
  RooFitResult *fitResult=0;
  fitResult = model.fitTo(theData,
                          RooFit::Extended(),
                          RooFit::Strategy(2),
                          RooFit::Save());
    
  RooPlot *plot = time.frame(C::MAXSAMPLES);
  theData.plotOn(plot);
  model.plotOn(plot);
  //  model.plotOn(plot,RooFit::Components("extPed"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
  model.plotOn(plot,RooFit::Components("ext1mPulse"),RooFit::LineStyle(kDashed),RooFit::LineColor(kBlue));
  model.plotOn(plot,RooFit::Components("extSigPulse"),RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
    
  TCanvas c1("c1","c1",600,600);
  plot->Draw();
  c1.SaveAs("fit.pdf");
    
  // --- get the parameters of interest
  std::cout << "=== Resulting parameters: ====" << std::endl;
  RooRealVar *timeBias = (RooRealVar*)fitResult->floatParsFinal().find("meanGauss");
  float timeMax = timeBias->getVal();
  float timeMaxErr = timeBias->getError();
  std::cout << "time bias = " << timeMax * 25. << " +/- " << timeMaxErr * 25. << " ns." << std::endl;
    
  RooArgSet obs(time);
  float val = extSigPulse.getVal(&obs);
  RooRealVar *N0fit = (RooRealVar*)fitResult->floatParsFinal().find("N0");
  float norm = N0fit->getVal();
  float amplitude = val * norm;
  std::cout << "==> amplitude = " << amplitude << std::endl;

/*   float pedval = extPed.getVal(&obs); */
/*   RooRealVar *Npedfit = (RooRealVar*)fitResult->floatParsFinal().find("Nped"); */
/*   float pedNorm = Npedfit->getVal(); */
/*   float pedestalVal = pedval * pedNorm; */
/*   std::cout << "==> pedestal = " << pedestalVal << std::endl; */

  //    if ( std::string(gMinuit->fCstatu.Data()) == std::string("CONVERGED ") ) {
  if(1) {

    amplitude_ = amplitude;
    pedestal_  = pedestal;
    jitter_    = timeMax;
    chi2_ = 1.; // successful fit
    if (isSaturated) flag = EcalUncalibratedRecHit::kSaturated;
    /*
      std::cout << "separate fits\nA: " <<  amplitude_value << ", Ped: " << pedestal_value
      << ", t0: " << jitter_ << ", tp: " << pulseShape.GetParameter(1)
      << ", alpha: " << pulseShape.GetParameter(2)
      << std::endl;
    */

  }

  return EcalUncalibratedRecHit( dataFrame.id(), amplitude_, pedestal_, jitter_ - 6, chi2_, flag);
}

template<class C> void EcalUncalibRecHitRecAnalFitAlgo<C>::SetMinAmpl( double ampl){
  MinAmpl_ = ampl;
}
template<class C> void EcalUncalibRecHitRecAnalFitAlgo<C>::SetDynamicPedestal (bool p){
  dyn_pedestal = p;
}

#endif
