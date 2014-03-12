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
#include "RooNLLVar.h"
#include "RooMinuit.h"

template<class C> class EcalUncalibRecHitRecAnalFitAlgo : public EcalUncalibRecHitRecAbsAlgo<C>
{


 private:
  double MinAmpl_;
  bool dyn_pedestal;
  TH1D *shape_;
  std::vector<TH1D*> oot_shapes_;
  bool savePlot_;

 public:
  // constructor
  EcalUncalibRecHitRecAnalFitAlgo<C>(){
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
  void SetInTimeShape(TH1D* shape);
  void SetOutOfTimeShapes(std::vector<TH1D*> shapes);
  void SavePlot(bool saveplot) { savePlot_ = saveplot; }

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

  // amplitude too low for fit to converge 
  // timing set correctly is assumed 
  if(maxsample < MinAmpl_) {
    amplitude_ = frame[5];
    double sumA    = frame[5]+frame[4]+frame[6];
    if(sumA != 0) { jitter_ = 5+(frame[6]-frame[4])/sumA; }
    else{ jitter_ = -999; }
    pedestal_  = pedestal;
    chi2_ = -100.;
    return EcalUncalibratedRecHit( dataFrame.id(), amplitude_, pedestal_, jitter_ - 6, chi2_, flag);
  }

  // if timing very badly off, that just use max sample
  else if(imax <1 || imax > 7) {    
    amplitude_ = maxsample;
    pedestal_  = pedestal; 
    jitter_ = imax;
    chi2_ = -200.;  
    return EcalUncalibratedRecHit( dataFrame.id(), amplitude_, pedestal_, jitter_ - 6, chi2_, flag);
  }

  else {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooRealVar time("time","time",0,10);
    RooArgSet *varSet = new RooArgSet(time);
    RooDataHist theData("data","data",*varSet,&histo); 
    
    // --- real pulse
    RooRealVar meanGauss("meanGauss","meanGauss", 0, -2.5, 2.5);
    RooRealVar sigmaGauss("sigmaGauss","sigmaGauss", 0.1, 0, 1);
    ECALShapeConvGaussian pulse("0T",time,&(*shape_),&meanGauss,&sigmaGauss);

    // --- OOT -1 bx pulse
    RooDataHist shape1m("shape1m","",*varSet,oot_shapes_[9]);
    RooHistPdf pulse1m("1mT","",time,shape1m,8);

    // --- OOT -2 bx pulse
    RooDataHist shape2m("shape2m","",*varSet,oot_shapes_[8]);
    RooHistPdf pulse2m("2mT","",time,shape2m,8);

    // --- OOT -3 bx pulse
    RooDataHist shape3m("shape3m","",*varSet,oot_shapes_[7]);
    RooHistPdf pulse3m("3mT","",time,shape3m,8);

    RooRealVar N0("N0","N0",0,0,5e+6);
    RooExtendPdf extSigPulse("extSigPulse","extSigPulse",*(pulse.model),N0);
  
    RooRealVar N1m("N1m","N1m",10,0,10000);
    RooExtendPdf ext1mPulse("ext1mPulse","ext1mPulse",pulse1m,N1m);

    RooRealVar N2m("N2m","N2m",1,0,10000);
    RooExtendPdf ext2mPulse("ext2mPulse","ext2mPulse",pulse2m,N2m);

    RooRealVar N3m("N3m","N3m",1,0,10000);
    RooExtendPdf ext3mPulse("ext3mPulse","ext3mPulse",pulse3m,N3m);
    
    RooAddPdf model("model","model",RooArgList(extSigPulse,ext1mPulse,ext2mPulse,ext3mPulse));

    RooFitResult *fitResult=0;
    fitResult = model.fitTo(theData,
                            RooFit::Extended(),
                            RooFit::Strategy(0),
                            RooFit::Minimizer("Minuit2","migrad"),
                            RooFit::Hesse(kFALSE),
                            RooFit::Verbose(kFALSE),
                            RooFit::PrintLevel(-1),
                            RooFit::Warnings(kFALSE),
                            RooFit::Save());

    if(fitResult->covQual()>1) {

      // --- get the parameters of interest
      //  std::cout << "=== Resulting parameters: ====" << std::endl;
      RooRealVar *timeBias = (RooRealVar*)fitResult->floatParsFinal().find("meanGauss");
      float timeMax = timeBias->getVal();
      //    float timeMaxErr = timeBias->getError();
      //  std::cout << "time bias = " << timeMax * 25. << " +/- " << timeMaxErr * 25. << " ns." << std::endl;
    
      RooArgSet *obs = new RooArgSet(time);
      float val = extSigPulse.getVal(obs);
      RooRealVar *N0fit = (RooRealVar*)fitResult->floatParsFinal().find("N0");
      float norm = N0fit->getVal();
      float amplitude = val * norm;
      //  std::cout << "==> amplitude = " << amplitude << std::endl;

      amplitude_ = amplitude;
      pedestal_  = pedestal;
      jitter_    = timeMax;
    
      RooPlot *plot = time.frame(C::MAXSAMPLES);
      theData.plotOn(plot,RooFit::Name("data"));
      model.plotOn(plot,RooFit::Name("model"));
      model.plotOn(plot,RooFit::Components("extSigPulse"),RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
      model.plotOn(plot,RooFit::Components("ext1mPulse"),RooFit::LineStyle(kDashed),RooFit::LineColor(kBlue+2));
      model.plotOn(plot,RooFit::Components("ext2mPulse"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen+2));
      model.plotOn(plot,RooFit::Components("ext3mPulse"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed+2));
    
      if(savePlot_) {
        TCanvas c1("c1","c1",600,600);
        plot->Draw();
        c1.SaveAs("fit.pdf");
      }
      int ndof=6;
      chi2_ = plot->chiSquare("model","data",ndof);
      if (isSaturated) flag = EcalUncalibratedRecHit::kSaturated;
      delete obs;
    } else {
      amplitude_ = frame[5];
      double sumA    = frame[5]+frame[4]+frame[6];
      if(sumA != 0) { jitter_ = 5+(frame[6]-frame[4])/sumA; }
      else{ jitter_ = -999; }
      pedestal_  = pedestal;
      chi2_ = -1.;
      return EcalUncalibratedRecHit( dataFrame.id(), amplitude_, pedestal_, jitter_ - 6, chi2_, flag);
    }
    delete varSet;
    return EcalUncalibratedRecHit( dataFrame.id(), amplitude_, pedestal_, jitter_ - 6, chi2_, flag);
  }
}

template<class C> void EcalUncalibRecHitRecAnalFitAlgo<C>::SetMinAmpl( double ampl){
  MinAmpl_ = ampl;
}
template<class C> void EcalUncalibRecHitRecAnalFitAlgo<C>::SetDynamicPedestal (bool p){
  dyn_pedestal = p;
}
template<class C> void EcalUncalibRecHitRecAnalFitAlgo<C>::SetInTimeShape(TH1D* shape){
  shape_ = shape;
}
template<class C> void EcalUncalibRecHitRecAnalFitAlgo<C>::SetOutOfTimeShapes(std::vector<TH1D*> shapes){
  for(unsigned int i=0; i<shapes.size(); ++i) oot_shapes_.push_back(shapes[i]);
}

#endif
