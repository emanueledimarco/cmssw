
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooPolynomial.h"

#include "RecoLocalCalo/EcalRecAlgos/src/PulseShapeTemplate.cc"

void testTemplateFit() {

  // Get the default ECAL templates
  TFile *fileTemplates = TFile::Open("RecoLocalCalo/EcalRecAlgos/data/EcalShapes.root");
  TH1D *shape = (TH1D*)fileTemplates->Get("EcalBarrelShape");

  TFile *fileToy = TFile::Open("RecoLocalCalo/EcalRecAlgos/test/gensample.root");
  RooDataSet *theData = (RooDataSet*)fileToy->Get("theData");
  RooRealVar *time = new RooRealVar("time","time",0,10);


  //--- pedestal
  RooRealVar pedSlope("pedSlope","pedSlope",0);
  RooPolynomial *pedestal = new RooPolynomial("pedestal","pedestal",*time,pedSlope);
  
  // --- real pulse
  RooRealVar *meanGauss  = new RooRealVar("meanGauss","meanGauss", 0, -5, 5);
  RooRealVar *sigmaGauss = new RooRealVar("sigmaGauss","sigmaGauss", 1, 0, 5);
  ECALShapeConvGaussian *pulse = new ECALShapeConvGaussian("0T",*time,shape,meanGauss,sigmaGauss);

  // --- OOT -1 bx pulse
  RooRealVar *meanGauss1m  = new RooRealVar("meanGauss1m","meanGauss1m", -10, -15, -5);
  RooRealVar *sigmaGauss1m = new RooRealVar("sigmaGauss1m","sigmaGauss1m", 1, 0, 5);
  ECALShapeConvGaussian *pulse1m = new ECALShapeConvGaussian("1mT",*time,shape,meanGauss1m,sigmaGauss1m);

  RooRealVar N0("N0","N0",0,0,10000);
  RooExtendPdf *extSigPulse = new RooExtendPdf("extSigPulse","extSigPulse",*(pulse->model),N0);
  
  RooRealVar N1m("N1m","N1m",0,0,10000);
  RooExtendPdf *ext1mPulse = new RooExtendPdf("ext1mPulse","ext1mPulse",*(pulse1m->model),N1m);

  RooRealVar Nped("Nped","Nped",0,0,10000);
  RooExtendPdf *extPed = new RooExtendPdf("extPed","extPed",*pedestal,Nped);
  
  RooAddPdf *model = new RooAddPdf("model","model",RooArgList(*extSigPulse,*ext1mPulse,*extPed));
  
  RooFitResult *fitResult=0;
  fitResult = model->fitTo(*theData,
			   RooFit::Extended(),
			   RooFit::Strategy(2),
			   RooFit::Save());
  
}
