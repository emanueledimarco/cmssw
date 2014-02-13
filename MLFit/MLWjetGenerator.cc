#include "MLFit/MLWjetGenerator.hh"
#include "MLFit/MLSpecies.hh"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooRandom.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "TFile.h"

ClassImp(MLWjetGenerator);

MLWjetGenerator::MLWjetGenerator(const MLFit &theFit, const char* model, Bool_t doBVETO, Bool_t doCHARGE, Bool_t doTruth) :
  MLGenerator(theFit, model, doTruth)
{
  _doBVETO = doBVETO;
  _doCHARGE = doCHARGE;

}


RooDataSet* MLWjetGenerator::generatePrototype(const MLSpecies &spec, RooArgSet &genVars, Int_t ngen, Bool_t verbose)
{
  RooDataSet *cData ;

  // Prototype variables to generate:
  // BVETO

  TString speciesname = spec.GetName();
  
  //===================================
  // Charge ONLY
  //==================================
  if(!_doBVETO && _doCHARGE) {
    RooRealVar *BVETO       = _theFit->RealObservable("BVETO"); 
    cData = new RooDataSet("cData", "BVeto Data", *BVETO);
    if (speciesname == "sig_a" || speciesname == "ttbar_a" || speciesname == "bkg_a") {
      BVETO->setVal(+1.);
    }else if(speciesname == "sig_r" || speciesname == "ttbar_r" || speciesname == "bkg_r") {
      BVETO->setVal(-1.);
    }
    for (Int_t i = 0; i < ngen; i++) cData->add(*BVETO);
  }

  //===================================
  // BVETO ONLY
  //==================================

  else if(_doBVETO && !_doCHARGE) {
    RooRealVar *charge       = _theFit->RealObservable("charge"); 
    cData = new RooDataSet("cData", "charge Data", *charge);
    if (speciesname == "W+" || speciesname == "qcd+" || speciesname == "other+") {
      charge->setVal(+1.);
    }else if(speciesname == "W-" || speciesname == "qcd-" || speciesname == "other-") {
      charge->setVal(-1.);
    }
    for (Int_t i = 0; i < ngen; i++) cData->add(*charge);
  }

  //===================================
  // BVETO AND CHARGE
  //==================================
  else if(_doBVETO && !_doCHARGE) {
    // TO IMPLEMENT
  }
  return cData;
}

