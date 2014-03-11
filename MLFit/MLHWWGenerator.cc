#include "MLFit/MLHWWGenerator.hh"
#include "MLFit/MLSpecies.hh"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooRandom.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "TFile.h"

ClassImp(MLHWWGenerator);

MLHWWGenerator::MLHWWGenerator(const MLFit &theFit, const char* model, Bool_t doTruth) :
  MLGenerator(theFit, model, doTruth)
{}


RooDataSet* MLHWWGenerator::generatePrototype(const MLSpecies &spec, RooArgSet &genVars, Int_t ngen, Bool_t verbose)
{
  RooDataSet *cData ;

  // Prototype variables to generate:
 
  TString speciesname = spec.GetName();
  
  //===================================
  // The sub-channel
  //==================================
  RooRealVar *finalstate       = _theFit->RealObservable("finalstate"); 
  cData = new RooDataSet("cData", "Finalstate Data", *finalstate);
  if(speciesname == "sig_ee" || speciesname == "WW_ee" || speciesname == "ttbar_ee" || speciesname == "other_ee") {
    finalstate->setVal(0);
  } else if(speciesname == "sig_mm" || speciesname == "WW_mm" || speciesname == "ttbar_mm" || speciesname == "other_mm") {
    finalstate->setVal(1);
  } else if(speciesname == "sig_em" || speciesname == "WW_em" || speciesname == "ttbar_em" || speciesname == "other_em") {
    finalstate->setVal(2);
  }
  for (Int_t i = 0; i < ngen; i++) cData->add(*finalstate);

  return cData;
}

