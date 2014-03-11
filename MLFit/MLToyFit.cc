/**********************************
 *MLToyFit
 *Nick Danielson
 ***********************************/

#include "MLFit/MLToyFit.hh"
#include "RooDataSet.h"
#include "RooErrorVar.h"
#include "RooRealVar.h"

ClassImp(MLToyFit);

MLToyFit::MLToyFit(const RooAbsPdf &pdf, const RooArgSet &dependents, const RooArgSet &extraParams) :
  TNamed(pdf.GetName(),pdf.GetTitle()),
  _pdf((RooAbsPdf*)&pdf),
  _dependents(dependents),
  _datacuts("")
{
  // Get fit parameters...
  _fitParams = pdf.getParameters(&dependents);

  // Initial parameters...
  _fitInitParams = (RooArgSet*)_fitParams->snapshot(kTRUE);
 
  // Create extra variables you want stored...
  _nllVar = new RooRealVar("NLL","-log(likelihood)",0);
  _status = new RooRealVar("status", "Fit Status",0);
  _covQual = new RooRealVar("covQual", "Covariance Quality",0);

  _extraParams.add(*(RooArgSet*)extraParams.snapshot(kTRUE));

  // Create data set containing parameter values, errors, and extra variables...
  RooArgSet fitparams(*_fitParams);
  fitparams.add(RooArgSet(*_nllVar,*_status,*_covQual));
  fitparams.add(_extraParams);
  fitparams.setAttribAll("StoreError",kTRUE);
  fitparams.setAttribAll("StoreAsymError",kTRUE);
  _fitParData = new RooDataSet("fitParData", "Fit Parameters Dataset", fitparams);
  fitparams.setAttribAll("StoreError",kFALSE);
  fitparams.setAttribAll("StoreAsymError",kFALSE);
}

MLToyFit::~MLToyFit()
{
}

void MLToyFit::addErrVars()
{
  //Adds the fit errors to the fitpar dataset...
  
  TIterator *iter = _fitParams->createIterator();
  while (RooRealVar *par = (RooRealVar*)iter->Next()) {

    RooErrorVar* err = par->errorVar() ;
    _fitParData->addColumn(*err) ;
  }
  delete iter;
}

void MLToyFit::initializeParameters() 
{
  *_fitParams = *_fitInitParams;
  
  // And now make sure everything's constant that should be...
  TIterator *parIter = _fitParams->createIterator();
  while (RooRealVar* par = (RooRealVar*)parIter->Next()) {
    par->setAttribute("Constant",((RooRealVar*)_fitInitParams->find(par->GetName()))->getAttribute("Constant"));
  }
  
}
