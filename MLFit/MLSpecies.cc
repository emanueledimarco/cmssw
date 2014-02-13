/**********************************************************************
 * MLSpecies
 * Nick Danielson / Amir Farbin
 *********************************************************************/

#include "TString.h"
#include "MLFit/MLSpecies.hh"
#include "RooAbsCategoryLValue.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooSuperCategory.h"

ClassImp(MLSpecies);

MLSpecies::MLSpecies(const char* name, const char* title) :
  TNamed(name, title),
  _yieldVar(0),
  _fracVar(0),
  _formula("")
{
}


MLSpecies::~MLSpecies()
{
}

RooAbsReal* MLSpecies::getFracVar()
{
  if (_fracVar == 0) _fracVar = new RooRealVar(TString("e_")+GetName(),TString("Epsilon ")+GetTitle(),1);
  return _fracVar;
}

RooAbsReal* MLSpecies::buildYieldVar(const RooArgList &args)
{
  if (args.getSize() != _protoParams.getSize()) {
    cout << GetName() << " Error: Trying to build yield var with " << args.getSize() << " when " << _protoParams.getSize() << " are expected:" << endl;
    args.Print();
    return 0;
  }
  
  if (_formula == "") return new RooRealVar(TString("N_")+GetName(),TString("N ")+GetTitle(),0);

  return new RooFormulaVar(TString("N_")+GetName(),TString("N ")+GetTitle(),_formula, args);
}

void MLSpecies::setYieldVar(const char* formula, RooArgList args) 
{
  _formula = formula;
  _protoParams.add(args);
}

RooArgSet MLSpecies::getProtoParams() 
{
  return _protoParams;
}

Bool_t MLSpecies::replace (TString oldparname, RooAbsArg &newpar, Bool_t verbose)
{
  Bool_t replaced(kFALSE);
  RooAbsArg *oldpar = _protoParams.find(oldparname);
  if (oldpar != NULL) {
    _protoParams.replace(*oldpar, newpar);
    replaced = kTRUE;
  }
  
  if (verbose && !replaced) {
    cout << "MLSpecies " << GetName() << ": Can't replace variable " << oldparname << " because it's not a parameter of this species." << endl;
  }
  return replaced;
  
}
