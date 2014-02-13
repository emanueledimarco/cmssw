#include "MLFit/MLPdfRecipe.hh"
#include "MLFit/MLPdfFactory.hh"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooCatType.h"
#include "RooSimultaneous.h"
#include "RooStringVar.h"
#include "RooSuperCategory.h"

ClassImp(MLPdfRecipe);

// Normal constructor...
MLPdfRecipe::MLPdfRecipe(const char* name, const char* title, const char* type, const RooArgList &obs, const RooArgList &partial, const RooArgList &parameters, const TList &args) :
  TNamed(name, title),
  _type(type),
  _obs(obs),
  _partialList(partial),
  //  _parameters(parameters),
  _protoParams(parameters),
  _masterSplitCat(0),
  _simPdfCat(0),
  _rf(0)
{
  TIterator *argIter = args.MakeIterator();
  while(TObject *arg = argIter->Next()) _args.Add(arg);
  delete argIter;
}

// Sim PDF constructor
MLPdfRecipe::MLPdfRecipe(const char* name, const char* title, const char* obs, const RooAbsCategoryLValue &cat) :
  TNamed(name, title),
  _type("RooSimultaneous"),
  _obs(obs),
  _masterSplitCat(0),
  _simPdfCat((RooAbsCategoryLValue*)&cat),
  _rf(0)
{
}


MLPdfRecipe::~MLPdfRecipe()
{
}

Bool_t MLPdfRecipe::replace(const char* oldparname, RooAbsArg &newpar, Bool_t verbose)
{
  // Give the name of a parameter that this PDF depends on (before splitting), and I will replace it with the new one.
  // There are three kinds of paramters that a PDF can depend on:
  // * Normal parameter, like mean of a gaussian
  // * A parameter of the convolution function
  // * If the PDF is a simultaneous PDF, a parameter of one of its PDF components.

  Bool_t replaced(kFALSE);
  // Replace normal parameters
  RooAbsArg* oldpar = _protoParams.find(oldparname);
  if (oldpar != 0) {
    _protoParams.replace(*oldpar, newpar);
    replaced = kTRUE;
  }

  // Replace resolution function parameters
  if (isConvoluted()) {
    Bool_t r = _rf->replace(oldparname, newpar, verbose);
    replaced = replaced || r;
  }
  
  // Replace sim pdf parameters.
  if (isSimultaneous()) {
    TIterator *recIter = _simPdfRecipeList.MakeIterator();
    while (MLPdfRecipe *rec = (MLPdfRecipe*)recIter->Next()){
      Bool_t r = rec->replace(oldparname, newpar, verbose);
      replaced = replaced || r;
    }
    delete recIter;
  }
  
  if (verbose && !replaced) {
    cout << "MLPdfRecipe " << GetName() << ": Can't replace variable " << oldparname << " because it's not a parameter of this PDF." << endl;
  }
  
  return replaced;
}


RooArgSet MLPdfRecipe::getProtoParams()
{
  // First add the parameters that this PDF directly depends on...
  RooArgSet protoParams(_protoParams);

  // Now the parameters from its resolution function...
  if (isConvoluted()) protoParams.add(_rf->getProtoParams(),kTRUE);
  
  // Now the parameters from its component pdfs if it's a sim pdf...
  if (isSimultaneous()) {
    TIterator *recIter = _simPdfRecipeList.MakeIterator();
    while (MLPdfRecipe* rec = (MLPdfRecipe*)recIter->Next()) protoParams.add(rec->getProtoParams(),kTRUE);
    delete recIter;
  }
  
  return protoParams;
}

TString MLPdfRecipe::getSplitName() const
{
  if (!isSplit()) return GetName();
  //cout << "Getsplitname in " << GetName() << endl;
  //_masterSplitCat->Print();
  //cout << "done print" << endl;
  return TString(GetName()) + "_" + _masterSplitCat->getLabel();
}

MLPdfRecipe* MLPdfRecipe::getSimComponent() const 
{
  // If this is a recipe for a simultaneous pdf, this returns the recipe for the correct PDF based on the current
  // value of simPdfCat.

  if (!isSimultaneous()) {
    cout << GetName() << " isn't a simultaneous pdf, but getSimComponent() was called!" << endl;
    return 0;
  }
  
  RooStringVar *catLabel = (RooStringVar*)_simPdfLabels.find(_simPdfCat->getLabel());
  if (catLabel == 0) {  
    cout << GetName() << ": There is no label " << _simPdfCat->getLabel() << " in label list: " << endl;
    _simPdfLabels.Print();
    return 0;
  }
  
  MLPdfRecipe *rec = (MLPdfRecipe*)_simPdfRecipeList.At(_simPdfLabels.index(catLabel));
  if (rec == 0) {
    cout << GetName() << "::getSimComponent() Can't match label " << _simPdfCat->getLabel() << " with a pdf recipe" << endl;
    cout << "Simultaneous PDF recipe component list: " << endl;
    _simPdfRecipeList.Print();
    cout << "Labels: " << endl;
    _simPdfLabels.Print();
    return 0;
  }  
  
 return rec;
 
}


void MLPdfRecipe::createMasterSplitCat(const RooArgList &splitCatList) 
{
  _masterSplitCat = new RooMultiCategory("masterSplitCat", "Master Splitting Category",splitCatList) ;
  
//   RooArgList fundcats ;
//   TIterator* it = splitCatList.createIterator() ;
//   while( RooAbsArg* arg = static_cast<RooAbsArg*>(it->Next()))  
//     if( arg->InheritsFrom("RooAbsCategoryLValue") ) fundcats.add(*arg) ;
//   delete it ;
//   if( fundcats.getSize()!=0 ) 
//     _masterSplitCat = new RooSuperCategory("masterSplitCat", "Master Splitting Category",fundcats);
}
 

