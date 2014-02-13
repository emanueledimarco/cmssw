#include "MLFit/MLModel.hh"
#include "RooStringVar.h"

ClassImp(MLModel);

MLModel::MLModel(const char* name, const char* title) :
  TNamed(name, title),
  _masterSplitCat(0)
{
}

MLModel::~MLModel()
{
}

void MLModel::addPdf(const char *species, const char *obsLabel, const char* pdfname)
{
  TString key = makeKey(species, obsLabel);
  _pdfLookupTable.add(*(new RooStringVar(key, key, pdfname)));
  if (_labelList.FindObject(obsLabel) == NULL) _labelList.Add(obsLabel);
}

Int_t MLModel::numSpecies() const
{
  return _speciesList.GetSize(); 
}

void MLModel::createMasterSplitCat(const RooArgList &splitCatList)
{
  _masterSplitCat = new RooSuperCategory("masterSplitCat", "Master Splitting Category",splitCatList);
}

Bool_t MLModel::isSplit() const
{
  if (_masterSplitCat == 0) return kFALSE;
  if (_masterSplitCat->inputCatList().getSize() == 0) return kFALSE;
  return kTRUE;
}

void MLModel::addNoNormVars(RooArgSet vars)
{
  // Adds vars to the list of variables not to normalize over.
  // Nothing is done with these.  It's just for easy bookkeeping if
  // you want to use it.
  _noNormVars.add(vars);
}

TString MLModel::getPdfName(const char *species, const char *obs) const
{
  RooStringVar *pdfname = (RooStringVar*)_pdfLookupTable.find(makeKey(species,obs));
  if (pdfname == NULL) {
    cout << GetName() << ": Can't find species " << species << " and observable " << obs << " in lookup table:" << endl;
    _pdfLookupTable.Print("V");
    return "";
  }
  return TString(pdfname->getVal());
}

