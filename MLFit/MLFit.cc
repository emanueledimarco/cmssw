/**********************************************************************
 * MLFit
 * Nick Danielson / Amir Farbin
 *********************************************************************/
#include "TClass.h"
#include <fstream>
#include "TIterator.h"
#include "TFile.h"
#include "RooAddPdf.h"
#include "RooCustomizer.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooSimultaneous.h"
#include "RooStringVar.h"
#include "RooSuperCategory.h"
#include "RooUnblindUniform.h"
#include "RooUnblindCPAsymVar.h"
#include "RooArgusBG.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooRandom.h"
#include "MLFit/MLFit.hh"
#include "MLFit/MLModel.hh"
#include "MLFit/MLPdfRecipe.hh"
#include "MLFit/MLPdfFactory.hh"
#include "MLFit/MLSpecies.hh"
#include "MLFit/MLStrList.hh"

ClassImp(MLFit);

MLFit::MLFit(const char *name, const char* title) : 
  TNamed(name, title),
  _initialized(kFALSE)
{
  _recipeList = new TList();
  DataSets = new TList();
}

MLFit::~MLFit()
{
  // A lot of things should be deleted here... I've tried to delete all temporary items when they aren't needed any more.
  // but there are a lot of lists of items that need to exist as long as the MLFit object exists.  They should be deleted
  // here... As long as you're only creating a few MLFit objects this shouldn't be a problem.
  delete _recipeList;
  delete DataSets;
}

//===================================
// Flat Files/Observables...
//=================================
void MLFit::AddFlatFileColumn(RooAbsArg *obj,int column)
{
  if(column==0)
    _flatFileVariables.add(*obj);
  else
    if(_flatFileVariables.getSize()>column){
      RooAbsArg *blah=_flatFileVariables.at(column-1);
      if(blah){
	_flatFileVariables.replace(*blah,*obj);
	delete blah;
      }
    } else {
      while(_flatFileVariables.getSize()<column-1) {
	RooRealVar *blah=new RooRealVar("dummy","dummy",0,"dummy");
	_flatFileVariables.add( *blah);
      }  
      
      _flatFileVariables.add(*obj);
    }
  _observables.add(*obj);
}

RooAbsArg *MLFit::Observable(const char *obs) const
{
  RooAbsArg *ret=_observables.find(obs);
  if(ret)
    return ret;
  
  cout<< fName <<" Error: cannot find observable "<<obs<<"."<<endl;
  return 0;
}

RooRealVar *MLFit::RealObservable(const char *obs)
{
  RooAbsArg *ret=_observables.find(obs);
  if(ret) 
    if (ret->IsA()->InheritsFrom("RooRealVar"))
      return (RooRealVar*)ret;
    else {
      cout<< fName <<" Error: observable "<<obs<<" is different type."<<endl;
      return 0;
    }
    
  cout << fName <<" Error: cannot find observable "<<obs<<"."<<endl;
  return 0;
}

RooCategory *MLFit::CatObservable(const char *obs) const
{
  RooAbsArg *ret=_observables.find(obs);
  if(ret) 
    if (ret->IsA()->InheritsFrom("RooAbsCategory"))
      return (RooCategory*)ret;
    else {
      cout<<fName<<" Error: observable "<<obs<<" is different type."<<endl;
      return 0;
    }
  
  cout<<fName<<" Error: cannot find observable "<<obs<<"."<<endl;
  return 0;
}

RooArgSet MLFit::getParSet(const MLStrList &strList) const
{
  // Sometimes you want to easily get a RooArgList of parameters and
  // don't want to have to type myfit.getRealPar("mypar") for each one.
  // Use this function: getParSet(MLStrList("mean", "sigma", "kappa", "blah"))
  RooArgList arglist;
  TIterator *strIter = strList.MakeIterator();
  while (TObjString *str = (TObjString*)strIter->Next()) {
    RooAbsArg *arg = getRealPar(str->GetName());
    if (arg == 0) continue;
    arglist.add(*arg);
  }
  return arglist;
}

RooArgList MLFit::getObsList(const char* obs1, const char* obs2, const char* obs3, const char* obs4, const char* obs5, const char* obs6, const char* obs7, const char* obs8) const 
{
  MLStrList strlist;
  if (obs1 != 0) strlist.Add(obs1);
  if (obs2 != 0) strlist.Add(obs2);
  if (obs3 != 0) strlist.Add(obs3);
  if (obs4 != 0) strlist.Add(obs4);
  if (obs5 != 0) strlist.Add(obs5);
  if (obs6 != 0) strlist.Add(obs6);
  if (obs7 != 0) strlist.Add(obs7);
  if (obs8 != 0) strlist.Add(obs8);
  return getObsList(strlist);
}


RooArgList MLFit::getObsList(const MLStrList &strList) const
{

  RooArgList arglist;
  TIterator *strIter = strList.MakeIterator();
  while (TObjString *str = (TObjString*)strIter->Next()) {
    RooAbsArg *arg = Observable(str->GetName());
    if (arg == 0) continue;
    arglist.add(*arg);
  }
  delete strIter;
  return arglist;
}


//===============================
// Species...
//================================

void MLFit::addSpeciesFromModel(const char* model, const char* othermodel) 
{
  // This takes the species used by othermodel and adds them to model.
  MLModel *theOtherModel = getModel(othermodel);
  if (theOtherModel == 0) return;
  MLModel *theModel = getModel(model);
  if (theModel == 0) return;

  TIterator *specIter = theOtherModel->_speciesList.MakeIterator();
  while (MLSpecies* spec = (MLSpecies*)specIter->Next()){
    theModel->_speciesList.Add(spec);
  }
  delete specIter;
  
}
void MLFit::addSpecies(const char* model, const char* sname, const char* stitle, RooAbsReal *svar)
{
  // Adds a species to the MLFit.  But sometimes you might want to use some arbitrary variable or function
  // to measure the yields instead of a simple f_sig or N_sig.  So create the yield variable you want to associate with
  // this species and pass that as svar.

  // Make sure the model exists.  If it's "default" and doesn't exist, create it, unless you already have a model created somewhere.
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if ( (theModel == 0 && TString(model) != "default") ||
       (theModel == 0 && TString(model) == "default" && _modelList.GetSize() != 0)) {
    cout << "Model: " << model << " doesn't exist!" << endl;
    return;
  }
  
  if (theModel == 0 && TString(model) == "default") {
    addModel("default", "Default Model");
    theModel = (MLModel*)_modelList.FindObject(model);
  }
  
  // Add this species to the species list
  MLSpecies *s = (MLSpecies*)theModel->_speciesList.FindObject(sname);
  if (s == 0) {
    s = new MLSpecies(sname,stitle);
    _speciesList.Add(s);
  }
  theModel->_speciesList.Add(s);
  //RooAbsReal *yieldvar = s->getYieldVar();
  //_parameterSet.add(*yieldvar);
}

/*
void MLFit::addSpecies(const char* model, const char* sname, const char* stitle)
{
  // Make sure the model exists.  If it's "default" and doesn't exist, create it.
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0 && TString(model) != "default") {
    cout << "Model: " << model << " doesn't exist!" << endl;
    return;
  }
  
  if (theModel == 0 && TString(model) == "default") {
    addModel("default", "Default Model");
    theModel = (MLModel*)_modelList.FindObject(model);
  }
  
  // Add this species to the species list
  MLSpecies *s = (MLSpecies*)theModel->_speciesList.FindObject(sname);
  if (s == 0) {
    s = new MLSpecies(sname,stitle);
    _speciesList.Add(s);
  }
  theModel->_speciesList.Add(s);

}
*/
void MLFit::addSpecies(const char* sname, const char* stitle) 
{
  addSpecies("default", sname, stitle);
}

MLSpecies* MLFit::getSpecies(const char *sname) const
{
  MLSpecies *species = (MLSpecies*)_speciesList.FindObject(sname);
  if (species == 0) {
    cout << fName << " Error: cannot find species " << sname << "." << endl;
    cout << "Available species: " << endl;
    _speciesList.Print();
    return 0;
  }

  return species;
  
}
RooAbsReal* MLFit::getCoefficient(const char *model, TString species)
{
  MLModel *theModel = getModel(model);
  if (model == NULL) return 0;
  return getCoefficient(*theModel, species);
}

RooAbsReal* MLFit::getCoefficient(MLModel &theModel, TString species) 
{
  TString name;
  if (!theModel.isSplit()) {
    name = "N_" + species;
  } else {
    name = "coef_" + species + "_" + TString(theModel._masterSplitCat->getLabel());  
  }

  RooAbsReal *coef = (RooAbsReal*)_fracList.find(name);
  if (coef == 0) {
    cout << GetName() << ": Can't find coefficient " << name << " in coefficient list:" << endl;
    _fracList.Print("V");
    return 0;
  }
  return coef;
}

Int_t MLFit::numSpecies() 
{
  return _speciesList.GetSize();
}


//===================
// Models...
//====================
void MLFit::addModel(const char* name, const char* title)
{
  _modelList.Add(new MLModel(name, title));
}

void MLFit::addNoNormVars(const char* model, MLStrList &vars)
{
  MLModel *theModel = getModel(model);
  if (theModel == 0) return;
  RooArgList arglist(getObsList(vars));
  theModel->addNoNormVars(arglist);
}

RooArgSet MLFit::getNoNormVars(const char* model) const
{
  MLModel *theModel = getModel(model);
  if (theModel == 0) return RooArgSet();
  return theModel->getNoNormVars();
}

RooArgSet MLFit::getNormVars(TString model) const
{
  MLModel *theModel = getModel(model);
  if (theModel == 0) return RooArgSet();
  RooArgSet noNormVars(theModel->getNoNormVars());
  RooArgSet normVars(_flatFileVariables);
  normVars.remove(noNormVars);
  return normVars;
}


MLModel* MLFit::getModel(const char *model) const
{
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0) cout << fName << " Error: model " << model << " doesn't exist!" << endl;
  return theModel;
}


//=============================
// Recipes...
//==================================

void MLFit::addPdfWName(const char* species, const char *obs, TString pdftype, TString pdfbase) 
{
  addPdfWName("default", species, MLStrList(obs), pdftype, new TList(), pdfbase);
}
void MLFit::addPdfWName(const char* model, const char* species, const char *obs, TString pdftype, TString pdfbase) 
{
  addPdfWName(model, species, MLStrList(obs), pdftype, new TList(), pdfbase);
}

void MLFit::addPdf(const char* species, const char *obs, TString pdftype) 
{
  addPdf("default", species, obs, pdftype);
}

void MLFit::addPdf(const char *model, const char* species, const char *obs, TString pdftype) 
{
  if (Observable(obs) == 0) return;
  addPdf(model, species, MLStrList(obs), pdftype, new TList());
}

void MLFit::addPdf(const char* species, const MLStrList &obs, TString pdftype) 
{
  addPdf("default", species, obs, pdftype);
}

void MLFit::addPdf(const char* model, const char* species, const MLStrList &obs, TString pdftype) 
{
  addPdf(model, species, obs, pdftype, new TList());
}

void MLFit::addPdf(const char* species, const MLStrList &obs, TString pdftype, const TList& args) 
{
  addPdf("default", species, obs, pdftype, args);
}

void MLFit::addPdf(const char* model, const char* species, const MLStrList &obs, TString pdftype, const TList& args) 
{
  addPdfWName(model, species, obs, pdftype, args, "");
}

void MLFit::add2DPdf(const char* model, const char* species, const MLStrList &obs, TString pdftype, const TList& args) 
{
  add2DPdfWName(model, species, obs, pdftype, args, "");
}

void MLFit::addPdfWName(const char* species, const MLStrList &obs, TString pdftype, const TList &args, TString pdfbase) 
{
  addPdfWName("default", species, obs, pdftype, args, pdfbase);
}

void MLFit::addPdfWName(const char* model, const char* species, const MLStrList &obs, TString pdftype, const TList &args, TString pdfbase) 
{
  // Make sure the species and all the observables exist...
  if (getSpecies(species) == 0) return;
  RooArgList obsList;
  TIterator *obsIter = obs.MakeIterator();
  while (TObjString *arg = (TObjString*)obsIter->Next()){
    RooAbsArg *obsArg = Observable(arg->GetName());
    if (obsArg == 0) return;
    obsList.add(*obsArg);
  }
  delete obsIter;

  // Make sure the model exists.  If it's "default" and doesn't exist, create it.
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0 && TString(model) != "default") return;
  if (theModel == 0 && TString(model) == "default") {
    addModel("default", "Default Model");
    theModel = (MLModel*)_modelList.FindObject(model);
  }

  // Give the PDF a good name.
  RooAbsArg *firstObs = obsList.first();
  if (pdfbase == "") pdfbase = TString(species)+"_"+firstObs->GetName();
  TString pdfname = pdfbase + "_" + pdftype;

  // Get the parameters needed by this PDF.
  RooArgList parameters(MLPdfFactory::makeParameters(pdftype,pdfbase,args));
  
  // Add entries to the various tables...
  theModel->addPdf(species, firstObs->GetName(), pdfname);
  _recipeList->Add(new MLPdfRecipe(pdfname, pdfname, pdftype, obsList, RooArgSet(), parameters, args));

}
void MLFit::add2DPdfWName(const char* model, const char* species, const MLStrList &obs, TString pdftype, const TList &args, TString pdfbase) 
{
  // Make sure the species and all the observables exist...
  if (getSpecies(species) == 0) return;
  RooArgList obsList;
  TIterator *obsIter = obs.MakeIterator();
  while (TObjString *arg = (TObjString*)obsIter->Next()){
    RooAbsArg *obsArg = Observable(arg->GetName());
    if (obsArg == 0) return;
    obsList.add(*obsArg);
  }
  delete obsIter;

  // Make sure the model exists.  If it's "default" and doesn't exist, create it.
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0 && TString(model) != "default") return;
  if (theModel == 0 && TString(model) == "default") {
    addModel("default", "Default Model");
    theModel = (MLModel*)_modelList.FindObject(model);
  }

  // Give the PDF a good name.
  //RooAbsArg *firstObs = obsList.first();
  RooAbsArg *firstObs = obsList.at(0);
  RooAbsArg *secondObs = obsList.at(1);
  if (pdfbase == "") pdfbase = TString(species)+"_"+firstObs->GetName();
  TString pdfname = pdfbase + "_" + pdftype;

  // Get the parameters needed by this PDF.
  RooArgList parameters(MLPdfFactory::makeParameters(pdftype,pdfbase,args));
  
  // Add entries to the various tables...
  theModel->addPdf(species, firstObs->GetName(), pdfname);
  theModel->addPdf(species, secondObs->GetName(), pdfname);
  _recipeList->Add(new MLPdfRecipe(pdfname, pdfname, pdftype, obsList, RooArgSet(), parameters, args));

}

void MLFit::addSimPdf(const char* model, const char* species, const char *obs, const char* catvar, TString pdfbase)
{
  // Make sure the species exists...
  if (getSpecies(species) == 0) return;

  // Make sure category exists
  RooCategory *cat = CatObservable(catvar);
  if (cat == 0) return;
 
  // Make sure the model exists.  If it's "default" and doesn't exist, create it.
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0 && TString(model) != "default") return;
  if (theModel == 0 && TString(model) == "default") {
    addModel("default", "Default Model");
    theModel = (MLModel*)_modelList.FindObject(model);
  }
  
  // Give the PDF a good name.
  if (pdfbase == "") pdfbase = TString(species)+"_"+obs;
  TString pdfname = pdfbase;

  // Add entries to the various tables...
  theModel->addPdf(species, obs, pdfname);
  _recipeList->Add(new MLPdfRecipe(pdfname, pdfname, obs, *cat));
}

void MLFit::addPdfToSimPdf(const char *model, const char* species, const char* simpdfname, const MLStrList &obs, TString catLabel, TString pdftype, const TList &args, TString pdfbase)
{
  // Make sure the species and all the observables exist...
  if (getSpecies(species) == 0) return;
  RooArgList obsList;
  TIterator *obsIter = obs.MakeIterator();
  while (TObjString *arg = (TObjString*)obsIter->Next()){
    RooAbsArg *obsArg = Observable(arg->GetName());
    if (obsArg == 0) return;
    obsList.add(*obsArg);
  }
  delete obsIter;
  
  // Make sure the model exists.
  MLModel *theModel = getModel(model);
  if (theModel == 0) return;
  
  // Give the PDF a good name.
  RooAbsArg *firstObs = obsList.first();
  if (pdfbase == "") pdfbase = TString(species)+"_"+firstObs->GetName();
  TString pdfname = pdfbase + "_" + pdftype;
  
  // Get the parameters needed by this PDF.
  RooArgList parameters(MLPdfFactory::makeParameters(pdftype,pdfbase,args));
  
  // Add entries to the various tables...
  // This PDF recipe should already have been created (as a RooSimultaneous), so let's get it.
  MLPdfRecipe *theRecipe = getRecipe(model, species, simpdfname);
  if (theRecipe == 0) {
    cout << "MLFit: Trying to add pdf " << pdfname << " to simultaneous pdf for " << model << ", " << species << ", " << firstObs->GetName() << endl;
    cout << "       but you haven't created the sim pdf yet!" << endl;
    return;
  }
  theRecipe->_simPdfRecipeList.Add(new MLPdfRecipe(pdfname, pdfname, pdftype, obsList, RooArgSet(), parameters, args));
  theRecipe->_simPdfLabels.add(*(new RooStringVar(catLabel, catLabel, catLabel)));
  theRecipe->_protoParams.add(parameters);
  
}

void MLFit::addRF(const MLStrList &obs, TString rftype, const TList &args, TString rfbase)
{
  // Add a resolution function for later use.  This command doesn't actually add the resolution function to your model.  It just
  // makes it available so you can add convoluted pdfs with it.

  // Check that the observables exist.
  RooArgList obsList;
  TIterator *obsIter = obs.MakeIterator();
  while (TObjString *arg = (TObjString*)obsIter->Next()){
    RooAbsArg *obsArg = Observable(arg->GetName());
    if (obsArg == 0) return;
    obsList.add(*obsArg);
  }
  delete obsIter;
    
  // Give the RF a good name.
  RooAbsArg *firstObs = obsList.first();
  if (rfbase == "") rfbase = firstObs->GetName();
  TString pdfname = rfbase + "_RF";

   // Get the parameters needed by this PDF.
  RooArgList parameters(MLPdfFactory::makeParameters(rftype,rfbase,args));
  
  // Add entries to the various tables...
  _rfRecipeList.Add(new MLPdfRecipe(pdfname, pdfname, rftype, obsList, RooArgSet(), parameters, args));


}

void MLFit::addConvolutedPdfToSimPdf(const char *model, const char* species, const char* simpdfname, const MLStrList &obs, TString catLabel, TString pdftype, const char *rf, const TList &args, TString pdfbase)
{
  // Make sure the species and all the observables exist...
  if (getSpecies(species) == 0) return;
  RooArgList obsList;
  TIterator *obsIter = obs.MakeIterator();
  while (TObjString *arg = (TObjString*)obsIter->Next()){
    RooAbsArg *obsArg = Observable(arg->GetName());
    if (obsArg == 0) return;
    obsList.add(*obsArg);
  }
  delete obsIter;
  
  // Make sure the model exists.
  MLModel *theModel = getModel(model);
  if (theModel == 0) return;
  
  // Make sure the resolution function exists.
  MLPdfRecipe *theRF = getRF(rf);
  if (theRF == 0) return;
  
  // Give the PDF a good name.
  RooAbsArg *firstObs = obsList.first();
  if (pdfbase == "") pdfbase = TString(species)+"_"+firstObs->GetName();
  TString pdfname = pdfbase + "_" + pdftype;
  
  // Get the parameters needed by this PDF.Making workdir/MLFit//RhhDoubleGaussian.d
  RooArgList parameters(MLPdfFactory::makeParameters(pdftype,pdfbase, new TList()));
  
  // Add entries to the various tables...
  // This PDF recipe should already have been created (as a RooSimultaneous), so let's get it.
  MLPdfRecipe *theRecipe = getRecipe(model, species, simpdfname);
  if (theRecipe == 0) {
    cout << "MLFit: Trying to add pdf " << pdfname << " to simultaneous pdf for " << model << ", " << species << ", " << firstObs->GetName() << endl;
    cout << "       but you haven't created the sim pdf yet!" << endl;
    return;
  }
  MLPdfRecipe *rec = new MLPdfRecipe(pdfname, pdfname, pdftype, obsList, RooArgSet(), parameters, args);
  rec->setRF(*theRF);
  //  _recipeList->Add(rec);
  theRecipe->_simPdfRecipeList.Add(rec);
  theRecipe->_simPdfLabels.add(*(new RooStringVar(catLabel, catLabel, catLabel)));
  theRecipe->_protoParams.add(parameters);
  
}

void MLFit::addConvolutedPdf(const char *model, const char *species, const MLStrList &obs, TString pdftype, const char *rf, const TList &args, TString pdfbase)
{
  // Make sure the species and all the observables exist...
  if (getSpecies(species) == 0) return;
  RooArgList obsList;
  TIterator *obsIter = obs.MakeIterator();
  while (TObjString *arg = (TObjString*)obsIter->Next()){
    RooAbsArg *obsArg = Observable(arg->GetName());
    if (obsArg == 0) return;
    obsList.add(*obsArg);
  }
  delete obsIter;
  
  // Make sure the model exists.  If it's "default" and doesn't exist, create it.
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0 && TString(model) != "default") return;
  if (theModel == 0 && TString(model) == "default") {
    addModel("default", "Default Model");
    theModel = (MLModel*)_modelList.FindObject(model);
  }

  // Make sure the resolution function exists.
  MLPdfRecipe *theRF = getRF(rf);
  if (theRF == 0) return;
  
  // Give the PDF a good name.
  RooAbsArg *firstObs = obsList.first();
  if (pdfbase == "") pdfbase = TString(species)+"_"+firstObs->GetName();
  TString pdfname = pdfbase + "_" + pdftype;
 
  // Get the parameters needed by this PDF.
  RooArgList parameters(MLPdfFactory::makeParameters(pdftype,pdfbase,args));

  // Add entries to the various tables...
  theModel->addPdf(species, firstObs->GetName(), pdfname);
  MLPdfRecipe *rec = new MLPdfRecipe(pdfname, pdfname, pdftype, obsList, RooArgSet(), parameters, args);
  rec->setRF(*theRF);
  _recipeList->Add(rec);

}


void MLFit::addPdfFromModel(const char *model, const char* species, const char* obs, const char *othermodel)
{
  // You have some model, say, myMesModel.  And you are making another model, say, myCPModel.  You want these two models to
  // share their signal mes PDFs.  Then use this function.
  MLModel *theModel = getModel(model);
  if (theModel == 0) return;

  MLModel *theOtherModel = getModel(othermodel);
  if (theOtherModel == 0) return;

  TString pdfname = theOtherModel->getPdfName(species, obs);
  theModel->addPdf(species, obs, pdfname);
}

void MLFit::addPdfsFromModel(const char* model, const char* obs, const char* othermodel) 
{
  // Similar to addPdfFromModel, but this adds all the pdfs for a given observable at once (signal, bkg, etc..).
  MLModel *theModel = getModel(model);
  if (theModel == 0) return;

  TIterator *specIter = theModel->_speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()){
    addPdfFromModel(model, spec->GetName(), obs, othermodel);
  }
  delete specIter;
}

void MLFit::blindUniform(TString argname, TString blindstring, Double_t scale)
{
  // Blind a parameter using RooUnblindUniform

  // Get arg to blind
  RooRealVar *arg = (RooRealVar*)getProtoParam(argname,kFALSE);
  if (arg == 0) {
    // There is no normal parameter with this name.  But it might be a name of a yield variable.
    // Yield variables don't get created until you buildModel().  So we'll add this to a list, then
    // do the blinding after building the yield variables.  Annoying? Yes.
    _blindedYieldNames.add(*(new RooRealVar(argname,blindstring,scale)));
    return;
  }
  blindUniform(*arg, blindstring, scale);
}
void MLFit::blindUniform(RooRealVar &arg, TString blindstring, Double_t scale)
{
  // Replace it with the blinded paramter, and add it to a list of
  // blinded parameters just in case you want to access it later
  // for some reason
  _blindParSet.add(arg);
  RooCategory *bs = new RooCategory(TString(arg.GetName())+"_bs",TString(arg.GetName())+"_bs");
  bs->defineType("Blind");
  bs->defineType("Unblind");
  RooUnblindCPAsymVar *blindvar = new RooUnblindCPAsymVar(TString(arg.GetName())+"_blinded", TString(arg.GetName())+"_blinded", blindstring, arg);
  replace(arg.GetName(), *blindvar);
}


void MLFit::split(const char* cat, const char* splitarg)
{
  RooAbsArg *arg = getProtoParam(splitarg);
  if (arg == 0) return;
  split(cat, *arg);
}


void MLFit::split(const char* cat, RooArgSet splitSet)
{
  RooAbsCategory *category = dynamic_cast<RooAbsCategory*>(CatObservable(cat));
  if (category == 0) return;
  // First check that all the splitPars are really parameters.  Then add the par and the category to the
  // appropriate lists...
  TIterator *parIter = splitSet.createIterator();
  while (RooAbsArg* par = (RooAbsArg*)parIter->Next()){
    if (getProtoParam(par->GetName()) == 0) {
      cout << "MLFit: can't split parameter" << par->GetName() 
	   << " because it's not a parameter of " << fName;
      continue;
    }
    _splitArgList.add(*par);
    _splitCatList.add(*category);
  }
  delete parIter;
}

void MLFit::noSplitFrac(TString catname)
{
  // If you want to split by something like run, you probably don't want
  // e_sig_Run1, e_bkg_Run2, etc in front of your yields.
  // Call this function after splitting.
  
  RooCategory *cat = dynamic_cast<RooCategory*>(CatObservable(catname));
  if (cat == NULL) return;
  
  _noSplitFracSet.add(*cat);
}



void MLFit::addPdfCopy(const char* species, const char* obs, const char* otherspec)
{
  addPdfCopy("default", species, obs, otherspec);
}

void MLFit::addPdfCopy(const char* model, const char* species, const char* obs, const char* otherspec)
{
  MLModel *theModel = getModel(model);
  if (theModel == 0) return;
  if (getSpecies(species) == 0 || getSpecies(otherspec) == 0) return;
  TString pdfname = theModel->getPdfName(otherspec, obs);
  theModel->addPdf(species,obs,pdfname);  
}

TList* MLFit::getRecipeList() 
{
  return _recipeList;
}

MLPdfRecipe* MLFit::getRecipe(const char* species, const char* obs)
{
  return getRecipe("default", species, obs);
}


MLPdfRecipe* MLFit::getRecipe(const char* model, const char* species, const char* obs)
{
  if (getSpecies(species) == 0) return 0;
  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0) return 0;
  
  TString pdfname = theModel->getPdfName(species, obs);
  MLPdfRecipe *rec = ((MLPdfRecipe*)_recipeList->FindObject(pdfname));
  if (rec == 0) {
    cout << GetName() << ": Can't get recipe for species " << species << ", observable " << obs << ", and model " << model << endl;
    cout << "I'm looking for pdfname: " << pdfname << ", but can't find it in" << endl;
    _recipeList->Print();
    return 0;
  }
  return rec;
}

//========================================
// PDFs...
//========================================

RooAbsPdf* MLFit::getPdf(const char* pname) const 
{
  RooAbsPdf *pdf = (RooAbsPdf*)_pdfList.find(pname);
  if (pdf == 0) cout << "MLFit " << fName << ": Can't find pdf " << pname << endl;
  return pdf;
}

RooAbsPdf* MLFit::getPdfComponent(const char* model, const char* species) const
{
  MLModel *theModel = getModel(model);
  if (theModel == 0) return 0;
  return getPdfComponent(*theModel, species);
}


RooAbsPdf* MLFit::getPdfComponent(const MLModel &theModel, const char* species) const
{
  TString pdfname;
  if (theModel.numSpecies() == 1 && !theModel.isSplit()) {
    pdfname = theModel.GetName();
  } else if (theModel.isSplit()) {
    pdfname = TString(theModel.GetName())+"_"+species+"_"+theModel._masterSplitCat->getLabel();
  } else {
    pdfname = TString(theModel.GetName())+"_"+species;
  }
  RooAbsPdf *pdf = (RooAbsPdf*)_pdfList.find(pdfname);

  if (pdf == 0) cout << "MLFit::getPdfComponent can't find the pdf " << pdfname << endl;
  return pdf;
}

RooAbsPdf* MLFit::getPdfComponent(const char* model, const char* species, const char* obs) const
{
  MLModel *theModel = getModel(model);
  if (theModel == 0) return 0;
  return getPdfComponent(*theModel, species, obs);
}

RooAbsPdf* MLFit::getPdfComponent(const MLModel &theModel, const char* species, const char* obs) const
{
  if (getSpecies(species) == 0 || Observable(obs) == 0) return 0;
  TString pdfname = theModel.getPdfName(species, obs);

  MLPdfRecipe *rec = (MLPdfRecipe*)_recipeList->FindObject(pdfname);
  TString name = rec->getSplitName();
  RooAbsPdf *pdf = (RooAbsPdf*)_pdfList.find(name);
  if (pdf == 0) cout << GetName() << ": Can't find pdf " << name << endl;
  return pdf;
}



TList* MLFit::getRecipeListBySpecies(MLModel &theModel, TString species)
{
  if (getSpecies(species) == 0) return 0;
  TList *recipeList = new TList();
  MLStrList obsList(theModel.getLabelList());
  TIterator *obsIter = obsList.MakeIterator();
  while (TObjString *obs = (TObjString*)obsIter->Next()){
    TString pdfname = theModel.getPdfName(species, obs->GetName());
    recipeList->Add((MLPdfRecipe*)getRecipe(theModel.GetName(), species, obs->GetName()));
  }
  return recipeList;
}


//=========================================
// Resolution functions...
//=========================================
MLPdfRecipe* MLFit::getRF(const char* rf) const
{
  MLPdfRecipe *theRF = (MLPdfRecipe*)_rfRecipeList.FindObject(rf);
  if (theRF == 0) {
    cout << GetName() << ": Can't find resolution function " << rf << endl
	 << "Please, choose from: " ;
    TIter iter(&_rfRecipeList) ;
    TObject* obj ;
    while( (obj=iter.Next()) ) cout << obj->GetName() << ", " ;
    cout << endl ;
  }
  return theRF;
}


//=========================================
// Parameters...
//=========================================

RooArgSet MLFit::getProtoParams()
{
  // Before we build any models, we have a list of parameters.  These parameters might be split
  // later when we build the models, but for now they're just place-holders.  They're my
  // prototype parameters.  This function loops over all the PDF recipes and gets each of their
  // prototype parameters.
  RooArgSet params;
  TIterator *recIter = _recipeList->MakeIterator();
  while (MLPdfRecipe *rec = (MLPdfRecipe*)recIter->Next()){
    params.add(rec->getProtoParams(),kTRUE);
    //if (rec->isConvoluted()) params.add(rec->_rf->getProtoParams());
  }
  delete recIter;

  // The species yield coefficients also have parameters which might be split later.  Add
  // these prototype parameters to the list
  TIterator *specIter = _speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()){
    params.add(spec->getProtoParams(),kTRUE);
  }
  delete specIter;
  
  return params;
}

RooRealVar* MLFit::getProtoParam(const char* parname, Bool_t verbose)
{
  RooRealVar *var = (RooRealVar*)getProtoParams().find(parname);
  if (var == 0 && verbose) {
    cout << fName << " Error: Can't find prototype parameter: " << parname << " in list" << endl;
    getProtoParams().Print("V");
  }
  return var;
  
}

RooRealVar *MLFit::getSplitParameter(RooAbsReal *arg) const {
  
  TString name = getSplitParameterName(arg->GetName());
  RooRealVar *var = (RooRealVar*)_parameterSet.find(getSplitParameterName(arg->GetName()));
  if (var == NULL) {
    cout << GetName() << ": Can't find parameter: " << name << " in list " << endl;
    _parameterSet.Print("V");
  }
  return var;
}

TString MLFit::getSplitParameterName(TString name) const
{
  // Takes the name of an unsplit parameter and returns its split name (assuming the current values of
  // the splitting categories).  Example: getSplitParameterName("mean") = "mean_Cat3"
  RooAbsArg *arg = _splitArgList.find(name);
  if (arg == 0) return name;

  // HUGE HACK for simultaneous BReco fit
  RooAbsCategory *splitCat = (RooAbsCategory*)_splitCatList.at(_splitArgList.index(arg));
  return name + "_" + splitCat->getLabel();
}

RooAbsArg* MLFit::getPar(TString parname) const
{
  RooAbsArg *p = _parameterSet.find(parname);
  if (p == 0) {
    cout << GetName() << " Error: " << parname << " isn't a parameter of this fit!" << endl;
    return 0;
  }
  return p;
}


RooRealVar* MLFit::getRealPar(const char* parname) const
{
  RooArgSet params(_parameterSet);
  params.add(_blindParSet);
  RooAbsArg *p = params.find(parname);
  if (p == 0) {
    cout << GetName() << " Error: " << parname << " isn't a parameter of this fit!" << endl;
    return 0;
  }
  if (!p->IsA()->InheritsFrom("RooRealVar")){
    cout << GetName() << " Error: " << parname << " is a parameter, but isn't a RooRealVar!" << endl;
  }

  return (RooRealVar*)p;
}


void MLFit::bind(const MLStrList &parList, const char* name, const char* title)
{

  RooAbsArg *firstPar = getProtoParam(((TObjString*)parList.First())->GetName());
  if (firstPar == 0) return;
  
  firstPar->SetName(name);
  firstPar->SetTitle(title);
  
  TIterator *parIter = parList.MakeIterator();
  while (TObjString* par = (TObjString*)parIter->Next()) {
    //if (getProtoParam(par->GetName()) == 0) continue;
    replace(par->GetName(), *firstPar);
  }
  delete parIter;

}

void MLFit::bindFractions(const MLStrList &parList, const char* name, const char* title)
{
  TString firstname  = ((TObjString*)parList.First())->GetName();
  RooAbsReal *firstPar(0);
  TIterator *specIter = _speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()){
    if (TString(spec->getFracVar()->GetName()) == firstname) {
      firstPar = spec->getFracVar();
      break;
    }
  }
  delete specIter;
  if (firstPar == 0) {
    cout << GetName() << " Error: can't find any fraction named " << firstname << endl;
    return;
  }
  firstPar->SetName(name);
  firstPar->SetTitle(title);

  TIterator *parIter = parList.MakeIterator();
  while (TObjString* par = (TObjString*)parIter->Next()) {
    replaceFraction(par->GetName(), *firstPar);
  }
  delete parIter;

}

void MLFit::rename(const char *oldname, const char* newname, const char* newtitle) 
{
  RooAbsReal *p = getProtoParam(oldname);
  p->SetName(newname);
  p->SetTitle(newtitle);
}


void MLFit::replace(TString oldpar, RooAbsArg &newpar)
{
  // Replace a parameter with a new one.
  
  // We have to loop over all PDFs to find the ones that depend on this parameter, and
  // replace all instances of this parameter.
  TIterator *recIter = _recipeList->MakeIterator();
  while (MLPdfRecipe* rec = (MLPdfRecipe*)recIter->Next()) {
    rec->replace(oldpar, newpar,kFALSE);
  }
  delete recIter;

  // oldpar might be a parameter *not* of a pdf, but of a *yield variable*, like an asymmetry A_kpi.
  TIterator *specIter = _speciesList.MakeIterator();
  while (MLSpecies *spec = dynamic_cast<MLSpecies*>(specIter->Next())) {
    spec->replace(oldpar, newpar, kFALSE);
  }

  // oldpar might be a yield variable (like N_sig), so also loop over all the yield vars...
  // (This is not really useful to the user, since the yield variables aren't created
  // until you buildModel(), at which point, you can't replace parameters any more!
  // This is used internally by MLFit for blinding yield parameters.  It calls
  // this function right after it creates the yield variables, but before it creates the PDFs.)
  TIterator *fracIter = _fracList.createIterator();
  while (RooAbsReal *yieldvar = dynamic_cast<RooAbsReal*>(fracIter->Next())) {
    if (TString(yieldvar->GetName()) == oldpar) {
      _fracList.replace(*yieldvar,newpar);
      _parameterSet.replace(*yieldvar,newpar);
      break;
    }
  }
  delete fracIter;
}

void MLFit::replaceFraction(const char *oldf, RooAbsReal &newf)
{
  TIterator *specIter = _speciesList.MakeIterator();
  while (MLSpecies* spec = (MLSpecies*)specIter->Next()){
    if (spec->getFracVar()->GetName() == TString(oldf)) {
      spec->_fracVar = &newf;
    }    
  }
  delete specIter;
}


void MLFit::addYieldVars(RooArgSet &yieldVars) 
{
  _protoFracSet.add(yieldVars);
}


//=====================
// Datasets...
//=======================
void MLFit::addDataSetFromAsciiFile(TString dataname, const char *file, const char* commonPath, const char* opt)
{
  RooDataSet* data = RooDataSet::read(file, _flatFileVariables, opt, commonPath);
  if (data == 0) {
    cout << GetName() << ": Can't get data set in file " << file << endl;
    return;
  }
  addDataSet(dataname, *data);
}

void MLFit::addDataSetFromRootFile(TString dataname, const char *nameinfile, const char* file)
{
  TFile *f = new TFile(file);
  if (f == 0) {
    cout << GetName() << ": Can't open file " << file << endl;
    return;
  }
  RooDataSet *data = (RooDataSet*)f->Get(nameinfile);
  if (data == 0) {
    cout << GetName() << ": Can't get data set " << nameinfile << " from file " << file << endl;
    cout << "File contains: " << endl;
    f->ls();
    return;
  }
  addDataSet(dataname, *data);
}

void MLFit::addDataSetFromRootFile(TString dataname, const char *nameinfile, const char* file, const RooArgSet& vars)
{
  TFile *f = new TFile(file); 

  f->ls(); 

  if (f == 0) {
    cout << GetName() << ": Can't open file " << file << endl;
    return;
  }

  
  TTree *tree = (TTree*)f->Get(nameinfile);

  RooDataSet *data = new RooDataSet("dataset2", "dataset2",tree,vars);

  //data->Print();
  if (data == 0) {
    cout << GetName() << ": Can't get data set " << nameinfile << " from file " << file << endl;
    cout << "File contains: " << endl;
    f->ls();
    return;
  }
  addDataSet(dataname, *data);
}

void MLFit::addDataSet(TString dataname, RooDataSet &data)
{
  // Add a dataset.
  data.SetName(dataname);
  DataSets->Add(&data);

  // If you wanted to add some columns to the dataset that are functions of other columns in the dataset, first you must have
  // used addDataSetFormula.  Now we actually calculate those columns.
  TIterator *formulaIter = _dataSetFormulas.createIterator();
  while (RooAbsReal *formula = (RooAbsReal*)formulaIter->Next()) {addDataSetColumn(dataname,*formula, kTRUE);}
  delete formulaIter;
}

RooDataSet* MLFit::getDataSet(const char *dataname)
{
  RooDataSet* data = (RooDataSet*)DataSets->FindObject(dataname);
  if (data == 0) cout << fName << ": Can't find data set: " << dataname << endl;
  return data;
}

void MLFit::addDataSetColumn(const char *dataname, RooAbsArg &argtoadd, Bool_t silent) 
{
  RooDataSet *data = getDataSet(dataname);
  if (data == 0) return;
  RooAbsArg *par = data->addColumn(argtoadd);
  _observables.add(*par,silent);
}

/*
void MLFit::addDataSetColumn(TString dataname, TString argtoadd)
{
  // This adds a column to your dataset, calculated from some formula variable (argtoadd)
  // which must have been added to the observables list with addObservable() earlier.
  RooDataSet *data = getDataSet(dataname);
  if (data == 0) return;
  RooAbsArg *arg = _dataSetFormulas.find(argtoadd);
  if (arg == 0) {
    cout << GetName() << " Error: Can't find formula " << argtoadd << " in " << endl;
    _dataSetFormulas.Print();
    return;
  }
  
  RooAbsArg *par = data->addColumn(*arg);
  _observables.replace(*_observables.find(par->GetName()),*par);
}
*/
//====================================
// Models...
//================================
RooAbsPdf* MLFit::buildModel(const char* model, Bool_t verbose) 
{
  if (verbose) cout << "Building model: " << model << endl;

  MLModel *theModel = (MLModel*)_modelList.FindObject(model);
  if (theModel == 0) return 0;


  // Split all the parameters that need splitting... Only do this the first time we build a model.
  if (verbose) cout << "Checking if parameters have been split." << endl;
  if (!_initialized) {
    if (verbose) cout << "Parameters have not been split yet.  Splitting..." << endl;
    // First make a copy of prototype parameters; loop through them, splitting anything that needs to be split.
    RooArgSet protoParams(getProtoParams());
    TIterator *iter = protoParams.createIterator();
    while (RooAbsArg* arg = (RooAbsArg*)iter->Next()) {
      if (_splitArgList.find(arg->GetName()) == 0) {
	// Unsplit parameter
	_parameterSet.add(*arg);
      } else {
	// Split parameter
	if (verbose) cout << "Splitting parameter: " << arg->GetName() << endl;
	RooAbsCategory *splitCat = (RooAbsCategory*)_splitCatList.at(_splitArgList.index(arg));
	TIterator *stateIter = splitCat->typeIterator();
	while (RooCatType *state = (RooCatType*)stateIter->Next()){
	  TString name = TString(arg->GetName()) + "_" + state->GetName();
	  TString title = TString(arg->GetTitle()) + " " + state->GetName();
	  RooAbsArg *clone = (RooAbsArg*)arg->Clone(name);
	  clone->SetTitle(title);
	  _parameterSet.add(*clone);
	  if (verbose) cout << "Created split parameter: " << clone->GetName() << endl;
	}
	delete stateIter;	
      }//endif
    }
    delete iter; 

    // Let's keep track of which categories split each PDF.  This is useful for naming PDFs and for creating only the PDFs we need to
    // (no need to create 5 copies of the fisher pdf if it's not split by tagging category).  This means setting its _masterSplitCat
    // correctly..
    if (verbose) cout << "Begin splitting the PDF recipes" << endl;
    TIterator *recIter = _recipeList->MakeIterator();
    while (MLPdfRecipe *rec = (MLPdfRecipe*)recIter->Next()) splitPdfRecipe(*rec, verbose);
    
    delete recIter;

    _initialized = kTRUE;
  }//endsplitting

  // Create this model's master splitting category...
  // Loop over all the split parameters and see if this model depends on that them...
  TIterator *argIter = _splitArgList.createIterator();
  RooArgSet splitCatSet;
  while( RooAbsArg* arg = static_cast<RooAbsArg*>(argIter->Next()) )  {
    // And loop over the PDFs in this model to see if any depend on this parameter...
    TIterator *pdfIter = theModel->_pdfLookupTable.createIterator();
    while (RooStringVar *pdfname = (RooStringVar*)pdfIter->Next()) {
      // If this PDF recipe depends on the split parameter, 
      // add the splitting category to the list of this model's splitting categories...
      MLPdfRecipe *rec = (MLPdfRecipe*)_recipeList->FindObject(pdfname->getVal());
      if ( rec->getProtoParams().find(arg->GetName()) != 0) {
	
	RooAbsCategory* cat = static_cast<RooAbsCategory*>(_splitCatList.at(_splitArgList.index(arg))) ;
	if( !splitCatSet.find(cat->GetName()) ) {
	  if( cat->InheritsFrom("RooAbsCategoryLValue") ) 
	    splitCatSet.add(*cat);
	  else {
	    bool foundFundamental = false ;
	    TIterator *catIter = _splitCatList.createIterator();
	    TObject* fundcat ;
	    while( (fundcat = catIter->Next()) && !foundFundamental)
	      if( fundcat->InheritsFrom("RooAbsCategoryLValue") )
		foundFundamental = cat->dependsOn(*(RooAbsCategoryLValue*)fundcat) ;
	    if(!foundFundamental) {
	      // FIX ME: memory leak
	      splitCatSet.add(*cat->createFundamental());
	    }
	  }
	}
      }
    }
  }

  // Also, if this PDF recipe is a simultaneous PDF, add the category that the sim pdf depends on.
  TIterator *pdfIter = theModel->_pdfLookupTable.createIterator();
  while (RooStringVar *pdf = (RooStringVar*)pdfIter->Next()) {
    MLPdfRecipe *rec = (MLPdfRecipe*)_recipeList->FindObject(pdf->getVal());
    if (rec->_simPdfCat != 0) splitCatSet.add(*rec->_simPdfCat);
  }
  delete pdfIter;

  if (verbose) {
    cout << "Creating " << theModel->GetName() << "'s master splitting category from" << endl;
    splitCatSet.Print();
  }
  
  theModel->createMasterSplitCat(splitCatSet);

  // Create the yield variables...
  createYieldVariables(*theModel, verbose);

  // And build all the basic PDFs from the PdfRecipes
  // Loop over the protype pdf names in this model
  if (verbose) cout << "Finally building the PDFs" << endl;
  TIterator *recIter = theModel->_pdfLookupTable.createIterator();
  while (RooStringVar *pdfname = (RooStringVar*)recIter->Next()){

    // Get the PDF recipe
    MLPdfRecipe *rec = (MLPdfRecipe*)_recipeList->FindObject(pdfname->getVal());
    if (verbose) cout << "Will build PDF: " << pdfname->GetName() << ", " << rec->GetName() <<endl;

    // If split, we need to loop over all splitting states...
    if (theModel->isSplit()) {
      TIterator *mcIter = theModel->_masterSplitCat->typeIterator();
      while (RooCatType *mcState = (RooCatType*)mcIter->Next()){

	// Set the master splitting category...
	theModel->_masterSplitCat->setLabel(mcState->GetName());

	if (_pdfList.find(rec->getSplitName()) == 0) {
	  RooAbsPdf *pdf = buildSinglePdf(*theModel, *rec, verbose);
	  _pdfList.add(*pdf);
	}
      }
      delete mcIter;
    } else {
      // If not split, just build the pdf...
      if (_pdfList.find(rec->getSplitName()) == 0) {
	RooAbsPdf *pdf = buildSinglePdf(*theModel, *rec, verbose);
	_pdfList.add(*pdf);
      }
    }// end if is/isn't split
  } //end pdfname loop
  delete recIter;
  
  // And now put all these pdfs into a...
  // RooSimPdf of RooAddPdfs of RooProdPdfs... bleh.
  RooAbsPdf *simPdf = buildSimPdf(*theModel);
  return simPdf;
  
}

void MLFit::splitPdfRecipe(MLPdfRecipe &theRecipe, Bool_t verbose)
{
  // This takes a PDF recipe and splits it by the categories in _splitCatList.
  // That means we set its _masterSplitCat correctly.

  if (verbose) cout << "Splitting PDF recipe: " << theRecipe.GetName() << endl;

  // First off, if this has a resolution function, we must split that
  if (theRecipe.isConvoluted()) {
    if (verbose) cout << "PDF Recipe is convoluted.  Will first split convolution function." << endl;
    splitPdfRecipe(*theRecipe._rf,verbose);
  }
  
  // First make a list of splitting categories that actually split this PDF (taken from _splitCatList)
  RooArgSet splitCats;
  TIterator *parIter = theRecipe._protoParams.createIterator();
  while (RooAbsArg *par = (RooAbsArg*)parIter->Next()){
    if (_splitArgList.find(par->GetName())){
      RooAbsCategory *cat = (RooAbsCategory*)_splitCatList.at(_splitArgList.index(par));
      splitCats.add(*cat);
    }
  }
  delete parIter;

  // If this is a simPdf, we need to add its sim category to the list.
  if (theRecipe._simPdfCat != 0) splitCats.add(*theRecipe._simPdfCat);

  // If it has a resolution function, we need to add the categories that split the
  // RF...  Note that this function must already have been called on the RF or we won't
  // know the RF is split yet!
  if (theRecipe.isConvoluted()) {
    if (theRecipe._rf->isSplit()) splitCats.add(theRecipe._rf->_masterSplitCat->inputCatList());
  }
  
  // So now if the PDF is, in fact, split by some categories, give it a master splitting category...
  if (splitCats.getSize() != 0) theRecipe.createMasterSplitCat(splitCats);

  if (verbose) {
    cout << "Will split recipe: " << theRecipe.GetName() << " by " << endl;
    splitCats.Print();
  }
  
  // Finally, if this is a simultaneous pdf, we have to split all the component pdf recipes...
  if (theRecipe._simPdfCat == 0) return;
  
  if (verbose) cout << theRecipe.GetName() << " is a simultaneous PDF, we must split each component..." << endl;
  TIterator *recIter = theRecipe._simPdfRecipeList.MakeIterator();
  while (MLPdfRecipe *rec = (MLPdfRecipe*)recIter->Next()) splitPdfRecipe(*rec);
  delete recIter;
}

void MLFit::getSplitParameters(const MLModel &theModel, const MLPdfRecipe &theRecipe, RooArgList &splitpars) const
{
  // If you have a PDF that's split and you want a list of all the (maybe split, maybe not) parameters
  // use this.  Pass it the model and recipe and an empty RooArgList that I'll fill.
  if (splitpars.getSize() != 0) {
    cout << GetName() << ": You gave getSplitParameters a non-empty list to fill!" << endl;
    return;
  }
  
  if (theModel.isSplit()){
    TIterator *argIter = theRecipe._protoParams.createIterator();
    while(RooAbsReal *arg = (RooAbsReal*)argIter->Next()){
      RooAbsReal *splitarg(0);
      splitarg = getSplitParameter(arg);
      
      if (splitarg == NULL) {
	cout << GetName() << ": Couldn't locate arg " << arg->GetName() << endl;
	return;
      }
      splitpars.add(*splitarg);
    }
    delete argIter;
  } else {
    // It's not split, so just copy the "prototype parameters"
    TIterator *argIter = theRecipe._protoParams.createIterator();      
    while(RooAbsReal *arg = (RooAbsReal*)argIter->Next()){
      RooAbsReal *splitarg = getSplitParameter(arg);
      if (splitarg == NULL) {
	cout << GetName() << ": Couldn't locate arg " << arg->GetName() << endl;
	return;
      }
      splitpars.add(*splitarg);
    }
    delete argIter;
  }
  
}

RooAbsPdf* MLFit::buildSinglePdf(const MLModel &theModel, const MLPdfRecipe &theRecipe, Bool_t verbose)
{
  // This builds a single PDF given a recipe (and a model that the recipe belongs to).

  TString newname = theRecipe.getSplitName();
  if (verbose) cout << "buildSinglePdf: " << newname << endl;

  // Get the parameters for this PDF (which may or may not be split)
  RooArgList newparameters;
  getSplitParameters(theModel, theRecipe, newparameters);
  
  // If this is a normal (non-simultaneous) pdf, just build it and return
  if (theRecipe._simPdfCat == 0 && !theRecipe.isConvoluted()) {
    if (verbose) cout << "Building normal single PDF: " << newname << "..." << endl;
    RooAbsPdf *pdf = MLPdfFactory::makePdf(theRecipe._type, newname, theRecipe._obs, newparameters, theRecipe._args);
    if (verbose) cout << "                                                       ... Done!" << endl;
    return pdf;
  }
  
  // If it's a convoluted pdf, we have to build its RF first.
  if (theRecipe._simPdfCat == 0 && theRecipe.isConvoluted()) {
    MLPdfRecipe *rfRec = (MLPdfRecipe*)_rfRecipeList.FindObject(theRecipe._rf->GetName());
    if (rfRec == NULL) {
      cout << GetName() << ": Trying to build " << newname << ", but its resolution function " << theRecipe._rf->GetName() << " can't be found in" << endl;
      _rfRecipeList.Print();
      return 0;
    }
    TString rfname = rfRec->getSplitName();
    
    // Get its (split or unsplit) parameters.
    RooArgList rfparameters;    
    getSplitParameters(theModel, *rfRec, rfparameters);

    // And build it (if it isn't already built)
    RooResolutionModel *rf = (RooResolutionModel*)_rfSet.find(rfname);
    if (rf == 0) {
      if (verbose) cout << "  Building RF: " << rfname << endl;
      rf = MLPdfFactory::makeRF(rfRec->_type, rfname, rfRec->_obs, rfparameters, rfRec->_args);
      _rfSet.add(*rf);
    }
    // And return the convoluted pdf...
    if (verbose) cout << "  Building convoluted PDF: " << newname << endl;
    return MLPdfFactory::makeConvolutedPdf(theRecipe._type, newname, theRecipe._obs, newparameters, *rf, theRecipe._args);
  }
  
  
  // Uh oh, if we get this far, it's a simultaneous pdf.  
  MLPdfRecipe *componentRec = theRecipe.getSimComponent();
  
  // Now build this pdf (by calling this very function again!)
  RooAbsPdf *pdf = (RooAbsPdf*)_pdfList.find(componentRec->getSplitName());
  if (pdf == 0){
    if (verbose) cout << "  Will build component of sim pdf" << endl;
    pdf = buildSinglePdf(theModel, *componentRec, verbose);
    _pdfList.add(*pdf);
  }
     
  RooAbsPdf *simPdf = (RooAbsPdf*)pdf->clone(newname);
  simPdf->SetTitle(newname);
  return simPdf;
}



RooAbsPdf* MLFit::buildSimPdf(MLModel &theModel)
{
  // If it's not split, we don't have to build a RooSimultaneous PDF...
  if (!theModel.isSplit()) {
    RooAbsPdf *pdf = buildAddPdf(theModel);
    pdf->SetName(theModel.GetName());
    pdf->SetTitle(theModel.GetTitle());
    return pdf;
  }
  
  RooSimultaneous *simPdf = new RooSimultaneous(theModel.GetName(), theModel.GetTitle(), *theModel._masterSplitCat);
  TIterator *mcIter = theModel._masterSplitCat->typeIterator();
  while (RooCatType *mcState = (RooCatType*)mcIter->Next()){
    // Build the full PDF for this one state
    theModel._masterSplitCat->setLabel(mcState->GetName());
    RooAbsPdf *pdf(0);
    // HUUGE HACK for simultaneous BReco fit...
    if (TString(theModel.GetName()) == "BRecoTagMixFitXXX") {
      TString oldlabel = CatObservable("moriondCat")->getLabel();
      CatObservable("moriondCat")->setLabel(CatObservable("tagCatMN")->getLabel());
      pdf = buildAddPdf(theModel);
      CatObservable("moriondCat")->setLabel(oldlabel);
      simPdf->addPdf(*pdf,theModel._masterSplitCat->getLabel());    
    } else {
      pdf = buildAddPdf(theModel);
      simPdf->addPdf(*pdf,theModel._masterSplitCat->getLabel());    
    }
    
  }
  delete mcIter;
  _pdfList.add(*simPdf);
  return simPdf;
}

RooAbsPdf* MLFit::buildAddPdf(MLModel &theModel)
{

  // If there's only one species, no need to create a RooAddPdf...
  if (theModel.numSpecies() == 1) {
    return buildProdPdf(theModel, _speciesList.First()->GetName());
  }
  
  RooArgList addList;
  RooArgList coefficients;
  // Iterate over species to build the RooAddPdf
  TIterator *specIter = theModel._speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()){
    RooAbsPdf *pdf = buildProdPdf(theModel, spec->GetName());
    addList.add(*pdf);
    RooAbsArg *coef = getCoefficient(theModel, spec->GetName());
    if (coef == 0) return 0;
    coefficients.add(*getCoefficient(theModel, spec->GetName()));
  }
  delete specIter;

  //  RooAbsArg *lastcoef = coefficients.at(coefficients.getSize()-1);
  //coefficients.remove(*lastcoef);  
  //_parameterSet.remove(*lastcoef);
  TString name, title;
  if (!theModel.isSplit()) {
    name = theModel.GetName();
    title = theModel.GetTitle();
  } else {
    name  = TString(theModel.GetName()) + "_" + theModel._masterSplitCat->getLabel();
    title = TString(theModel.GetTitle()) + " (" + theModel._masterSplitCat->getLabel() + ")";
  }

  RooAddPdf *addPdf = new RooAddPdf(name, title, addList, coefficients);
  _pdfList.add(*addPdf);
  
  return addPdf;
}

RooAbsPdf* MLFit::buildProdPdf(MLModel &theModel, const char* species) 
{

  // First determine what we're going to call it.
  TString name, title;
  if (!theModel.isSplit()) {
    name = TString(theModel.GetName())+"_"+getSpecies(species)->GetName();
    title = TString(theModel.GetName())+" "+getSpecies(species)->GetName() + " PDF";
  } else {
    name = TString(theModel.GetName())+"_"+getSpecies(species)->GetName() + "_" + TString(theModel._masterSplitCat->getLabel());
    title = TString(theModel.GetName())+" "+getSpecies(species)->GetName() + " (" + TString(theModel._masterSplitCat->getLabel()) + " )";
  }

  // A list of PDFs to make the RooProdPdf out of...
  //RooLinkedList prodList;
  RooArgList prodList;

  // Iterate over PDFs (recipes) for this species (mes, de, etc...) and create the RooProdPdf
  TIterator *recIter = getRecipeListBySpecies(theModel, species)->MakeIterator();
  while (MLPdfRecipe *rec = (MLPdfRecipe*)recIter->Next()){
    RooArgList newparameters;
    
    // Loop over this PDF's prototype parameters, and if they're split, substitute the appropriate ones...
    TIterator *argIter = rec->_protoParams.createIterator();
    while(RooAbsReal *arg = (RooAbsReal*)argIter->Next()){
      newparameters.add(*getSplitParameter(arg));
    }
    delete argIter;

    // Build the PDF
    RooAbsPdf *pdf = (RooAbsPdf*)_pdfList.find(rec->getSplitName());

    // If there's only one pdf recipe, there's no need to make a prod pdf.  Just make the
    // one pdf and return.
    if (getRecipeListBySpecies(theModel, species)->GetSize() == 1) {
      delete recIter;
      RooAbsPdf *prodpdf = (RooAbsPdf*)_pdfList.find(name);
      if (prodpdf == 0) {
	prodpdf = (RooAbsPdf*)pdf->clone(name);
	prodpdf->SetTitle(title);
	_pdfList.add(*prodpdf);
      }
      return prodpdf;
    }
    
    prodList.add(*pdf);

    /*
    if (rec->isPartial()) {
      prodList.Add(Partial(*pdf,rec->_partialList).Clone());
    } else {
      prodList.Add(Full(*pdf).Clone());
    }
    */
  }//end recipe loop
  delete recIter;
  
  RooProdPdf *pdf = new RooProdPdf(name,title,prodList);
  _pdfList.add(*pdf);
  return pdf;
}

void MLFit::defineYield(const char* species, const char* formula, RooArgList args) 
{
  MLSpecies *spec = (MLSpecies*)getSpecies(species);
  if (spec == 0) return;
  spec->setYieldVar(formula, args);
}

void MLFit::fitWithEff(TString species1, TString species2, TString namebase)
{

  if (getSpecies(species1) == 0) return;
  if (getSpecies(species2) == 0) return;

  RooRealVar *total = new RooRealVar("N_"+namebase,"N_"+namebase, 0);
  RooRealVar *eff  = new RooRealVar("eff_"+namebase,"eff_"+namebase,0);

  defineYield(species1, "@0*@1", RooArgList(*total, *eff));
  defineYield(species2, "@0*(1.-@1)", RooArgList(*total, *eff));

}

void MLFit::fitAsymmetry(TString species1, TString species2, TString namebase) 
{
  // Since you often want to fit for, say Total(Kpi) and the asymmetry between K+pi- and K-pi+, then
  // this is a function that redefines these yield variables.  You could do this in the macro, but
  // since it's so common, I'll make a special function for it.

  if (getSpecies(species1) == 0) return;
  if (getSpecies(species2) == 0) return;
  
  RooRealVar *total = new RooRealVar("N_"+namebase+"_total","N_"+namebase+"_total", 0);
  RooRealVar *asym  = new RooRealVar("A_"+namebase,"A_"+namebase,0);

  defineYield(species1, "@0*(0.5+@1/2.0)", RooArgList(*total, *asym));
  defineYield(species2, "@0*(0.5-@1/2.0)", RooArgList(*total, *asym));
  
}

void MLFit::fitInclusiveRatio(const MLStrList & species, TString namebase, int firstbin) {
  // this assumes that we start from >= 1j
  RooRealVar *total = new RooRealVar("N_"+namebase+"_total","N_"+namebase+"_total", 0);
  RooRealVar *alphaR = new RooRealVar("alphaR_"+namebase,"alphaR_"+namebase, 1);
  RooRealVar *betaR = new RooRealVar("betaR_"+namebase,"betaR_"+namebase, 0);

  int Nspecies = species.GetSize();
  TString scale("@0");
  for(int i=0; i<Nspecies; i++) {
    TString spec = ((TObjString*)species.At(i))->GetString();
    if (getSpecies(spec.Data()) == 0) return;
    TString addScale(""); 
    char nj[256];
    if(i < Nspecies-1) {
      sprintf(nj,"*%i)",i+firstbin);
      addScale = TString("-")+scale+TString("/(@1+@2")+TString(nj);
    }
    TString thisScale = scale+addScale;
    std::cout << thisScale.Data() << std::endl;
    defineYield(spec.Data(), thisScale.Data(), RooArgList(*total, *alphaR, *betaR));
    sprintf(nj,"*%i)",i+firstbin);
    scale = scale+TString("/(@1+@2")+TString(nj);
  }
}

void MLFit::fitInclusiveRatioPoly(const MLStrList & species, TString namebase, int firstbin, TArrayD efficiency) {
  // this assumes that we start from >= 1j
  RooRealVar *total = new RooRealVar("N_"+namebase+"_total","N_"+namebase+"_total", 0);
  RooArgList vars(*total);
  for(int i =0 ;i <5 ;i++){
    char numstring[10];
    sprintf(numstring,"%d",i);
    TString num(numstring);
    RooRealVar *alphaR = new RooRealVar("alpha"+num+"_"+namebase,"alpha"+num+"_"+namebase, 1);
    vars.add(*alphaR);
  }
  
  int Nspecies = species.GetSize();
  TString scale("@0");
  for(int i=0; i<Nspecies; i++) {
    TString spec = ((TObjString*)species.At(i))->GetString();
    if (getSpecies(spec.Data()) == 0) return;
    TString addScale(""); 
    char nj[256];
    char effRatio[200];
    if(i < Nspecies-1) {
      sprintf(nj,"/(@1+@2*%d +@3*%d*%d + @4*%d*%d*%d + @5*%d*%d*%d*%d)",i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin);
      addScale = TString("-")+scale+TString(nj);
      float thisEff = 1.0;
      float normEff = 1.0;
      if(efficiency.GetSize()!=0) {
        cout << "USING EXTERNAL EFFICIENCIES" << endl;
        thisEff = (float)efficiency.GetAt(firstbin+i+1);
        normEff = (float)efficiency.GetAt(firstbin);
        sprintf(effRatio,"*%f/%f",thisEff,normEff);
      } else { sprintf(effRatio, "*1"); }
    }
    TString thisScale = TString("(")+scale+addScale+TString(")")+TString(effRatio);
    std::cout << "species i = "<< i << " : " << thisScale.Data() << std::endl;
    defineYield(spec.Data(), thisScale.Data(), vars);
    sprintf(nj,"/(@1+@2*%d +@3*%d*%d + @4*%d*%d*%d + @5*%d*%d*%d*%d)",i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin);
    scale = scale+TString(nj);
  }
}

void MLFit::fitInclusiveRatioPolySUSY(const MLStrList & species, TString namebase, int firstbin, RooRealVar* sigyield) {
  // this assumes that we start from >= 1j
  RooRealVar *total = new RooRealVar("N_"+namebase+"_total","N_"+namebase+"_total", 0);
  RooArgList vars(*total);
  for(int i =0 ;i <5 ;i++){
    char numstring[10];
    sprintf(numstring,"%d",i);
    TString num(numstring);
    RooRealVar *alphaR = new RooRealVar("alpha"+num+"_"+namebase,"alpha"+num+"_"+namebase, 1);
    vars.add(*alphaR);
  }
  vars.add(*sigyield);

  int Nspecies = species.GetSize();
  TString scale("@0");
  for(int i=0; i<Nspecies; i++) {
    TString spec = ((TObjString*)species.At(i))->GetString();
    if (getSpecies(spec.Data()) == 0) return;
    TString addScale(""); 
    char nj[256];
    if(i < Nspecies-1) {
      sprintf(nj,"/(@1+@2*%d +@3*%d*%d + @4*%d*%d*%d + @5*%d*%d*%d*%d)",i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin);
      addScale = TString("-")+scale+TString(nj);
    }
    TString thisScale = scale+addScale;
    if( i == Nspecies-1 )
      thisScale +=" + @6";
    std::cout << thisScale.Data() << std::endl;
    defineYield(spec.Data(), thisScale.Data(), vars);
    sprintf(nj,"/(@1+@2*%d +@3*%d*%d + @4*%d*%d*%d + @5*%d*%d*%d*%d)",i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin);
    scale = scale+TString(nj);
  }

}
void MLFit::fitInclusiveRatioPolyUnfold(const MLStrList & species, TString namebase, TMatrix unfma, TString firstlabel, TArrayD efficiency) {
  // this assumes that we start from >= 1j
  RooRealVar *total = new RooRealVar("N_"+namebase+"_total","N_"+namebase+"_total", 0);
  RooArgList vars(*total);
  for(int i =0 ;i <5 ;i++){
    char numstring[10];
    sprintf(numstring,"%d",i);
    TString num(numstring);
    RooRealVar *alphaR = new RooRealVar("alpha"+num+"_"+namebase,"alpha"+num+"_"+namebase, 1);
    vars.add(*alphaR);
  }

  int firstbin=1;
  
  RooRealVar *zeroyield = new RooRealVar("N_excl0", "N_excl0",0);
  vars.add(*zeroyield);

  MLStrList exclgen;
  exclgen.Add("@6"); //yield for zero jets

  //from inclusive yields to exclusive yields
  int Nspecies = species.GetSize();
  TString scale("@0");
  for(int i=0; i<Nspecies; i++) {
    TString spec = ((TObjString*)species.At(i))->GetString();
    if (getSpecies(spec.Data()) == 0) return;
    TString addScale(""); 
    char nj[256];
    char effRatio[200];
    if(i < Nspecies-1) {
      sprintf(nj,"/(@1+@2*%d +@3*%d*%d + @4*%d*%d*%d + @5*%d*%d*%d*%d)",i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin);
      addScale = TString("-")+scale+TString(nj);
    }
    float thisEff = 1.0;
    float normEff = 1.0;
    if(efficiency.GetSize()!=0) {
      cout << "USING EXTERNAL EFFICIENCIES" << endl;
      thisEff = (float)efficiency.GetAt(i+1);
      normEff = (float)efficiency.GetAt(0); // when unfolding, firstbin is always 0
      sprintf(effRatio,"*%f/%f",thisEff,normEff);
    } else { sprintf(effRatio, "*1"); }

    TString thisScale = TString("(")+scale+addScale+TString(")")+TString(effRatio);
    std::cout << thisScale.Data() << std::endl;
    exclgen.Add(thisScale.Data());
    sprintf(nj,"/(@1+@2*%d +@3*%d*%d + @4*%d*%d*%d + @5*%d*%d*%d*%d)",i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin,i+firstbin);
    scale = scale+TString(nj);
  }

  //from unfolded exlcusive yields to folded exclusive yields
  for(int i=0; i<Nspecies+1; i++) {
    TString spec("");
    if(i==0) spec = firstlabel;
    else spec = ((TObjString*)species.At(i-1))->GetString();
    if (getSpecies(spec.Data()) == 0) return;
    //prepare matrix string
    TString foldstring("");
    for(int j=0; j<Nspecies + 1; j++){
      foldstring += unfma(i,j);// TMatrix(row,col)
      foldstring +=  "*( ";
      foldstring +=  ((TObjString*)exclgen.At(j))->GetString();
      foldstring += " ) ";
      if( j <Nspecies )foldstring += "+";
    }
    
    std::cout << i << ":  " <<foldstring << std::endl;
    defineYield(spec.Data(),foldstring.Data(), vars);
  }

}


void MLFit::fitInclusive(const MLStrList & species, TString namebase, int firstbin) {
  // this assumes that we start from >= 1j
  RooArgList* N_incl=new RooArgList("InclYields");
  int Nspecies = species.GetSize();
  for(int i=0; i<Nspecies; i++) {
    char nj[256];
    sprintf(nj,"%ij",i+firstbin);
    RooRealVar *total = new RooRealVar( (TString("N_")+namebase+TString(nj)).Data(), 
					(TString("N_")+namebase+TString(nj)).Data(), 0);
    N_incl->add(*total);
  } 
  for(int i=0; i<Nspecies; i++) {
    TString spec = ((TObjString*)species.At(i))->GetString();
    if (getSpecies(spec.Data()) == 0) return;
    if(i < Nspecies-1)
      defineYield(spec.Data(), "@0-@1", RooArgList(*N_incl->at(i),*N_incl->at(i+1)));
    else
      defineYield(spec.Data(), "@0", RooArgList(*N_incl->at(i)));	
  }

}

void MLFit::fitInclusiveUnfold(const MLStrList & species, TString namebase, TMatrix unfma) {
  // unfoldding needs to include the NJet=0 bin to account for migrations
  RooArgList* N_incl=new RooArgList("InclYields");
  int Nspecies = species.GetSize();
  for(int i=0; i<Nspecies; i++) {
    char nj[256];
    sprintf(nj,"%ij",i);
    RooRealVar *total = new RooRealVar( (TString("N_")+namebase+TString(nj)).Data(), 
					(TString("N_")+namebase+TString(nj)).Data(), 0);
    N_incl->add(*total);
  } 
  N_incl->Print("v");
  for(int i=0; i<Nspecies; i++) {
    TString spec = ((TObjString*)species.At(i))->GetString();
    if (getSpecies(spec.Data()) == 0) return;
    //prepare matrix string
    TString foldstring("");
    for(int j=0; j<Nspecies; j++){
      if(j < Nspecies-1){
	foldstring += unfma(i,j);// TMatrix(row,col)
	foldstring +="*(@";
	foldstring +=j;
	foldstring +=" - @";
	foldstring +=j+1;
	foldstring +=" ) + ";
      }
      else{
	foldstring += unfma(i,j);
	foldstring +=  "*(@";
	foldstring += j ;
	foldstring += " )";
      }
    }
    std::cout << foldstring << std::endl;
    defineYield(spec.Data(),foldstring.Data(), *N_incl);
  }

}



void MLFit::createYieldVariables(MLModel &theModel, Bool_t verbose)
{
  // If there's just one species, just return...
  //if (theModel.numSpecies() == 1) return;

  if (verbose) cout << "Creating yield variables" << endl;
  // Okay, each PDF has a coefficient in front of it.
  // In this case, we're assuming we can factorize the coefficient into the product of some yield variable, say, N_sig, and
  // a bunch of "tagging fractions" e_sig,cat1.  One "tagging fraction" for each
  // splitting category.  With, say e_1+e_2+e_3=1.  Remember, N_sig itself might be a function of some other parameters, and
  // to make things worse, some fits (can we say "breco"?!) might want *those* parameters split, also.  This makes for a lot
  // of complicated things going on here.  See if you can follow it :)

  // First create the yields (N_sig).  We have up to n_species*n_categories different yield parameters.
  
  
  TIterator *specIter = theModel._speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()){
    
    // You only want to split the yield parameter if it depends on some other parameter that is split.
    if (spec->getProtoParams().overlaps(_splitArgList)) {
      TIterator *mcIter = theModel._masterSplitCat->typeIterator();
      while (RooCatType *mcState = (RooCatType*)mcIter->Next()){   
	// Set the master splitting category...
	theModel._masterSplitCat->setLabel(mcState->GetName());
	RooArgList newparameters;
	RooArgList protoparams(spec->getProtoParams());
	TIterator *argIter = protoparams.createIterator();
	while(RooRealVar* arg = (RooRealVar*)argIter->Next()){
	  newparameters.add(*getSplitParameter(arg));
	}
	delete argIter;
	TString newname  = TString("N_")+spec->GetName()+"_"+theModel._masterSplitCat->getLabel();
	TString newtitle = TString("N ")+spec->GetTitle()+" ("+theModel._masterSplitCat->getLabel()+")";
	RooAbsReal *yieldvar = (RooAbsReal*)_fracList.find(newname);
	if (yieldvar == 0){
	  yieldvar = spec->buildYieldVar(newparameters);
	  yieldvar->SetName(newname);
	  yieldvar->SetTitle(newtitle);
	  _fracList.add(*yieldvar);
	  if (yieldvar->isFundamental()) _parameterSet.add(*yieldvar,kTRUE);
	}
      }
      delete mcIter;
    } else {
      // If not split, we just need one yield var per species
      TString name = TString("N_")+spec->GetName();
      RooAbsReal *yieldvar = (RooAbsReal*)_fracList.find(name);
      if (yieldvar == 0){
	yieldvar = spec->buildYieldVar(spec->getProtoParams());
	_fracList.add(*yieldvar);
	if (yieldvar->isFundamental()) _parameterSet.add(*yieldvar,kTRUE);
      }     
    }
  }

  // Create the fractions and "tagging fractions" (e_sig,cat1).  One "tagging fraction" for each species-splitcat.
  if (theModel.isSplit()) {
    specIter->Reset();
    while (MLSpecies *spec = (MLSpecies*)specIter->Next()){
      // The "prototype" "tagging fraction" variable
      RooAbsReal *evar = spec->getFracVar();

      // Iterate over the different splitting categories that make up the master splitting category.
      TIterator *catIter = theModel._masterSplitCat->inputCatList().createIterator();
      while (RooAbsCategory* catvar = (RooAbsCategory*)catIter->Next()) {

	// If we don't want "tagging fractions" for this category, just skip it.
	if (_noSplitFracSet.find(catvar->GetName()) != NULL) continue;
	
	// The final "tagging fraction" is a function of the others (like 1-e_1-e_2-e_3 or something).
	// So we'll slowly build up the function as a TString, and keep track of the tagging fractions
	// that go into the function.
	Int_t numtypes = catvar->numTypes();
	Int_t ntype = 0;
	TString lastFracFunc;
	RooArgList lastFracList;

	// And now iterate over the different state types of this splitting category.
	TIterator *typeIter = catvar->typeIterator();
	while (RooCatType *type = (RooCatType*)typeIter->Next()){
	  RooAbsReal *epsilon;
	  TString ename = TString(evar->GetName())+"_"+type->GetName();
	  TString etitle = TString(evar->GetTitle())+" ("+type->GetName()+")";

	  // If this is the last type in this category, it's a function (i.e. 1-e_1-e_2-e_3)
	  ntype++;
	  if (ntype == numtypes) {
	    lastFracFunc.Chop();
	    lastFracFunc.Prepend("1-");
	    if (_fracList.find(ename) == 0){
	      epsilon = new RooFormulaVar(ename,etitle,lastFracFunc,lastFracList);
	      _fracList.add(*epsilon);
	    }
	  } else {
	    epsilon = new RooRealVar(ename,etitle,1);
	    if (_fracList.find(ename)==0) _fracList.add(*epsilon);
	    _parameterSet.add(*epsilon,kTRUE);
	    lastFracFunc.Append(ename+"-");
	    lastFracList.add(*epsilon);
	  }
	}// end type loop
	delete typeIter;
      }// end cat loop
      delete catIter;
    }//end if split
  }// end spec loop
  delete specIter;
  
  // Finally, build the coefficients, which are a product of the yield variable and the 
  // "tagging fractions" (N_sig*e_sig,lep*e_sig,good)
  const RooArgSet &masterCatList = theModel._masterSplitCat->inputCatList();
  TIterator *mcIter = theModel._masterSplitCat->typeIterator();
  while (RooCatType *mcState = (RooCatType*)mcIter->Next()){
    theModel._masterSplitCat->setLabel(mcState->GetName());
    
    specIter = theModel._speciesList.MakeIterator();
    while (MLSpecies *spec = (MLSpecies*)specIter->Next()){  
      // Get the yield variable for this splitting state.  If the yield variable isn't split
      // (i.e. doesn't depend on any of the splitting args) it's name is simple.  Otherwise,
      // append the splitting state.
      TString yieldvarname;
      if (spec->getProtoParams().overlaps(_splitArgList)) {
	yieldvarname = TString("N_")+spec->GetName()+"_"+mcState->GetName();
      } else {
	yieldvarname = TString("N_")+spec->GetName();
      }
      RooAbsReal *yieldvar = (RooAbsReal*)_fracList.find(yieldvarname);

      // The coefficient is a function, so let's start building up the function
      // as a TString, and keep track of the parameters that go into it.
      RooAbsReal *evar = spec->getFracVar();
      TString tagfracFunc = yieldvar->GetName();
      RooArgList tagfracList;
      tagfracList.add(*yieldvar);
      TString statename = "{";
      TIterator *catIter = masterCatList.createIterator();
      while (RooAbsCategory *catvar = (RooAbsCategory*)catIter->Next()) {

	// Don't make "tagging fractions" for this category if we've been specifically asked
	// not to.
	statename.Append(TString(catvar->getLabel())+";");
	if (_noSplitFracSet.find(catvar->GetName()) != NULL) continue;
	if ((TString(theModel.GetName()) == "BRecoTagMixFit") && TString(catvar->GetName()) == "moriondCat") continue;
	TString ename = TString(evar->GetName())+"_"+catvar->getLabel();
	tagfracFunc.Append(TString("*")+ename);
	RooAbsArg *epsilon = _fracList.find(ename);
	tagfracList.add(*epsilon);
      }
      delete catIter;
      statename.Chop();
      statename.Append("}");
      TString cname =  TString("coef_") + spec->GetName() + "_" + statename;
      TString ctitle = TString("Coef. ") + spec->GetName() + " (" + statename + ")";
      if (_fracList.find(cname) == 0) {
	RooFormulaVar *coef = new RooFormulaVar(cname,ctitle,tagfracFunc,tagfracList);
	_fracList.add(*coef);
      }
    }
    delete specIter;
  }
  delete mcIter;

  // One last thing:  If the user wanted to blind a yield parameter, they already called the blindUniform() function,
  // and the yield parameters didn't exist.  So we added the name to a list of _blindedYieldNames.  Now that we've created the
  // yields variables, we can blind them.
  TIterator *blindIter = _blindedYieldNames.createIterator();
  while (RooRealVar *blindname = dynamic_cast<RooRealVar*>(blindIter->Next())){
    RooRealVar *par = dynamic_cast<RooRealVar*>(_parameterSet.find(blindname->GetName()));
    if (par == NULL){
      cout << "Trying to blind " << blindname->GetName() << ", but there is no parameter of that type in " << endl;
      _parameterSet.Print();
      continue;
    }
    blindUniform(*par,blindname->GetTitle(),blindname->getVal());
  }
  delete blindIter;

}


void MLFit::createSingleSpeciesFractions(MLModel &theModel, TString namesuffix, TString titlesuffix)
{
  // Creates the fraction variables for fit using fractions instead of yields...
  
  MLSpecies *lastSpec = (MLSpecies*)theModel._speciesList.At(theModel._speciesList.GetSize()-1);
  TList firstSpecList;
  firstSpecList.AddAll(&theModel._speciesList);
  firstSpecList.Remove(lastSpec);
  RooArgList fundFrac;
  TString lastFormula;
  TIterator *specIter = firstSpecList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()){
    TString name  = "F_" + TString(spec->GetName()) + namesuffix;  
    TString title = "F_" + TString(spec->GetName()) + titlesuffix;
    if (_fracList.find(name) == 0) {
      
      RooRealVar *frac = new RooRealVar(name, title, 0);
      fundFrac.add(*frac);
      _fracList.add(*frac);
      _parameterSet.add(*frac);
    }
    lastFormula.Append(name + "-");
  }

  TString name = "F_" + TString(lastSpec->GetName()) + namesuffix;  
  TString title = "F_" + TString(lastSpec->GetName()) + titlesuffix;
  if (_fracList.find(name) == 0) {
    lastFormula.Chop();
    lastFormula.Prepend("1-");
    RooFormulaVar *frac = new RooFormulaVar(name, title, lastFormula, fundFrac);
    _fracList.add(*frac);
  }
  
}
//=======================================
// Misc...
//=======================================
void MLFit::printFloating() 
{
  _parameterSet.selectByAttrib("Constant",kFALSE)->Print("V");
}
void MLFit::printConstants() 
{
  _parameterSet.selectByAttrib("Constant",kTRUE)->Print("V");
}

void MLFit::Print(Option_t *option) const
{
  cout << "MLFit: " << fName << "   " << fTitle << endl;
  
  cout << "Species:" << endl;
  _speciesList.Print();
  
  cout << "Fractions:" << endl;
  _fracList.Print();
  
  cout << "Parameters:" << endl;
  _parameterSet.Print("V");
  
  cout << "PDFs:" << endl;
  _pdfList.Print();
  
  cout << "Models:" << endl;
  _modelList.Print();
}

void MLFit::initialize(const char* fileName) 
{
  // This loads in all parameters and observables from a config file
  RooArgSet paramsettmp, paramset;
  //paramset.add(_speciesList);
  paramsettmp.add(_fracList);
  paramsettmp.add(_parameterSet,kTRUE);
  TIterator *partmpIter = paramsettmp.createIterator();
  while (RooAbsArg* par = (RooAbsArg*)partmpIter->Next()){
    if (!par->IsA()->InheritsFrom("RooAbsHiddenReal")) paramset.add(*par);
  }
  delete partmpIter;
  paramset.add(_blindParSet);
  paramset.readFromFile(fileName, "READ", "Parameters");
  
  // Now check that _all_ the parameters were read in.
  TIterator *parIter = paramset.createIterator();
  while (RooAbsArg* par = (RooAbsArg*)parIter->Next()) {
    if (!par->getAttribute("READ") && par->isFundamental()) cout << "*** WARNING: " << par->GetName() << " was not initialized from file." << endl;
  }
  delete parIter;

  // Now read in the flat file variables...
  initializeFlatFileVars(fileName);

}

void MLFit::initializeFlatFileVars(TString fileName)
{
  // Initialize the flat file variables from a file.

  RooArgSet(_observables).readFromFile(fileName, "READ", "Observables");

  // Verify that everything was read in that was supposed to be.
  TIterator *obsIter = _observables.createIterator();
  while (RooAbsArg* obs = (RooAbsArg*)obsIter->Next()) {
    if (!obs->getAttribute("READ") && !obs->IsA()->InheritsFrom("RooAbsCategory")) cout << "*** WARNING: " << obs->GetName() << " was not initialized from file." << endl;
  }
  delete obsIter;

}


void MLFit::writeConfigFile(const char* fileName)
{
  //
  // Write out a configuration file (named fileName) containing all the
  // values, limits, etc. of the observables and parameters used
  // in this fit.  It outputs the current values of these variables.
  //
  ofstream ofs(fileName) ;
  if (ofs.fail()) {
    cout << "MLFit: error opening file " << fileName << endl ;
    return;
  }
  ofs << endl << "   [Parameters]" << endl;
  //RooArgSet(_speciesList).writeToStream(ofs, kFALSE);
  _fracList.writeToStream(ofs, kFALSE);
  _parameterSet.writeToStream(ofs, kFALSE);
  ofs << endl << "   [Observables]" << endl;
  RooArgSet obset(getObservables());
  obset.writeToStream(ofs,kTRUE);
}

RooArgSet MLFit::getParameters() const 
{
  RooArgSet theset(_parameterSet);
  theset.add(_blindParSet);
  return theset;
}

void MLFit::smearConstantParameters(const MLStrList &strList) const
{
  cout << ">>>>>>>>>>>>>>>>>>>" << endl;
  RooArgList arglist;
  TIterator *strIter = strList.MakeIterator();
  while (TObjString *str = (TObjString*)strIter->Next()) {
    RooRealVar *param = getRealPar(str->GetName());
    if (param == 0) continue;

    if (!param->hasError()) return;

    Double_t val=0;
    while (1) {
      val = RooRandom::randomGenerator()->Gaus( param->getVal(), param->getError() );
      if (val>param->getMin() && val<param->getMax()) break;
    }
    cout << "Smearing constant parameter " << param->GetName()
         << "  " << param->getVal() << " --> ";
    param->setVal( val);
    cout << param->getVal() << endl;

  }
  cout << "<<<<<<<<<<<<<<<<<<<<<<" << endl;
}

void MLFit::uniformVaryConstantParameters(const MLStrList &strList) const
{
  cout << ">>>>>>>>>>>>>>>>>>>" << endl;
  RooArgList arglist;
  TIterator *strIter = strList.MakeIterator();
  while (TObjString *str = (TObjString*)strIter->Next()) {
    RooRealVar *param = getRealPar(str->GetName());
    if (param == 0) continue;

    if (!param->hasError()) return;

    Double_t val=0;
    while (1) {
      // vary parameter between boundaries uniformly: error in this case has to assume the meaning of boundary
      val = RooRandom::randomGenerator()->Uniform( param->getVal() - param->getError(),  param->getVal() + param->getError());
      if (val>param->getMin() && val<param->getMax()) break;
    }
    cout << "Varying uniformly constant parameter " << param->GetName()
         << "  " << param->getVal() << " --> ";
    param->setVal( val);
    cout << param->getVal() << endl;

  }
  cout << "<<<<<<<<<<<<<<<<<<<<<<" << endl;
}

