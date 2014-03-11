
#include "MLFit/MLOptions.hh"
#include "RooCategory.h"
#include "RooStringVar.h"
#include "RooRealVar.h"

ClassImp(MLOptions);


MLOptions::MLOptions()
{
}

MLOptions::MLOptions(const MLOptions &opts) :
  _boolOptions(opts._boolOptions),
  _stringOptions(opts._stringOptions),
  _realOptions(opts._realOptions)
{
}

MLOptions& MLOptions::operator= (const MLOptions &other)
{
  // Check for self-assignment
  if (this == &other) return *this;
  
  _boolOptions.removeAll();
  _stringOptions.removeAll();
  _realOptions.removeAll();
  _boolOptions.add(other._boolOptions);
  _stringOptions.add(other._stringOptions);
  _realOptions.add(other._realOptions);
  return *this;
}


MLOptions::~MLOptions()
{
}

void MLOptions::addBoolOption(const char* name, const char* title, Bool_t def) 
{
  RooCategory *opt = new RooCategory(name, title);
  opt->defineType("kFALSE", 0);
  opt->defineType("kTRUE" , 1);
  if (def == kTRUE) {
    opt->setLabel("kTRUE");
  } else {
    opt->setLabel("kFALSE");
  }
  _boolOptions.add(*opt);
}

void MLOptions::addStringOption(const char* name, const char* title, const char* def)
{
  RooStringVar *opt = new RooStringVar(name, title, def);
  _stringOptions.add(*opt);
}

void MLOptions::addRealOption(const char* name, const char* title, double def)
{
  RooRealVar *opt = new RooRealVar(name, title, def);
  _realOptions.add(*opt);
}

void MLOptions::readEnv( Bool_t explicitMode )
{
  Print();
  // This reads all variables from environmental variables
  TIterator *boptIter = _boolOptions.createIterator();
  while (RooCategory *bopt = (RooCategory*)boptIter->Next()) {
    TString name = bopt->GetName();
    if( explicitMode ){
       TString val = getenv(name);
       // only environmenttal variables with values are considered
       if (val.Length() != 0){
	  if (val.Index("true")>=0 || val.Index("TRUE")>=0 )
	    bopt->setLabel("kTRUE");
	  else
	    bopt->setLabel("kFALSE");
       }
    }else{
       // If the env variable is set to anything, it's true
       if (getenv(name) != 0)
	 bopt->setLabel("kTRUE");
       else
	 bopt->setLabel("kFALSE");
    }
  }
  delete boptIter;

  Print();
  TIterator *strIter = _stringOptions.createIterator();
  while (RooStringVar *opt = (RooStringVar*)strIter->Next()) {
    TString name = opt->GetName();
    const char* val = getenv(name);
    if (val != 0) opt->setVal(val);
  }
  delete strIter;

  Print();
  TIterator *realIter = _realOptions.createIterator();
  while (RooRealVar *opt = (RooRealVar*)realIter->Next()){
    TString name = opt->GetName();
    const char *val = getenv(name);
    if (val != 0) opt->setVal(atof(val));
  }
  delete realIter;
}

void MLOptions::setBoolVal(TString opt, Bool_t value)
{
  RooCategory *arg = (RooCategory*)_boolOptions.find(opt);
  if (arg == NULL) {
    cout << "MLOptions Error: Can't find option " << opt << endl;
    assert(0) ;
    return;
  }
  if (value == kTRUE){
    arg->setLabel("kTRUE");
  } else {
    arg->setLabel("kFALSE");
  }
}


Bool_t MLOptions::getBoolVal(const char* opt) const
{
  RooCategory *arg = (RooCategory*)_boolOptions.find(opt);
  if (arg == NULL) {
    cout << "MLOptions Error: Can't find option " << opt << ".  Returning kFALSE" << endl;
    assert(0) ;
    return kFALSE;
  }
  
  if (*arg == "kTRUE") {
    return kTRUE;
  } else {
    return kFALSE;
  }
}

void MLOptions::setStringVal(TString opt, TString value)
{
  RooStringVar *arg = (RooStringVar*)_stringOptions.find(opt);
  if (arg == 0) {
    cout << "MLOptions Error: Can't find string option " << opt << endl;
    assert(0) ;
    return;
  }
  arg->setVal(value);
}


TString MLOptions::getStringVal(const char* opt) const
{
  RooStringVar *arg = (RooStringVar*)_stringOptions.find(opt);
  if (arg == 0) {
    cout << "MLOptions Error: Can't find string option " << opt << ".  Returning empty string." << endl;
    assert(0) ;
    return "";
  }
  return arg->getVal();
}

void MLOptions::setRealVal(TString opt, double value)
{
  RooRealVar *arg = (RooRealVar*)_realOptions.find(opt);
  if (arg == 0) {
    cout << "MLOptions Error: Can't find real-valued option " << opt << endl;
    assert(0) ;
    return;
  }
  arg->setVal(value);
}


double MLOptions::getRealVal(const char* opt) const
{
  RooRealVar *arg = (RooRealVar*)_realOptions.find(opt);
  if (arg == 0){
    cout << "MLOptions Error: Can't find real-valued option " << opt << ".  Returning 0." << endl;
    assert(0) ;
    return 0;
  }
  return arg->getVal();
}


void MLOptions::Print(Option_t *option) const 
{
   cout << "## Logical options: " <<endl;
   _boolOptions.Print("V");
   cout << "## Text options: " <<endl;
   _stringOptions.Print("V");
   cout << "## Float-point options: " <<endl;
   _realOptions.Print("V");
}


void MLOptions::setOptions(const MLOptions& rhs)
{
  TObject* obj ;
  TIterator *iter ;
  
  iter = rhs._boolOptions.createIterator();
  while ( (obj = iter->Next()) ) 
    setBoolVal( obj->GetName(), rhs.getBoolVal(obj->GetName()) ) ;
  delete iter ;

  iter = rhs._stringOptions.createIterator();
  while ( (obj = iter->Next()) ) 
    setStringVal( obj->GetName(), rhs.getStringVal(obj->GetName()) ) ;
  delete iter ;

  iter = rhs._realOptions.createIterator();
  while ( (obj = iter->Next()) ) 
    setRealVal( obj->GetName(), rhs.getRealVal(obj->GetName()) ) ;
  delete iter ;
}
