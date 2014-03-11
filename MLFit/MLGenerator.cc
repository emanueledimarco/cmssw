#include "MLFit/MLGenerator.hh"
#include "MLFit/MLModel.hh"
#include "MLFit/MLSpecies.hh"
#include <TFile.h>
#include <RooAbsPdf.h>
#include <RooCategory.h>
#include <RooRandom.h>
#include <algorithm>

ClassImp(MLGenerator);

MLGenerator::MLGenerator(const MLFit &theFit, const char* model, Bool_t doTruth) :
  _theFit((MLFit*)&theFit),
  _doTruth(doTruth)
{
  _theModel = theFit.getModel(model);
  _masterPdf = theFit.getPdf(model);

 // If you want to store the MC truth, create a category variable that shows which species was generated.
  if (_doTruth){
    _truth = new RooCategory("truth", "MC Truth");
    TIterator *specIter = _theModel->_speciesList.MakeIterator();
    while (MLSpecies *spec = (MLSpecies*)specIter->Next()) _truth->defineType(spec->GetName());
    delete specIter;
  }

  // Keep track of the actual number of events of each species that we wanted to generate.
  TIterator *specIter = _theModel->_speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()) {
    TString name = TString(spec->GetName())+"_actual";
    TString title = TString(spec->GetTitle())+" (Actual)";
    _extraParams.add(*(new RooRealVar(name, title, 0)));
  }
 
}

MLGenerator::~MLGenerator()
{
}

RooDataSet* MLGenerator::generate(const RooArgSet &whatVars, Int_t ngen, Bool_t verbose, TString species)
{
  // This generates a dataset from the MLFit object _theFit.
  // It generates each species and category separately, then adds them together at the end.
  // In the bizzare event that you want to generate ONLY one species, use the species argument.
  // Otherwise, I assume you want to generate some of all species.


  cout << "Will generate " << ngen << " events." << endl;
  // An array containing the number of events of each species to generate.  I.e. yieldSet[pipi:{Cat1}] = 12
  map<string,int> yieldSet; 
  
  calculateYields(ngen, yieldSet, verbose, species);
   cout << "In main" << endl;
   for (map<string,int>::iterator iter = yieldSet.begin(); iter != yieldSet.end(); iter++){
    cout << "Key: " << iter->first << ", value: " << iter->second << endl;
  }
  
  // Now generate each of the data sets.  One for each species/splitting category
  RooArgSet genVars(whatVars);
  genVars.add(*_truth);
  RooDataSet *theData = new RooDataSet("theData", "The Data", genVars);
  
  TIterator *specIter = _theModel->_speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next()) {
    if (_theModel->isSplit()) {
      TIterator *mcIter = _theModel->_masterSplitCat->typeIterator();
      while (RooCatType *mcState = (RooCatType*)mcIter->Next()){
	_theModel->_masterSplitCat->setLabel(mcState->GetName());
	string key = string(spec->GetName())+":"+string(mcState->GetName());
	//==Int_t n = (Int_t)((RooAbsReal*)(&yieldSet[key]))->getVal();
	int n = yieldSet[key];
	// Nick
	cout << "Should generate " << n << " events here..." << endl;
	generateComponent(*spec, whatVars, n, *theData, verbose);	
      }
      delete mcIter;
    } else {
      // not split...
      string key = spec->GetName();
      int n = yieldSet[key];

      Bool_t mockThisSpecies = kFALSE;
      TString specName = spec->GetName();
      TString fileName = TString("");
      TString dataName = TString("");
      TIterator *mockSpecIter = _mockSpecies.MakeIterator();
      TIterator *mockFilesIter = _mockFiles.MakeIterator();
      TIterator *mockDatasetIter = _mockDatasetNames.MakeIterator();
      while (TObjString *str = (TObjString*)mockSpecIter->Next()) {
        TObjString *filename = (TObjString*)mockFilesIter->Next();
        TObjString *dataname = (TObjString*)mockDatasetIter->Next();
        if (str->GetString().Contains(specName)) {
          mockThisSpecies = kTRUE;
          fileName = filename->GetString();
          dataName = dataname->GetString();
        }
      }

      if (mockThisSpecies) {
        // Emanuele
        cout << "Should mock " << n << " events from file " << fileName << " here..." << endl;
        mockComponent(fileName, dataName, *spec, whatVars, n, *theData, verbose);
      } else {
        // Nick
        cout << "Should generate " << n << " events here..." << endl;
        generateComponent(*spec, whatVars, n, *theData, verbose);
      }
    }
    
  }
  delete specIter;

  yieldSet.clear();

  return theData;
}

void MLGenerator::calculateYields(Int_t ngen, map<string,int> &yieldSet, Bool_t verbose, TString species)
{
  // First get the sum of all the coefficients.  For an extended likelihood fit, this should just
  // be Ntot
  Double_t sumCoef = 0;
  TIterator *specIter = _theModel->_speciesList.MakeIterator();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next())
  {
    // If the model is split, iterate over the splitting category states
    if (_theModel->isSplit()) {
      TIterator *mcIter = _theModel->_masterSplitCat->typeIterator();
      while (RooCatType *mcState = (RooCatType*)mcIter->Next()){
	_theModel->_masterSplitCat->setLabel(mcState->GetName());
	if (species == "" || (species == spec->GetName())) sumCoef += _theFit->getCoefficient(*_theModel, spec->GetName())->getVal();
      }
      delete mcIter;
    } else {
      // Not split..
      if (species == "" || (species == spec->GetName())) sumCoef += _theFit->getCoefficient(*_theModel, spec->GetName())->getVal();
    }    
  }
  
  if (verbose) cout << "MLGenerator: Sum of coefficients: " << sumCoef << endl;

  // Now determine how many of each species/splitcat to generate
  for (Int_t ievent = 0; ievent < ngen; ievent++){
    //cout << "New Event" << endl;
    Double_t r = RooRandom::uniform();
    Double_t max = 0;
    Double_t min = 0;
    specIter->Reset();
    while(MLSpecies *spec = (MLSpecies*)specIter->Next()){
      if (_theModel->isSplit()){
	TIterator *mcIter = _theModel->_masterSplitCat->typeIterator();
	while (RooCatType *mcState = (RooCatType*)mcIter->Next()){
	  _theModel->_masterSplitCat->setLabel(mcState->GetName());
	  string key = string(spec->GetName())+":"+string(mcState->GetName());
	  if (yieldSet.find(key) == yieldSet.end()) yieldSet[key] = 0;

	  // Skip all this if for some reason you want to generate one species (and this isn't the one).
	  if (species != "" && (species != spec->GetName())) continue;

	  // Fraction of events of this species in this splitting category
	  Double_t f = _theFit->getCoefficient(*_theModel, spec->GetName())->getVal()/sumCoef;
	  max += f;
	  //==if (r >= min && r < max) yieldSet[key] = ((RooAbsReal*)(&yieldSet[key]))->getVal() + 1;
	  if (r >= min && r < max) yieldSet[key]++;
	  min += f;
	}
	delete mcIter;
      } else {
	// Not split...
	string key = spec->GetName();
	if (yieldSet.find(key) == yieldSet.end()) yieldSet[key] = 0;
	// Skip all this if for some reason you want to generate one species (and this isn't the one).
	if (species != "" && (species != spec->GetName())) continue;

	Double_t f = _theFit->getCoefficient(*_theModel, spec->GetName())->getVal()/sumCoef;
	max += f;
	//==if (r >= min && r < max) yieldSet[key] = ((RooAbsReal*)(&yieldSet[key]))->getVal() + 1;
	if (r >= min && r < max) yieldSet[key]++;
	min += f;
      }//end if split
    }// end spec loop
  }//endfor

  
  // We'll want to store the number of each species that was actually generated, so
  // let's make the variables to do that.
  specIter->Reset();
  while (MLSpecies *spec = (MLSpecies*)specIter->Next())
  {
    TString name = TString(spec->GetName())+"_actual";
    RooRealVar *specsum = (RooRealVar*)_extraParams.find(name);
    // If the model is split, iterate over the splitting category states
    if (_theModel->isSplit()) {
      TIterator *mcIter = _theModel->_masterSplitCat->typeIterator();
      while (RooCatType *mcState = (RooCatType*)mcIter->Next()){
	string key = string(spec->GetName())+":"+string(mcState->GetName());
	specsum->setVal(specsum->getVal()+yieldSet[key]);
      }
      delete mcIter;
    } else {
      // Not split..
      string key = spec->GetName();
      specsum->setVal(specsum->getVal()+yieldSet[key]);
    }    
  }

     cout << "In calc" << endl;
   for (map<string,int>::iterator iter = yieldSet.begin(); iter != yieldSet.end(); iter++){
    cout << "Key: " << iter->first << ", value: " << iter->second << endl;
  }

  /*
  double total = 0;
  TIterator *yIter = yieldSet.createIterator();
  while (RooRealVar* y = (RooRealVar*)yIter->Next()) {
    total += y->getVal();
  }
  cout << "Total sum: " << total << endl;
  */
}



void MLGenerator::generateComponent(TString species, const RooArgSet &whatVars, Int_t ngen, RooDataSet &theData, Bool_t verbose)
{
  MLSpecies *spec = _theFit->getSpecies(species);
  if (spec == NULL) return;
  generateComponent(*spec, whatVars, ngen, theData, verbose);
}


void MLGenerator::generateComponent(MLSpecies &spec, const RooArgSet &whatVars, Int_t ngen, RooDataSet &theData, Bool_t verbose)
{

  if (ngen == 0) return;
  
  // Get the PDF to generate with...
  RooAbsPdf *pdf = _theFit->getPdfComponent(*_theModel, spec.GetName());
  
  // Get the prototype dataset.  Default is no prototype.
  RooArgSet genVars(whatVars);
  RooDataSet *prototype = generatePrototype(spec, genVars, ngen, verbose);

  // Generate the data...
  RooDataSet *data = 0;
  if (prototype == 0){
    data = pdf->generate(genVars, ngen, verbose);
  } else {
    if (verbose) {
      cout << "Generating from prototype data set:" << endl;
      prototype->Print("V");
    }
    data = pdf->generate(genVars, *prototype, ngen, verbose);
  }
  delete prototype;
  
  // Generate the splitting categories if there are any
  RooDataSet *catData(0);
  if (_theModel->isSplit()) {
    RooArgSet catSet(_theModel->_masterSplitCat->inputCatList());
    catData = new RooDataSet("catData", "catData", catSet);
    for (Int_t i = 0; i < ngen; i++) catData->add(catSet);
    data->merge(catData);
    delete catData;
  }
  

  // Generate the MC truth if requested.
  RooDataSet *truthData = 0;
  if (_doTruth){
    truthData = new RooDataSet("truthData", "MC Truth Data", RooArgSet(*_truth));
    _truth->setLabel(spec.GetName());
    for (Int_t i = 0; i < ngen; i++) truthData->add(*_truth);
    data->merge(truthData);
    delete truthData;
  }

  Int_t nentries = theData.numEntries();
  theData.append(*data);
  if (nentries == 0) {
    // If this is the first time, theData might not have all the data columns (from the prototype dataset).
    // Call merge to add these columns if they're not there yet.
    theData.merge(data);
  }

  delete data;
  
}

void MLGenerator::mockComponent(TString filename, TString dataname, MLSpecies &spec, const RooArgSet &whatVars, Int_t ngen, RooDataSet &theData, Bool_t verbose)
{

  if (ngen == 0) return;

  TFile *mockfile = TFile::Open(filename);
  RooDataSet *mockDataSetFull = (RooDataSet*)mockfile->Get(dataname);

  if (!mockDataSetFull) {
    cout << "MLGenerator ERROR: dataset with name " << dataname << " doesn't exist for file " 
         << filename << " used to mock the species " << spec.GetName() 
         << " so not adding this species in the generated dataset. " << endl;
    return;
  }

  Int_t nentries = mockDataSetFull->numEntries();

  std::vector<Int_t> usedEvents;
  usedEvents.clear();
  
  for (int i=0; i<ngen; i++) {

    Int_t ntrials = 0;
    Int_t ievt = -1;
    while (ntrials <= nentries) {
      Double_t r = RooRandom::uniform();
      Int_t ievtmp = (Int_t) (r * nentries);
      if (find(usedEvents.begin(),usedEvents.end(),ievtmp) == usedEvents.end()) {
        ievt = ievtmp;
        usedEvents.push_back(ievt);
        break;
      }
      ntrials++;
    }

    if (ievt<0) {
      cout << "ERROR, mock dataset in file " << filename << " for species " << spec.GetName() 
           << " has too few events to generate even a single experiment. " << endl
           << "Not adding any event for this species to the generated dataset. " << endl;
      return;
    }

    char buf[100];
    sprintf(buf,"event==%d",ievt);
    RooDataSet *oneEvent = (RooDataSet*)mockDataSetFull->reduce(buf);
    if (oneEvent && oneEvent->numEntries()>0) theData.append(*oneEvent);
    else cout << "MLGenerator ERROR: event variable not found in the mock dataset?" << endl;
  }

}

RooDataSet* MLGenerator::generatePrototype(const MLSpecies &spec, RooArgSet &genVars, Int_t ngen, Bool_t verbose)
{
  return 0;
}

void MLGenerator::initializeRandomSC() {}

void MLGenerator::setActualNsig(int actualNsig){}

void MLGenerator::setMockSpecies(MLStrList specList, MLStrList specFileList, MLStrList specDatasetNameList) {

  _mockSpecies.Add(specList);
  _mockFiles.Add(specFileList);
  _mockDatasetNames.Add(specDatasetNameList);

}

