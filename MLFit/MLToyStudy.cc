/**************************************
 *MLToyStudy
 *Nick Danielson
 ***************************************/

#include <TTree.h>
#include "MLFit/MLToyStudy.hh"
#include "MLFit/MLToyFit.hh"
#include "MLFit/MLFit.hh"
#include "MLFit/MLGenerator.hh"
#include "MLFit/MLStrList.hh"
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <TFile.h>
#include <TBranch.h>

#include <vector>

ClassImp(MLToyStudy);

MLToyStudy::MLToyStudy(const MLGenerator &theGenerator, const RooArgSet &dependents, const char* genOptions, const char* fitOptions, const RooDataSet *genProtoData, const RooArgSet& projDeps, const RooArgSet& extConstr) :
  _theGenerator((MLGenerator*)&theGenerator),
  _dependents(dependents),
  _genProtoData(genProtoData),
  _projDeps(projDeps),
  _fitOptions(fitOptions),
  _extConstr(extConstr)
{
  // Decode generator options
  TString genOpt(genOptions);
  genOpt.ToLower();
  _verboseGen = genOpt.Contains("v");
  _extendedGen = genOpt.Contains("e");
  
  
  // Get the set of parameters that were requested to be generated
  RooArgSet *genparams = theGenerator.getMasterPdf()->getParameters(&dependents);
  _genParams = (RooArgSet*)genparams->snapshot(kFALSE);
  delete genparams;
}


MLToyStudy::MLToyStudy(const MLGenerator &theGenerator, RooAbsPdf &theFit, const RooArgSet &dependents, const char* genOptions, const char* fitOptions, const RooDataSet *genProtoData, const RooArgSet& projDeps, const RooArgSet &extConstr) :
  _theGenerator((MLGenerator*)&theGenerator),
  _dependents(dependents),
  _genProtoData(genProtoData),
  _projDeps(projDeps),
  _fitOptions(fitOptions),
  _extConstr(extConstr)

{
  // Decode generator options
  TString genOpt(genOptions);
  genOpt.ToLower();
  _verboseGen = genOpt.Contains("v");
  _extendedGen = genOpt.Contains("e");
  

  // Get the set of parameters that were requested to be generated
  RooArgSet *genparams = theGenerator.getMasterPdf()->getParameters(&dependents);
  _genParams = (RooArgSet*)genparams->snapshot(kFALSE);
  delete genparams;

  // Add theFit to the list of fits to do
  addFit(theFit,RooArgSet(),RooArgSet());
  
}

MLToyStudy::~MLToyStudy()
{
}

void MLToyStudy::addFit(const RooAbsPdf &pdf, const RooArgSet &fixset, const RooArgSet &floatset, TString extraCut)
{
  //_theGenerator->getExtraParams().Print("V");
  MLToyFit *toyfit = new MLToyFit(pdf, _dependents, _theGenerator->getExtraParams());
  _fitList.Add(toyfit);
  toyfit->_fixSet.add(fixset);
  toyfit->_floatSet.add(floatset);
  toyfit->setDataSetCut(extraCut);
}

RooDataSet* MLToyStudy::generate(Int_t nEvt)
{
   *_theGenerator->getMasterPdf()->getParameters(_dependents) = *_genParams;
   _theGenerator->getMasterPdf()->getParameters(_dependents)->Print("V");
    
   cout << "Generating..." << endl;
   _theGenerator->setMockSpecies(_mockList,_mockFileList,_mockDatasetNameList);
   return _theGenerator->generate(_dependents,nEvt);
}

Bool_t MLToyStudy::generateAndFit(Int_t nSamples, Int_t nEvtPerSample, TString toy, TString signalFile)
{
  // Make sure you have a fit to perform!
  if (_fitList.GetSize() == 0) {
    cout << "MLToyStudy Error: You must first add a fit to perform with the addFit() command!" << endl;
    return kFALSE;
  }
  
  // Reset fitpardata;
  TIterator *tfIter = _fitList.MakeIterator();
  while(MLToyFit *toyfit = (MLToyFit*)tfIter->Next()) toyfit->_fitParData->reset();
  
  // Loop over samples...
  for (Int_t i = 1; i <= nSamples; i++) {
    cout << "MLToyStudy Generating experiment " << i << endl;
    RooDataSet *genSample(0);
    Int_t nEvt = nEvtPerSample;
    if (_extendedGen) nEvt = RooRandom::randomGenerator()->Poisson(nEvtPerSample);
  
    if (toy=="SCAN") {
       cout << "MLToyStudy::Randomly generating S and C" << endl;
      _theGenerator->initializeRandomSC();
      // Get the set of parameters that were requested to be generated once again
      RooArgSet *genparams = _theGenerator->getMasterPdf()->getParameters(_dependents);
      _genParams = (RooArgSet*)genparams->snapshot(kFALSE);
      delete genparams;
    }

    if(toy=="SCAN") {
      genSample = _theGenerator->generate(_dependents,nEvt);
    } else {
      cout << "MLToyStudy::Simple toy" << endl;
      genSample = generate(nEvt);
      genSample->Print("V");
      TFile* file = new TFile("gensample.root","RECREATE");
      genSample->Write();
      file->Close();
    }

    if(toy=="MCT") {
      TFile* file = new TFile(signalFile);
      RooDataSet *signalFromMc = (RooDataSet*) file->Get("dataset");      
      genSample->append(*signalFromMc);
      int nSigevents = signalFromMc->numEntries();
      _theGenerator->setActualNsig(nSigevents);
    }

    // And do each fit...
    tfIter->Reset();
    while(MLToyFit *toyfit = (MLToyFit*)tfIter->Next()) toyfit->initializeParameters();

    Int_t nfit = 0;
    tfIter->Reset();
    while (MLToyFit* toyfit = (MLToyFit*)tfIter->Next()) {
      cout << "Fitting sample " << i << " with " << toyfit->GetName() << " fit." << endl;
      if (toyfit->hasExtraDataSetCut()) {
	RooAbsData *dataToFit = genSample->reduce(toyfit->getDataSetCut());
	cout << "Cut: " << toyfit->getDataSetCut() << " leaves " << dataToFit->numEntries() << " of " << genSample->numEntries() << " entries." << endl;
	fitSample(dataToFit,*toyfit);
	delete dataToFit;
      } else {
	fitSample(genSample,*toyfit);
      }
      nfit++;
    } 
  }
  
  

  // All done.  Merge the various fit parameter data sets into one.
  mergeFitParData();
  delete tfIter;
  
  return kTRUE;
}



void MLToyStudy::fitSample(RooAbsData *data, MLToyFit &toyfit) 
{
  // Reset all fit parameters to their initial values
  //toyfit.initializeParameters();
  toyfit._fixSet.setAttribAll("Constant",kTRUE);
  toyfit._floatSet.setAttribAll("Constant",kFALSE);

  // Make sure when we're fitting we return a fit result...
  TString fitopt = TString(_fitOptions) + "r";
  
  // Fit!
  toyfit._pdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->Print("V");
  
  RooFitResult *fr(0);
  if (_projDeps.getSize()>0) {
    fr = toyfit._pdf->fitTo(*data,RooFit::ConditionalObservables(_projDeps),RooFit::FitOptions(fitopt.Data()),RooFit::ExternalConstraints(_extConstr));
  } else {
    fr = toyfit._pdf->fitTo(*data, RooFit::FitOptions(fitopt.Data()),RooFit::ExternalConstraints(_extConstr));
  }
      
  fr->Print("V");
  
  // Set values of extra variables...
  toyfit._nllVar->setVal(fr->minNll());
  toyfit._status->setVal(fr->status());
  toyfit._covQual->setVal(fr->covQual());
  cout << "Before i set themL " << endl;
  _theGenerator->getExtraParams().Print("V");
  toyfit._extraParams = _theGenerator->getExtraParams();
  toyfit._extraParams.Print("V");
  
  // Add an entry to the fit results data set...
  RooArgSet fitparams(*toyfit._fitParams);
  fitparams.add(RooArgSet(*toyfit._nllVar,*toyfit._status,*toyfit._covQual));
  fitparams.add(toyfit._extraParams);
  toyfit._fitParData->add(fitparams);
  
  delete fr;
}

void MLToyStudy::mergeFitParData() 
{

  // Merge the generated parameter data and the various fit parameter data into one data set:
  // _fitParData.

  // create the merged RooArgSet of variables
  RooArgSet totalset;
  Int_t nSamples=0;
  Int_t nfit = 0;
  TIterator *tfIter = _fitList.MakeIterator();
  while(MLToyFit *toyfit = (MLToyFit*)tfIter->Next()) {
    // take the number of toys from one of the fits (all have the same)
    nSamples = toyfit->_fitParData->numEntries();

    // First add error variables to all the fit parameters
    toyfit->addErrVars();

    TIterator *parIter = toyfit->_fitParData->get()->createIterator();
    while (RooRealVar *par = (RooRealVar*)parIter->Next()) {
      char cnfit[8];
      sprintf(cnfit,"%d",nfit);
      RooAbsReal *newPar = (RooAbsReal*)par->Clone(TString(par->GetName()) + "_" + TString(cnfit));
      totalset.add(*newPar);
    }
    nfit++;
  }  

  // create a empty RooDataset with those variables
  _fitParData = new RooDataSet("fitParData", "Fit Parameters Dataset", totalset);

  // fill the _fitParData
  for(int iexp=0; iexp<nSamples; iexp++) {
    // setting values
    nfit = 0;
    TIterator *tfIterEvent = _fitList.MakeIterator();
    while(MLToyFit *toyfit = (MLToyFit*)tfIterEvent->Next()) {
      TIterator *parIterEvent = toyfit->_fitParData->get(iexp)->createIterator();
      while (RooRealVar *par = (RooRealVar*)parIterEvent->Next()) {
        char cnfit[8];
        sprintf(cnfit,"%d",nfit);
        TString newName = TString(par->GetName()) + "_" + TString(cnfit);
        RooRealVar *newPar = (RooRealVar*)totalset.find(newName.Data());
        newPar->setVal(par->getVal());

      }
      nfit++;
    }
    _fitParData->add(totalset);
    
  } // loop over events

  // And add values of generated variables now...
  TIterator *parIter = _genParams->createIterator();
  while (RooRealVar *par = (RooRealVar*)parIter->Next()) {
    RooAbsReal *genPar = (RooAbsReal*)par->Clone("truth");
    TString name = TString(par->GetName()) + "gen";
    TString title = TString(par->GetName()) + " Generated";
    genPar->SetName(name);
    genPar->SetTitle(title);
    _fitParData->addColumn(*genPar);
    delete genPar;
  }
  delete parIter;

}

void MLToyStudy::addMockSpecies(const char* sname, const char* sfile, const char* sdatasetname) {

  _mockList.Add(sname);
  _mockFileList.Add(sfile);
  _mockDatasetNameList.Add(sdatasetname);

}
