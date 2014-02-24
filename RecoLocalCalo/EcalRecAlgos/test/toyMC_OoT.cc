// Set Fit Options
MLOptions GetDefaultOptions() {
  MLOptions opts;

  // Fit configuration
  opts.addBoolOption("useMass",           "Use Invariant Mass",     kTRUE);
  return opts;
}


void Generate(Int_t nexp = 1, UInt_t iseed = 65539, char* outfile= 0)
{

  // Various fit options...
  MLOptions opts = GetDefaultOptions();
  
  const char* e0 = "\033[44;37m";
  const char* en="\033[0m";
  
  TDatime *now = new TDatime();
  int today = now->GetDate();
  int clock = now->GetTime();
  int seed = today+clock+iseed;
  // Set the random number seed...
  RooRandom::randomGenerator()->SetSeed(seed);

  // define the structure of the dataset
  RooRealVar* time = new RooRealVar("time",  "time" , 0., 10., "samples");

  MLFit theFit;

  theFit.AddFlatFileColumn(time);

  // define a fit model
  theFit.addModel("myFit", "Pulse Shape fit");
  
  // define species
  theFit.addSpecies("myFit", "sigIT", "Signal In-Time Component");
  for(int i=1; i<5; ++i) {
    stringstream ootSpecName, ootSpecTitle;
    ootSpecName << "bkgOoT_minus" << i << "bx";
    ootSpecTitle << "Bkg Component OoT minus " << i << " bx";
    theFit.addSpecies("myFit", ootSpecName.str().c_str(), ootSpecTitle.str().c_str());
  }

  // mLL PDF
  theFit.addPdfWName("myFit", "sigIT" , "time",  "Pulse",  "sigIT_Time");
  for(int i=1; i<5; ++i) {
    stringstream ootSpecName, ootPdfName;
    ootSpecName << "bkgOoT_minus" << i << "bx";
    ootPdfName << "bkgOoT_minus" << i << "bx_Time";
    theFit.addPdfWName("myFit", ootSpecName.str().c_str() , "time",  "Pulse",    ootPdfName.str().c_str());    
  }


  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  theFit.initialize("toyMC_OoT.config");
  MLGenerator theGenerator(theFit, "myFit");
  
  Int_t ngen =
    theFit.getRealPar("N_sigIT")->getVal();
  for(int i=1; i<5; ++i) {
    stringstream yieldBkgName;
    yieldBkgName << "N_bkgOoT_minus" << i << "bx";
    ngen += theFit.getRealPar(yieldBkgName.str().c_str())->getVal();
  }

  // Generate...
  RooArgSet genVars(theFit.getObsList(MLStrList("time")));
  MLToyStudy theStudy(theGenerator, genVars, "E", "MTE", 0, theFit.getNoNormVars("myFit"));
  theStudy.addFit(*myPdf);

  theStudy.generateAndFit(nexp,ngen);
  
  theStudy._fitParData->write("results.dat");

  TFile varfile("variables.root","RECREATE");

  RooArgSet *variables = theStudy._fitParData->get();
  variables->setName("variables");
  variables->Write();
  varfile.Close();

  makeToyresultsTuple();

}

void makeToyresultsTuple(TString dir = "./", TString file = "results.dat", 
			 TString varfile = "variables.root") {


  dir.Append("/");
  file.Prepend(dir);
  varfile.Prepend(dir);

  cout << "Reading data from:      " << file << endl;
  cout << "Reading variables from: " << varfile << endl;


  TFile f(varfile);
  RooArgSet *variables = (RooArgSet*)f.Get("variables");
  RooDataSet *fitResData = RooDataSet::read(file, *variables);
  TreeFillerFromRooDataSet treeFiller(fitResData);
  treeFiller.addVar("N_sigIT_0");
  treeFiller.addVar("N_sigITgen");
  treeFiller.addVar("N_sigITerr_0");

  treeFiller.addVar("N_bkgOoT_minus1bx_0");
  treeFiller.addVar("N_bkgOoT_minus1bxgen");
  treeFiller.addVar("N_bkgOoT_minus1bxerr_0");

  treeFiller.addVar("N_bkgOoT_minus2bx_0");
  treeFiller.addVar("N_bkgOoT_minus2bxgen");
  treeFiller.addVar("N_bkgOoT_minus2bxerr_0");

  treeFiller.addVar("N_bkgOoT_minus3bx_0");
  treeFiller.addVar("N_bkgOoT_minus3bxgen");
  treeFiller.addVar("N_bkgOoT_minus3bxerr_0");

  treeFiller.addVar("N_bkgOoT_minus4bx_0");
  treeFiller.addVar("N_bkgOoT_minus4bxgen");
  treeFiller.addVar("N_bkgOoT_minus4bxerr_0");

  treeFiller.addVar("covQual_0");

  TFile *newfile = TFile::Open("toyresults.root","recreate");
  TTree *ntp = treeFiller.getTree();
  ntp->Write();
  newfile->Close();
}

