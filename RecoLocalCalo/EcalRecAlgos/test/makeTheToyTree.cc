void makeToyresultsTuple(TString dir = "./", TString file = "results.dat", 
			 TString varfile = "variables.root") {


  gSystem->Load("libRooFit");
  gSystem->Load("../../../MLFit/workdir/libMLFit.so");
  
  dir.Append("/");
  file.Prepend(dir);
  varfile.Prepend(dir);

  cout << "Reading data from:      " << file << endl;
  cout << "Reading variables from: " << varfile << endl;


  TFile f(varfile);
  RooArgSet *variables = (RooArgSet*)f.Get("variables");
  RooDataSet *fitResData = RooDataSet::read(file, *variables, "D");
  TreeFillerFromRooDataSet treeFiller(fitResData);
  treeFiller.addVar("N_sigIT_0");
  treeFiller.addVar("N_sigITgen");
  treeFiller.addVar("N_sigITerr_0");

  treeFiller.addVar("sigIT_Time_mean_0");
  treeFiller.addVar("sigIT_Time_meangen");
  treeFiller.addVar("sigIT_Time_meanerr_0");

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

  TFile *newfile = TFile::Open("toyresults2.root","recreate");
  TTree *ntp = treeFiller.getTree();
  ntp->Write();
  newfile->Close();
}

