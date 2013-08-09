#include <iostream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"




int main(int argc, char* argv[]) {

  if( argc < 2 ) {
    std::cout << "USAGE: ./create_rhoWeights [dataset]" << std::endl;
    std::cout << "where [dataset] can be either \"ZJet\" or \"JetPD\" or \"MB\"" << std::endl;
    exit(11);
  }


  std::string dataset(argv[1]);
  if( dataset!="MB" && dataset!="JetPD" && dataset!="ZJet" ) {
    std::cout << "Only supporting \"JetPD\" and \"MB\" and \"ZJet\" at the moment. Exiting." << std::endl;
    exit(13);
  }

  std::string mcDataset;
  std::string dataDataset;
  if( dataset=="MB" ) {
    mcDataset = "flatQCD_P6_Dijets_12Jul"; 
    dataDataset = "data2012ABCD_MBPD_12Jul"; 
  } else if( dataset=="JetPD" ) {
    mcDataset = "flatQCD_P6_Dijets_12Jul"; 
    dataDataset = "data2012ABCD_JetPD_12Jul"; 
  } else if( dataset=="ZJet" ) {
    mcDataset = "Zjets_12Jul";
    dataDataset = "data2012ABCD_MuPD_12Jul"; 
  }

  std::string selectionType = (dataset=="ZJet") ? "ZJet" : "DiJet";

  std::string mcFileName = "sunilFlat_" + selectionType + "_" + mcDataset + ".root";
  std::string dataFileName = "sunilFlat_" + selectionType + "_" + dataDataset + ".root";

  TFile* mcFile = TFile::Open(mcFileName.c_str());
  TFile* dataFile = TFile::Open(dataFileName.c_str());

  std::cout << "-> Opened MC file " << mcFileName << std::endl;
  std::cout << "-> Opened DATA file " << dataFileName << std::endl;

  TTree* mcTree = (TTree*)mcFile->Get("tree_passedEvents");
  TTree* dataTree = (TTree*)dataFile->Get("tree_passedEvents");

  TH1D* h1_rho_data = new TH1D("rho_data", "", 90, 0., 45.);
  TH1D* h1_rho_mc = new TH1D("rho_mc", "", 90., 0., 45.);

  mcTree->Project( "rho_mc", "rhoPF");
  dataTree->Project( "rho_data", "rhoPF");

 
  std::string weightFileName = "rhoWeights_" + mcDataset + ".root";
  TFile* file_weights = TFile::Open(weightFileName.c_str(), "recreate");
  file_weights->cd();

  TH1D* h1_weights = new TH1D(*h1_rho_data);
  h1_weights->Divide(h1_rho_mc);
  h1_weights->SetName("rho_weights");

  h1_weights->Write();
  file_weights->Close(); 

  std::cout << "-> Created weights file " << weightFileName << std::endl;

  return 0;

}
