#include <cstdlib>
#include <cmath>
#include "TH1D.h"
#include "../interface/QGSyst.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "DrawBase.h"
#include "fitTools.h"


bool doubleMin = true;
bool sunilTree = false;


void drawSinglePlot( DrawBase* db, const std::string& discrim, TH1D* h1_data, TH1D* h1_qgl, TH1D* h1_qglSyst, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax );
void drawEffVsPt( DrawBase* db, const std::string& discrim, float thresh1, float thresh3, float etaMin, float etaMax, TH1D* h1_eff_gluon_thresh1, TH1D* h1_eff_gluon_syst_thresh1, TH1D* h1_eff_quark_thresh3, TH1D* h1_eff_quark_syst_thresh3);
void drawMC_beforeAfter( DrawBase* db, const std::string& discrim, TH1D* h1_qglJet_quark, TH1D* h1_qglJetSyst_quark, TH1D* h1_qglJet_gluon, TH1D* h1_qglJetSyst_gluon, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax );

int main( int argc, char* argv[] ) {



  std::string discrim = "qglJet";
  if( argc>1 ) {
    std::string discrim_str(argv[1]);
    if( discrim_str!="qglJet" && discrim_str!="qgMLPJet" ) {
      std::cout << "-> Unknown disciminator. Only \"qglJet\" and \"qgMLPJet\" allowed. Exiting." << std::endl;
      exit(191);
    }
    discrim = discrim_str;
  }



  std::string outfilename = "checkSyst_doubleMin";
  if( argc>2 ) {
    std::string doubleMin_str(argv[2]);
    if( doubleMin_str != "false" ) {
      doubleMin = true;
      std::cout << "-> Switching ON double min" << std::endl;
      outfilename = "checkSyst_doubleMin";
    } else {
      doubleMin = false;
      std::cout << "-> Switching OFF double min" << std::endl;
      outfilename = "checkSyst";
    }
  }
  outfilename = outfilename + "_" + discrim + ".root";


  TFile* file_puWeights = TFile::Open("puWeights.root");
  TH1D* h1_puweights = (TH1D*)file_puWeights->Get("puweights");


  //TFile* file = TFile::Open("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_2.root");
  TFile* file = TFile::Open("sunilFlat_ZJet_Zjets_12Jul.root");
  TTree* tree = (TTree*)file->Get("tree_passedEvents");

  //TFile* file_data = TFile::Open("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DoubleMu-Run2012AB.root");
  TChain* tree_data = new TChain("tree_passedEvents");
  tree_data->Add("sunilFlat_ZJet_data2012ABCD_MuPD_12Jul.root");
  //tree_data->Add("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_Double*.root");
  //std::cout << "tot: " << tree_data->GetEntries() << std::endl;
  //TTree* tree_data = (TTree*)file_data->Get("tree_passedEvents");


  int njet;
  float pt[20];
  float eta[20];
  float phi[20];
  int pdgId[20];
  float rho;
  float rhoMLP;
  int nCharged[20];
  int nNeutral[20];
  float ptD[20];
  float ptD_QC[20];
  float axis2_QC[20];
  float axis1_QC[20];
  float rmsCand_QC[20];
  float R[20];
  float pull_QC[20];
  int nCharged_QC[20];
  int nNeutral_ptCut[20];
  float qglJet[20];

  tree->SetBranchAddress("nJet", &njet);
  tree->SetBranchAddress("ptJet", pt);
  tree->SetBranchAddress("etaJet", eta);
  tree->SetBranchAddress("pdgIdJet", pdgId);
  tree->SetBranchAddress("rhoPF", &rho);
  tree->SetBranchAddress("nChargedJet", nCharged);
  tree->SetBranchAddress("nNeutralJet", nNeutral);
  tree->SetBranchAddress("ptDJet", ptD);
  tree->SetBranchAddress("ptD_QCJet", ptD_QC);
  tree->SetBranchAddress("axis1_QCJet", axis1_QC);
  tree->SetBranchAddress("axis2_QCJet", axis2_QC);
  tree->SetBranchAddress("rmsCand_QCJet", rmsCand_QC);
  tree->SetBranchAddress("RJet", R);
  tree->SetBranchAddress("pull_QCJet", pull_QC);
  tree->SetBranchAddress("nChg_QCJet", nCharged_QC);
  tree->SetBranchAddress("nNeutral_ptCutJet", nNeutral_ptCut);
  tree->SetBranchAddress(discrim.c_str(), qglJet);


  QGSyst qgsyst;
  if( doubleMin )
    qgsyst.ReadDatabaseDoubleMin("../data/SystZJetHbb_2013_07_15.txt");
    //qgsyst.ReadDatabaseDoubleMin("../data/SystDatabase_doubleMin.txt");
    //qgsyst.ReadDatabaseDoubleMin("../data/SystDoubleZJetHbb.txt");
  else
    qgsyst.ReadDatabase("../data/SystDatabase.txt");

  std::string nameForSyst = (discrim=="qglJet") ? "QGLHisto" : "QGLMLP";
  qgsyst.SetTagger(nameForSyst);



  TFile* outfile = TFile::Open(outfilename.c_str(), "RECREATE");
  outfile->cd(); 

  TH1D* h1_ptZ = new TH1D("ptZ", "", 100, 0., 100.);
  h1_ptZ->Sumw2();
  TH1D* h1_rho = new TH1D("rho", "", 40, 0., 40.);
  h1_rho->Sumw2();


  float xMin = 0.;

  int nBins = 50;
  int nBinsFwd = 25;

  TH1D* h1_qglJet_pt3050_eta02_rho040 = new TH1D("qglJet_pt3050_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_pt3050_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt3050_eta02_rho040 = new TH1D("qglJetSyst_pt3050_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_pt3050_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt3050_eta02_rho040 = new TH1D("qglJet_quark_pt3050_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_quark_pt3050_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt3050_eta02_rho040 = new TH1D("qglJetSyst_quark_pt3050_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_quark_pt3050_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt3050_eta02_rho040 = new TH1D("qglJet_gluon_pt3050_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_gluon_pt3050_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt3050_eta02_rho040 = new TH1D("qglJetSyst_gluon_pt3050_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt3050_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_pt5065_eta02_rho040 = new TH1D("qglJet_pt5065_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_pt5065_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt5065_eta02_rho040 = new TH1D("qglJetSyst_pt5065_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_pt5065_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt5065_eta02_rho040 = new TH1D("qglJet_quark_pt5065_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_quark_pt5065_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt5065_eta02_rho040 = new TH1D("qglJetSyst_quark_pt5065_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_quark_pt5065_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt5065_eta02_rho040 = new TH1D("qglJet_gluon_pt5065_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_gluon_pt5065_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt5065_eta02_rho040 = new TH1D("qglJetSyst_gluon_pt5065_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt5065_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_pt80120_eta02_rho040 = new TH1D("qglJet_pt80120_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_pt80120_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt80120_eta02_rho040 = new TH1D("qglJetSyst_pt80120_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_pt80120_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt80120_eta02_rho040 = new TH1D("qglJet_quark_pt80120_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_quark_pt80120_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt80120_eta02_rho040 = new TH1D("qglJetSyst_quark_pt80120_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_quark_pt80120_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt80120_eta02_rho040 = new TH1D("qglJet_gluon_pt80120_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_gluon_pt80120_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt80120_eta02_rho040 = new TH1D("qglJetSyst_gluon_pt80120_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt80120_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_pt120250_eta02_rho040 = new TH1D("qglJet_pt120250_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_pt120250_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt120250_eta02_rho040 = new TH1D("qglJetSyst_pt120250_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_pt120250_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt120250_eta02_rho040 = new TH1D("qglJet_quark_pt120250_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_quark_pt120250_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt120250_eta02_rho040 = new TH1D("qglJetSyst_quark_pt120250_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_quark_pt120250_eta02_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt120250_eta02_rho040 = new TH1D("qglJet_gluon_pt120250_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJet_gluon_pt120250_eta02_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt120250_eta02_rho040 = new TH1D("qglJetSyst_gluon_pt120250_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt120250_eta02_rho040->Sumw2();



  TH1D* h1_qglJet_pt3050_eta35_rho040 = new TH1D("qglJet_pt3050_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_pt3050_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt3050_eta35_rho040 = new TH1D("qglJetSyst_pt3050_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_pt3050_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt3050_eta35_rho040 = new TH1D("qglJet_quark_pt3050_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_quark_pt3050_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt3050_eta35_rho040 = new TH1D("qglJetSyst_quark_pt3050_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_quark_pt3050_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt3050_eta35_rho040 = new TH1D("qglJet_gluon_pt3050_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_gluon_pt3050_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt3050_eta35_rho040 = new TH1D("qglJetSyst_gluon_pt3050_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt3050_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_pt5065_eta35_rho040 = new TH1D("qglJet_pt5065_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_pt5065_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt5065_eta35_rho040 = new TH1D("qglJetSyst_pt5065_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_pt5065_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt5065_eta35_rho040 = new TH1D("qglJet_quark_pt5065_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_quark_pt5065_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt5065_eta35_rho040 = new TH1D("qglJetSyst_quark_pt5065_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_quark_pt5065_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt5065_eta35_rho040 = new TH1D("qglJet_gluon_pt5065_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_gluon_pt5065_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt5065_eta35_rho040 = new TH1D("qglJetSyst_gluon_pt5065_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt5065_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_pt80120_eta35_rho040 = new TH1D("qglJet_pt80120_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_pt80120_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt80120_eta35_rho040 = new TH1D("qglJetSyst_pt80120_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_pt80120_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt80120_eta35_rho040 = new TH1D("qglJet_quark_pt80120_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_quark_pt80120_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt80120_eta35_rho040 = new TH1D("qglJetSyst_quark_pt80120_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_quark_pt80120_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt80120_eta35_rho040 = new TH1D("qglJet_gluon_pt80120_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_gluon_pt80120_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt80120_eta35_rho040 = new TH1D("qglJetSyst_gluon_pt80120_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt80120_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_pt120250_eta35_rho040 = new TH1D("qglJet_pt120250_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_pt120250_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_pt120250_eta35_rho040 = new TH1D("qglJetSyst_pt120250_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_pt120250_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_quark_pt120250_eta35_rho040 = new TH1D("qglJet_quark_pt120250_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_quark_pt120250_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_quark_pt120250_eta35_rho040 = new TH1D("qglJetSyst_quark_pt120250_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_quark_pt120250_eta35_rho040->Sumw2();

  TH1D* h1_qglJet_gluon_pt120250_eta35_rho040 = new TH1D("qglJet_gluon_pt120250_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJet_gluon_pt120250_eta35_rho040->Sumw2();
  TH1D* h1_qglJetSyst_gluon_pt120250_eta35_rho040 = new TH1D("qglJetSyst_gluon_pt120250_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_qglJetSyst_gluon_pt120250_eta35_rho040->Sumw2();





  //syst eff histos:
  float thresh1 = (discrim=="qgMLPJet") ? 0.7 : 0.3;
  float thresh2 = (discrim=="qgMLPJet") ? 0.5 : 0.5;
  float thresh3 = (discrim=="qgMLPJet") ? 0.3 : 0.8;

  Double_t ptBins[21];
  fitTools::getBins_int( 21, ptBins, 20., 2000.);

  TH1D* h1_effDenom_centr_quark_thresh1 = new TH1D("effDenom_centr_quark_thresh1", "", 20, ptBins);
  TH1D* h1_effDenom_centr_quark_thresh2 = new TH1D("effDenom_centr_quark_thresh2", "", 20, ptBins);
  TH1D* h1_effDenom_centr_quark_thresh3 = new TH1D("effDenom_centr_quark_thresh3", "", 20, ptBins);

  TH1D* h1_effNum_centr_quark_thresh1   = new TH1D("effNum_centr_quark_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_quark_thresh2   = new TH1D("effNum_centr_quark_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_quark_thresh3   = new TH1D("effNum_centr_quark_thresh3", "", 20, ptBins)  ;

  TH1D* h1_effNum_centr_quark_syst_thresh1   = new TH1D("effNum_centr_quark_syst_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_quark_syst_thresh2   = new TH1D("effNum_centr_quark_syst_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_quark_syst_thresh3   = new TH1D("effNum_centr_quark_syst_thresh3", "", 20, ptBins)  ;

  TH1D* h1_effDenom_centr_gluon_thresh1 = new TH1D("effDenom_centr_gluon_thresh1", "", 20, ptBins);
  TH1D* h1_effDenom_centr_gluon_thresh2 = new TH1D("effDenom_centr_gluon_thresh2", "", 20, ptBins);
  TH1D* h1_effDenom_centr_gluon_thresh3 = new TH1D("effDenom_centr_gluon_thresh3", "", 20, ptBins);

  TH1D* h1_effNum_centr_gluon_thresh1   = new TH1D("effNum_centr_gluon_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_gluon_thresh2   = new TH1D("effNum_centr_gluon_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_gluon_thresh3   = new TH1D("effNum_centr_gluon_thresh3", "", 20, ptBins)  ;

  TH1D* h1_effNum_centr_gluon_syst_thresh1   = new TH1D("effNum_centr_gluon_syst_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_gluon_syst_thresh2   = new TH1D("effNum_centr_gluon_syst_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_centr_gluon_syst_thresh3   = new TH1D("effNum_centr_gluon_syst_thresh3", "", 20, ptBins)  ;

                                                             

  TH1D* h1_effDenom_fwd_quark_thresh1 = new TH1D("effDenom_fwd_quark_thresh1", "", 20, ptBins);
  TH1D* h1_effDenom_fwd_quark_thresh2 = new TH1D("effDenom_fwd_quark_thresh2", "", 20, ptBins);
  TH1D* h1_effDenom_fwd_quark_thresh3 = new TH1D("effDenom_fwd_quark_thresh3", "", 20, ptBins);

  TH1D* h1_effNum_fwd_quark_thresh1   = new TH1D("effNum_fwd_quark_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_quark_thresh2   = new TH1D("effNum_fwd_quark_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_quark_thresh3   = new TH1D("effNum_fwd_quark_thresh3", "", 20, ptBins)  ;

  TH1D* h1_effNum_fwd_quark_syst_thresh1   = new TH1D("effNum_fwd_quark_syst_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_quark_syst_thresh2   = new TH1D("effNum_fwd_quark_syst_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_quark_syst_thresh3   = new TH1D("effNum_fwd_quark_syst_thresh3", "", 20, ptBins)  ;

  TH1D* h1_effDenom_fwd_gluon_thresh1 = new TH1D("effDenom_fwd_gluon_thresh1", "", 20, ptBins);
  TH1D* h1_effDenom_fwd_gluon_thresh2 = new TH1D("effDenom_fwd_gluon_thresh2", "", 20, ptBins);
  TH1D* h1_effDenom_fwd_gluon_thresh3 = new TH1D("effDenom_fwd_gluon_thresh3", "", 20, ptBins);

  TH1D* h1_effNum_fwd_gluon_thresh1   = new TH1D("effNum_fwd_gluon_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_gluon_thresh2   = new TH1D("effNum_fwd_gluon_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_gluon_thresh3   = new TH1D("effNum_fwd_gluon_thresh3", "", 20, ptBins)  ;

  TH1D* h1_effNum_fwd_gluon_syst_thresh1   = new TH1D("effNum_fwd_gluon_syst_thresh1", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_gluon_syst_thresh2   = new TH1D("effNum_fwd_gluon_syst_thresh2", "", 20, ptBins)  ;
  TH1D* h1_effNum_fwd_gluon_syst_thresh3   = new TH1D("effNum_fwd_gluon_syst_thresh3", "", 20, ptBins)  ;
                                                             


  int nentries = tree->GetEntries();

  for( unsigned iEntry = 0; iEntry < nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( njet<1 ) continue;

    if( rho>40. ) continue;
    //if( mZ<70. || mZ>110. ) continue;
    //if( deltaPhi_jet < 3.1415-0.5 ) continue;
    //if( ptD_QCJet <=0. ) continue;
    //if( fabs(etaJet)<2.5 && betaStarJet > 0.2 * log( (float)nvertex - 0.67) ) continue;
    //if( btagged ) continue;

    float puweight = h1_puweights->GetBinContent((int)rho);

    //h1_ptZ->Fill( ptZ, puweight );
    h1_rho->Fill( rho, puweight );

    bool smearJet = true;

    std::string type = "all";
    if( doubleMin ) {
      if( fabs(pdgId[0])>0 && fabs(pdgId[0])<6 ) type = "quark";
      else if( pdgId[0]==21 ) type = "gluon";
      else { // both 0 and -999
        //smearJet = false;
        type = "gluon";
        //std::cout << "Unknown jet PDG ID (" << pdgIdPartJet << "). Will use gluon instead." << std::endl;
      }
    }


    float qglJetSyst = (smearJet) ? qgsyst.Smear(pt[0], eta[0], rho, qglJet[0], type) : qglJet[0];

    bool discrimOK_thresh1 = (discrim=="qglJet" && qglJet[0]>thresh1) || (discrim=="qgMLPJet" && qglJet[0]<thresh1);
    bool discrimSystOK_thresh1 = (discrim=="qglJet" && qglJetSyst>thresh1) || (discrim=="qgMLPJet" && qglJetSyst<thresh1);
    bool discrimOK_thresh2 = (discrim=="qglJet" && qglJet[0]>thresh2) || (discrim=="qgMLPJet" && qglJet[0]<thresh2);
    bool discrimSystOK_thresh2 = (discrim=="qglJet" && qglJetSyst>thresh2) || (discrim=="qgMLPJet" && qglJetSyst<thresh2);
    bool discrimOK_thresh3 = (discrim=="qglJet" && qglJet[0]>thresh3) || (discrim=="qgMLPJet" && qglJet[0]<thresh3);
    bool discrimSystOK_thresh3 = (discrim=="qglJet" && qglJetSyst>thresh3) || (discrim=="qgMLPJet" && qglJetSyst<thresh3);


    if( fabs(eta[0])<2. ) {


      if( type=="quark" ) {

        h1_effDenom_centr_quark_thresh1->Fill( pt[0] );
        h1_effDenom_centr_quark_thresh2->Fill( pt[0] );
        h1_effDenom_centr_quark_thresh3->Fill( pt[0] );

        if( discrimOK_thresh1 ) h1_effNum_centr_quark_thresh1->Fill( pt[0] );
        if( discrimOK_thresh2 ) h1_effNum_centr_quark_thresh2->Fill( pt[0] );
        if( discrimOK_thresh3 ) h1_effNum_centr_quark_thresh3->Fill( pt[0] );

        if( discrimSystOK_thresh1 || (discrimOK_thresh1 && pt[0]<30.) ) h1_effNum_centr_quark_syst_thresh1->Fill( pt[0] );
        if( discrimSystOK_thresh2 || (discrimOK_thresh2 && pt[0]<30.) ) h1_effNum_centr_quark_syst_thresh2->Fill( pt[0] );
        if( discrimSystOK_thresh3 || (discrimOK_thresh3 && pt[0]<30.) ) h1_effNum_centr_quark_syst_thresh3->Fill( pt[0] );

      } else if( type=="gluon" ) {

        h1_effDenom_centr_gluon_thresh1->Fill( pt[0] );
        h1_effDenom_centr_gluon_thresh2->Fill( pt[0] );
        h1_effDenom_centr_gluon_thresh3->Fill( pt[0] );

        if( discrimOK_thresh1 ) h1_effNum_centr_gluon_thresh1->Fill( pt[0] );
        if( discrimOK_thresh2 ) h1_effNum_centr_gluon_thresh2->Fill( pt[0] );
        if( discrimOK_thresh3 ) h1_effNum_centr_gluon_thresh3->Fill( pt[0] );

        if( discrimSystOK_thresh1 || (discrimOK_thresh1 && pt[0]<30.) ) h1_effNum_centr_gluon_syst_thresh1->Fill( pt[0] );
        if( discrimSystOK_thresh2 || (discrimOK_thresh2 && pt[0]<30.) ) h1_effNum_centr_gluon_syst_thresh2->Fill( pt[0] );
        if( discrimSystOK_thresh3 || (discrimOK_thresh3 && pt[0]<30.) ) h1_effNum_centr_gluon_syst_thresh3->Fill( pt[0] );

      }


      if( pt[0] > 30. && pt[0] < 50. ) {

        h1_qglJet_pt3050_eta02_rho040->Fill( qglJet[0], puweight );
        h1_qglJetSyst_pt3050_eta02_rho040->Fill( qglJetSyst, puweight );
        if( type=="quark" ) {
          h1_qglJet_quark_pt3050_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_quark_pt3050_eta02_rho040->Fill( qglJetSyst, puweight );
        } else if( type=="gluon" ) {
          h1_qglJet_gluon_pt3050_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_gluon_pt3050_eta02_rho040->Fill( qglJetSyst, puweight );
        }

      } else if( pt[0] > 50. && pt[0] < 65. ) {

        h1_qglJet_pt5065_eta02_rho040->Fill( qglJet[0], puweight );
        h1_qglJetSyst_pt5065_eta02_rho040->Fill( qglJetSyst, puweight );
        if( type=="quark" ) {
          h1_qglJet_quark_pt5065_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_quark_pt5065_eta02_rho040->Fill( qglJetSyst, puweight );
        } else if( type=="gluon" ) {
          h1_qglJet_gluon_pt5065_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_gluon_pt5065_eta02_rho040->Fill( qglJetSyst, puweight );
        }

      } else if( pt[0] > 80. && pt[0] < 120. ) {

        h1_qglJet_pt80120_eta02_rho040->Fill( qglJet[0], puweight );
        h1_qglJetSyst_pt80120_eta02_rho040->Fill( qglJetSyst, puweight );
        if( type=="quark" ) {
          h1_qglJet_quark_pt80120_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_quark_pt80120_eta02_rho040->Fill( qglJetSyst, puweight );
        } else if( type=="gluon" ) {
          h1_qglJet_gluon_pt80120_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_gluon_pt80120_eta02_rho040->Fill( qglJetSyst, puweight );
        }

      } else if( pt[0] > 120. && pt[0] < 250. ) {

        h1_qglJet_pt120250_eta02_rho040->Fill( qglJet[0], puweight );
        h1_qglJetSyst_pt120250_eta02_rho040->Fill( qglJetSyst, puweight );
        if( type=="quark" ) {
          h1_qglJet_quark_pt120250_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_quark_pt120250_eta02_rho040->Fill( qglJetSyst, puweight );
        } else if( type=="gluon" ) {
          h1_qglJet_gluon_pt120250_eta02_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_gluon_pt120250_eta02_rho040->Fill( qglJetSyst, puweight );
        }

      } 

    } else if( fabs(eta[0])>3. && fabs(eta[0])<5. ) {

      if( type=="quark" ) {

        h1_effDenom_fwd_quark_thresh1->Fill( pt[0] );
        h1_effDenom_fwd_quark_thresh2->Fill( pt[0] );
        h1_effDenom_fwd_quark_thresh3->Fill( pt[0] );

        if( discrimOK_thresh1 ) h1_effNum_fwd_quark_thresh1->Fill( pt[0] );
        if( discrimOK_thresh2 ) h1_effNum_fwd_quark_thresh2->Fill( pt[0] );
        if( discrimOK_thresh3 ) h1_effNum_fwd_quark_thresh3->Fill( pt[0] );

        if( discrimSystOK_thresh1 || (discrimOK_thresh1 && pt[0]<30.) ) h1_effNum_fwd_quark_syst_thresh1->Fill( pt[0] );
        if( discrimSystOK_thresh2 || (discrimOK_thresh2 && pt[0]<30.) ) h1_effNum_fwd_quark_syst_thresh2->Fill( pt[0] );
        if( discrimSystOK_thresh3 || (discrimOK_thresh3 && pt[0]<30.) ) h1_effNum_fwd_quark_syst_thresh3->Fill( pt[0] );

      } else if( type=="gluon" ) {

        h1_effDenom_fwd_gluon_thresh1->Fill( pt[0] );
        h1_effDenom_fwd_gluon_thresh2->Fill( pt[0] );
        h1_effDenom_fwd_gluon_thresh3->Fill( pt[0] );

        if( discrimOK_thresh1 ) h1_effNum_fwd_gluon_thresh1->Fill( pt[0] );
        if( discrimOK_thresh2 ) h1_effNum_fwd_gluon_thresh2->Fill( pt[0] );
        if( discrimOK_thresh3 ) h1_effNum_fwd_gluon_thresh3->Fill( pt[0] );

        if( discrimSystOK_thresh1 || (discrimOK_thresh1 && pt[0]<30.) ) h1_effNum_fwd_gluon_syst_thresh1->Fill( pt[0] );
        if( discrimSystOK_thresh2 || (discrimOK_thresh2 && pt[0]<30.) ) h1_effNum_fwd_gluon_syst_thresh2->Fill( pt[0] );
        if( discrimSystOK_thresh3 || (discrimOK_thresh3 && pt[0]<30.) ) h1_effNum_fwd_gluon_syst_thresh3->Fill( pt[0] );

      }

      if( pt[0] > 30. && pt[0] < 50. ) {

        h1_qglJet_pt3050_eta35_rho040->Fill( qglJet[0], puweight );
        h1_qglJetSyst_pt3050_eta35_rho040->Fill( qglJetSyst, puweight );
        if( type=="quark" ) {
          h1_qglJet_quark_pt3050_eta35_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_quark_pt3050_eta35_rho040->Fill( qglJetSyst, puweight );
        } else if( type=="gluon" ) {
          h1_qglJet_gluon_pt3050_eta35_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_gluon_pt3050_eta35_rho040->Fill( qglJetSyst, puweight );
        }

      } else if( pt[0] > 50. && pt[0] < 65. ) {

        h1_qglJet_pt5065_eta35_rho040->Fill( qglJet[0], puweight );
        h1_qglJetSyst_pt5065_eta35_rho040->Fill( qglJetSyst, puweight );
        if( type=="quark" ) {
          h1_qglJet_quark_pt5065_eta35_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_quark_pt5065_eta35_rho040->Fill( qglJetSyst, puweight );
        } else if( type=="gluon" ) {
          h1_qglJet_gluon_pt5065_eta35_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_gluon_pt5065_eta35_rho040->Fill( qglJetSyst, puweight );
        }

      } else if( pt[0] > 80. && pt[0] < 120. ) {

        h1_qglJet_pt80120_eta35_rho040->Fill( qglJet[0], puweight );
        h1_qglJetSyst_pt80120_eta35_rho040->Fill( qglJetSyst, puweight );
        if( type=="quark" ) {
          h1_qglJet_quark_pt80120_eta35_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_quark_pt80120_eta35_rho040->Fill( qglJetSyst, puweight );
        } else if( type=="gluon" ) {
          h1_qglJet_gluon_pt80120_eta35_rho040->Fill( qglJet[0], puweight );
          h1_qglJetSyst_gluon_pt80120_eta35_rho040->Fill( qglJetSyst, puweight );
        }

      } 

    } // eta

  }


  // now switch to data:
  int event;

  tree_data->SetBranchAddress("event", &event );
  tree_data->SetBranchAddress("nJet", &njet);
  tree_data->SetBranchAddress("ptJet", pt);
  tree_data->SetBranchAddress("etaJet", eta);
  tree_data->SetBranchAddress("rhoPF", &rho);
  tree_data->SetBranchAddress("ptD_QCJet", ptD_QC);
  tree_data->SetBranchAddress(discrim.c_str(), qglJet);


  TH1D* h1_data_ptZ = new TH1D("data_ptZ", "", 100, 0., 100.);
  h1_data_ptZ->Sumw2();
  TH1D* h1_data_rho = new TH1D("data_rho", "", 40, 0., 40.);
  h1_data_rho->Sumw2();

  TH1D* h1_data_qglJet_pt3050_eta02_rho040 = new TH1D("data_qglJet_pt3050_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_data_qglJet_pt3050_eta02_rho040->Sumw2();

  TH1D* h1_data_qglJet_pt5065_eta02_rho040 = new TH1D("data_qglJet_pt5065_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_data_qglJet_pt5065_eta02_rho040->Sumw2();

  TH1D* h1_data_qglJet_pt80120_eta02_rho040 = new TH1D("data_qglJet_pt80120_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_data_qglJet_pt80120_eta02_rho040->Sumw2();

  TH1D* h1_data_qglJet_pt120250_eta02_rho040 = new TH1D("data_qglJet_pt120250_eta02_rho040", "", nBins, xMin, 1.0001);
  h1_data_qglJet_pt120250_eta02_rho040->Sumw2();


  TH1D* h1_data_qglJet_pt3050_eta35_rho040 = new TH1D("data_qglJet_pt3050_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_data_qglJet_pt3050_eta35_rho040->Sumw2();

  TH1D* h1_data_qglJet_pt5065_eta35_rho040 = new TH1D("data_qglJet_pt5065_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_data_qglJet_pt5065_eta35_rho040->Sumw2();

  TH1D* h1_data_qglJet_pt80120_eta35_rho040 = new TH1D("data_qglJet_pt80120_eta35_rho040", "", nBinsFwd, xMin, 1.0001);
  h1_data_qglJet_pt80120_eta35_rho040->Sumw2();


  int nentries_data = tree_data->GetEntries();

  for( unsigned iEntry = 0; iEntry < nentries_data; ++iEntry ) {

    tree_data->GetEntry(iEntry);

    if( njet<1 ) continue;
    if( rho>40. ) continue;

    //h1_data_ptZ->Fill( ptZ );
    h1_data_rho->Fill( rho );

    if( fabs(eta[0])<2. ) {

      if( pt[0] > 30. && pt[0] < 50. ) {

        h1_data_qglJet_pt3050_eta02_rho040->Fill( qglJet[0] );

      } else if( pt[0] > 50. && pt[0] < 65. ) {

//if( ptJet<65. )
//std::cout << "eventNumber" << event << std::endl;
        h1_data_qglJet_pt5065_eta02_rho040->Fill( qglJet[0] );

      } else if( pt[0] > 80. && pt[0] < 120. ) {

        h1_data_qglJet_pt80120_eta02_rho040->Fill( qglJet[0] );

      } else if( pt[0] > 120. && pt[0] < 250. ) {

        h1_data_qglJet_pt120250_eta02_rho040->Fill( qglJet[0] );

      } 

    } else if( fabs(eta[0])>3. && fabs(eta[0])<5. ) {

      if( pt[0] > 30. && pt[0] < 50. ) {

        h1_data_qglJet_pt3050_eta35_rho040->Fill( qglJet[0] );

      } else if( pt[0] > 50. && pt[0] < 65. ) {

        h1_data_qglJet_pt5065_eta35_rho040->Fill( qglJet[0] );

      } else if( pt[0] > 80. && pt[0] < 120. ) {

        h1_data_qglJet_pt80120_eta35_rho040->Fill( qglJet[0] );

      } 

    } // eta

  }


  DrawBase* db = new DrawBase("checkSystSyst");
  db->add_dataFile(file, "data"); //ok thats the MC file but who cares

  std::string outputdir = "checkSystPlots";
  db->set_outputdir(outputdir);
  std::string mkdir_command = "mkdir -p " + outputdir;
  system(mkdir_command.c_str());
  
  drawSinglePlot( db, discrim, h1_data_qglJet_pt3050_eta02_rho040, h1_qglJet_pt3050_eta02_rho040, h1_qglJetSyst_pt3050_eta02_rho040, 30., 50., 0., 2., 0., 40.);
  drawSinglePlot( db, discrim, h1_data_qglJet_pt5065_eta02_rho040, h1_qglJet_pt5065_eta02_rho040, h1_qglJetSyst_pt5065_eta02_rho040, 50., 65., 0., 2., 0., 40.);
  drawSinglePlot( db, discrim, h1_data_qglJet_pt80120_eta02_rho040,h1_qglJet_pt80120_eta02_rho040, h1_qglJetSyst_pt80120_eta02_rho040, 80., 120., 0., 2., 0., 40.);
  drawSinglePlot( db, discrim, h1_data_qglJet_pt120250_eta02_rho040,h1_qglJet_pt120250_eta02_rho040, h1_qglJetSyst_pt120250_eta02_rho040, 120., 250., 0., 2., 0., 40.);

  drawSinglePlot( db, discrim, h1_data_qglJet_pt3050_eta35_rho040, h1_qglJet_pt3050_eta35_rho040, h1_qglJetSyst_pt3050_eta35_rho040, 30., 50., 3., 5., 0., 40.);
  drawSinglePlot( db, discrim, h1_data_qglJet_pt5065_eta35_rho040, h1_qglJet_pt5065_eta35_rho040, h1_qglJetSyst_pt5065_eta35_rho040, 50., 65., 3., 5., 0., 40.);
  drawSinglePlot( db, discrim, h1_data_qglJet_pt80120_eta35_rho040,h1_qglJet_pt80120_eta35_rho040, h1_qglJetSyst_pt80120_eta35_rho040, 80., 120., 3., 5., 0., 40.);

  drawMC_beforeAfter( db, discrim, h1_qglJet_quark_pt3050_eta02_rho040, h1_qglJetSyst_quark_pt3050_eta02_rho040, h1_qglJet_gluon_pt3050_eta02_rho040, h1_qglJetSyst_gluon_pt3050_eta02_rho040, 30., 50., 0., 2., 0., 40.);
  drawMC_beforeAfter( db, discrim, h1_qglJet_quark_pt5065_eta02_rho040, h1_qglJetSyst_quark_pt5065_eta02_rho040, h1_qglJet_gluon_pt5065_eta02_rho040, h1_qglJetSyst_gluon_pt5065_eta02_rho040, 50., 65., 0., 2., 0., 40.);
  drawMC_beforeAfter( db, discrim, h1_qglJet_quark_pt80120_eta02_rho040, h1_qglJetSyst_quark_pt80120_eta02_rho040, h1_qglJet_gluon_pt80120_eta02_rho040, h1_qglJetSyst_gluon_pt80120_eta02_rho040, 80., 120., 0., 2., 0., 40.);
  drawMC_beforeAfter( db, discrim, h1_qglJet_quark_pt120250_eta02_rho040, h1_qglJetSyst_quark_pt120250_eta02_rho040, h1_qglJet_gluon_pt120250_eta02_rho040, h1_qglJetSyst_gluon_pt120250_eta02_rho040, 120., 250., 0., 2., 0., 40.);

  drawMC_beforeAfter( db, discrim, h1_qglJet_quark_pt3050_eta35_rho040, h1_qglJetSyst_quark_pt3050_eta35_rho040, h1_qglJet_gluon_pt3050_eta35_rho040, h1_qglJetSyst_gluon_pt3050_eta35_rho040, 30., 50., 3., 5., 0., 40.);
  drawMC_beforeAfter( db, discrim, h1_qglJet_quark_pt5065_eta35_rho040, h1_qglJetSyst_quark_pt5065_eta35_rho040, h1_qglJet_gluon_pt5065_eta35_rho040, h1_qglJetSyst_gluon_pt5065_eta35_rho040, 50., 65., 3., 5., 0., 40.);
  drawMC_beforeAfter( db, discrim, h1_qglJet_quark_pt80120_eta35_rho040, h1_qglJetSyst_quark_pt80120_eta35_rho040, h1_qglJet_gluon_pt80120_eta35_rho040, h1_qglJetSyst_gluon_pt80120_eta35_rho040, 80., 120., 3., 5., 0., 40.);

  // and now the syst smearing efficiencies:
  TH1D* h1_eff_centr_quark_thresh1 = new TH1D(*h1_effNum_centr_quark_thresh1);
  TH1D* h1_eff_centr_quark_thresh2 = new TH1D(*h1_effNum_centr_quark_thresh2);
  TH1D* h1_eff_centr_quark_thresh3 = new TH1D(*h1_effNum_centr_quark_thresh3);

  TH1D* h1_eff_centr_quark_syst_thresh1 = new TH1D(*h1_effNum_centr_quark_syst_thresh1);
  TH1D* h1_eff_centr_quark_syst_thresh2 = new TH1D(*h1_effNum_centr_quark_syst_thresh2);
  TH1D* h1_eff_centr_quark_syst_thresh3 = new TH1D(*h1_effNum_centr_quark_syst_thresh3);

  TH1D* h1_eff_centr_gluon_thresh1 = new TH1D(*h1_effNum_centr_gluon_thresh1);
  TH1D* h1_eff_centr_gluon_thresh2 = new TH1D(*h1_effNum_centr_gluon_thresh2);
  TH1D* h1_eff_centr_gluon_thresh3 = new TH1D(*h1_effNum_centr_gluon_thresh3);

  TH1D* h1_eff_centr_gluon_syst_thresh1 = new TH1D(*h1_effNum_centr_gluon_syst_thresh1);
  TH1D* h1_eff_centr_gluon_syst_thresh2 = new TH1D(*h1_effNum_centr_gluon_syst_thresh2);
  TH1D* h1_eff_centr_gluon_syst_thresh3 = new TH1D(*h1_effNum_centr_gluon_syst_thresh3);

  TH1D* h1_eff_fwd_quark_thresh1 = new TH1D(*h1_effNum_fwd_quark_thresh1);
  TH1D* h1_eff_fwd_quark_thresh2 = new TH1D(*h1_effNum_fwd_quark_thresh2);
  TH1D* h1_eff_fwd_quark_thresh3 = new TH1D(*h1_effNum_fwd_quark_thresh3);

  TH1D* h1_eff_fwd_quark_syst_thresh1 = new TH1D(*h1_effNum_fwd_quark_syst_thresh1);
  TH1D* h1_eff_fwd_quark_syst_thresh2 = new TH1D(*h1_effNum_fwd_quark_syst_thresh2);
  TH1D* h1_eff_fwd_quark_syst_thresh3 = new TH1D(*h1_effNum_fwd_quark_syst_thresh3);

  TH1D* h1_eff_fwd_gluon_thresh1 = new TH1D(*h1_effNum_fwd_gluon_thresh1);
  TH1D* h1_eff_fwd_gluon_thresh2 = new TH1D(*h1_effNum_fwd_gluon_thresh2);
  TH1D* h1_eff_fwd_gluon_thresh3 = new TH1D(*h1_effNum_fwd_gluon_thresh3);

  TH1D* h1_eff_fwd_gluon_syst_thresh1 = new TH1D(*h1_effNum_fwd_gluon_syst_thresh1);
  TH1D* h1_eff_fwd_gluon_syst_thresh2 = new TH1D(*h1_effNum_fwd_gluon_syst_thresh2);
  TH1D* h1_eff_fwd_gluon_syst_thresh3 = new TH1D(*h1_effNum_fwd_gluon_syst_thresh3);


  h1_eff_centr_quark_thresh1->Divide(h1_effDenom_centr_quark_thresh1);
  h1_eff_centr_quark_thresh2->Divide(h1_effDenom_centr_quark_thresh2);
  h1_eff_centr_quark_thresh3->Divide(h1_effDenom_centr_quark_thresh3);

  h1_eff_centr_quark_syst_thresh1->Divide(h1_effDenom_centr_quark_thresh1);
  h1_eff_centr_quark_syst_thresh2->Divide(h1_effDenom_centr_quark_thresh2);
  h1_eff_centr_quark_syst_thresh3->Divide(h1_effDenom_centr_quark_thresh3);

  h1_eff_centr_gluon_thresh1->Divide(h1_effDenom_centr_gluon_thresh1);
  h1_eff_centr_gluon_thresh2->Divide(h1_effDenom_centr_gluon_thresh2);
  h1_eff_centr_gluon_thresh3->Divide(h1_effDenom_centr_gluon_thresh3);

  h1_eff_centr_gluon_syst_thresh1->Divide(h1_effDenom_centr_gluon_thresh1);
  h1_eff_centr_gluon_syst_thresh2->Divide(h1_effDenom_centr_gluon_thresh2);
  h1_eff_centr_gluon_syst_thresh3->Divide(h1_effDenom_centr_gluon_thresh3);

  h1_eff_fwd_quark_thresh1->Divide(h1_effDenom_fwd_quark_thresh1);
  h1_eff_fwd_quark_thresh2->Divide(h1_effDenom_fwd_quark_thresh2);
  h1_eff_fwd_quark_thresh3->Divide(h1_effDenom_fwd_quark_thresh3);

  h1_eff_fwd_quark_syst_thresh1->Divide(h1_effDenom_fwd_quark_thresh1);
  h1_eff_fwd_quark_syst_thresh2->Divide(h1_effDenom_fwd_quark_thresh2);
  h1_eff_fwd_quark_syst_thresh3->Divide(h1_effDenom_fwd_quark_thresh3);

  h1_eff_fwd_gluon_thresh1->Divide(h1_effDenom_fwd_gluon_thresh1);
  h1_eff_fwd_gluon_thresh2->Divide(h1_effDenom_fwd_gluon_thresh2);
  h1_eff_fwd_gluon_thresh3->Divide(h1_effDenom_fwd_gluon_thresh3);

  h1_eff_fwd_gluon_syst_thresh1->Divide(h1_effDenom_fwd_gluon_thresh1);
  h1_eff_fwd_gluon_syst_thresh2->Divide(h1_effDenom_fwd_gluon_thresh2);
  h1_eff_fwd_gluon_syst_thresh3->Divide(h1_effDenom_fwd_gluon_thresh3);


  drawEffVsPt(db, discrim, thresh1, thresh1, 0., 2., h1_eff_centr_gluon_thresh1, h1_eff_centr_gluon_syst_thresh1, h1_eff_centr_quark_thresh1, h1_eff_centr_quark_syst_thresh1);
  drawEffVsPt(db, discrim, thresh1, thresh1, 3., 5., h1_eff_fwd_gluon_thresh1, h1_eff_fwd_gluon_syst_thresh1, h1_eff_fwd_quark_thresh1, h1_eff_fwd_quark_syst_thresh1);


  outfile->cd();
 
  h1_ptZ->Write();
  h1_rho->Write();

  h1_data_ptZ->Write();
  h1_data_rho->Write();

  h1_qglJet_pt3050_eta02_rho040->Write();
  h1_qglJetSyst_pt3050_eta02_rho040->Write();
  
  h1_qglJet_pt5065_eta02_rho040->Write();
  h1_qglJetSyst_pt5065_eta02_rho040->Write();
  
  h1_qglJet_pt80120_eta02_rho040->Write();
  h1_qglJetSyst_pt80120_eta02_rho040->Write();

  h1_qglJet_pt120250_eta02_rho040->Write();
  h1_qglJetSyst_pt120250_eta02_rho040->Write();

  h1_qglJet_pt3050_eta35_rho040->Write();
  h1_qglJetSyst_pt3050_eta35_rho040->Write();
  
  h1_qglJet_pt5065_eta35_rho040->Write();
  h1_qglJetSyst_pt5065_eta35_rho040->Write();
  
  h1_qglJet_pt80120_eta35_rho040->Write();
  h1_qglJetSyst_pt80120_eta35_rho040->Write();


  h1_data_qglJet_pt3050_eta02_rho040->Write();
  
  h1_data_qglJet_pt5065_eta02_rho040->Write();
  
  h1_data_qglJet_pt80120_eta02_rho040->Write();

  h1_data_qglJet_pt120250_eta02_rho040->Write();

  h1_data_qglJet_pt3050_eta35_rho040->Write();
  
  h1_data_qglJet_pt5065_eta35_rho040->Write();
  
  h1_data_qglJet_pt80120_eta35_rho040->Write();

  outfile->Close();
  


  return 0;
  
}


void drawSinglePlot( DrawBase* db, const std::string& discrim, TH1D* h1_data, TH1D* h1_qgl, TH1D* h1_qglSyst, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax ) {


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  h1_qgl->SetLineStyle(2);
  h1_qgl->SetLineWidth(2);
  h1_qglSyst->SetLineWidth(2);
  h1_data->SetMarkerStyle(20);

  float ymax = h1_data->GetMaximum();
  if( discrim=="qgMLPJet" )
    ymax *= 1.5;
  else
    ymax *= 1.3;

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., ymax);
  if( discrim=="qgMLPJet" ) 
    h2_axes->SetXTitle("Quark-Gluon MLP");
  else
    h2_axes->SetXTitle("Quark-Gluon Likelihood");
  h2_axes->SetYTitle("Events");

  h1_qglSyst->SetLineColor(kRed);

  float mcNorm = h1_data->Integral()/h1_qgl->Integral();
  h1_qgl->Scale(mcNorm);
  float mcSystNorm = h1_data->Integral()/h1_qglSyst->Integral();
  h1_qglSyst->Scale(mcSystNorm);

  h2_axes->Draw();

  char legendTitle[500];
  //sprintf( legendTitle, "p_{T}(%.0f-%.0f), #eta(%.1f-%.1f), #rho(%.0f-%.0f)", ptMin, ptMax, etaMin, etaMax, rhoMin, rhoMax);
  if( etaMin==0. ) 
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV,  |#eta| < %.0f", ptMin, ptMax, etaMax);
  else 
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV,  %.0f < |#eta| < %.0f", ptMin, ptMax, etaMin, etaMax);


  TLegend* legend = new TLegend(0.25, 0.66, 0.7, 0.9, legendTitle);
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_data, "Data", "p" );
  legend->AddEntry( h1_qgl, "MC Before Smearing", "L" );
  legend->AddEntry( h1_qglSyst, "MC After Smearing", "L" );
  legend->Draw("same");


  h1_qgl->Draw("Histo same");
  h1_qglSyst->Draw("Histo same");
  h1_data->Draw("P same");

  std::string doubleMin_str = "";
  if( doubleMin ) doubleMin_str = "DM";

  char canvasName[1023];
  sprintf( canvasName, "%s/checkSyst%s_%s_pt%.0f_%.0f_eta%.0f_%.0f_rho%.0f_%.0f.eps", db->get_outputdir().c_str(), discrim.c_str(), doubleMin_str.c_str(), ptMin, ptMax, etaMin, etaMax, rhoMin, rhoMax);

  TPaveText* label_top = db->get_labelTop();
  label_top->Draw("same");

  c1->SaveAs(canvasName);

  delete c1;
  delete legend;
  delete h2_axes;

}



void drawEffVsPt( DrawBase* db, const std::string& discrim, float thresh1, float thresh3, float etaMin, float etaMax, TH1D* h1_eff_gluon_thresh1, TH1D* h1_eff_gluon_syst_thresh1, TH1D* h1_eff_quark_thresh3, TH1D* h1_eff_quark_syst_thresh3) {


  float ptMax = 250.;
  char etaText[300];

  if( etaMin>0. )
    sprintf( etaText, "%.0f < |#eta| < %.0f", etaMin, etaMax );
  else
    sprintf( etaText, "|#eta| < %.0f", etaMax );


  if( etaMin>2.5 )
    ptMax = 135.;



  TH2D* h2_axes = new TH2D("axes", "", 10, 20., ptMax, 10, 0., 1.);
  h2_axes->SetXTitle("Jet p_{T} [Gev]");
  h2_axes->SetYTitle("Efficiency");
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetXaxis()->SetNoExponent();


  h1_eff_gluon_thresh1->SetMarkerStyle(20);
  h1_eff_gluon_thresh1->SetMarkerSize(1.6);
  h1_eff_gluon_thresh1->SetMarkerColor(46);

  h1_eff_gluon_syst_thresh1->SetMarkerStyle(24);
  h1_eff_gluon_syst_thresh1->SetMarkerSize(1.6);
  h1_eff_gluon_syst_thresh1->SetMarkerColor(46);

  h1_eff_quark_thresh3->SetMarkerStyle(21);
  h1_eff_quark_thresh3->SetMarkerSize(1.6);
  h1_eff_quark_thresh3->SetMarkerColor(38);

  h1_eff_quark_syst_thresh3->SetMarkerStyle(25);
  h1_eff_quark_syst_thresh3->SetMarkerSize(1.6);
  h1_eff_quark_syst_thresh3->SetMarkerColor(38);

  TPaveText* labelTop = db->get_labelTop();
  TPaveText* labelEta = new TPaveText(0.7, 0.8, 0.9, 0.9, "brNDC");
  labelEta->AddText(etaText);

  std::string discrim_text = (discrim=="qgMLPJet") ? "MLP" : "LD";
  std::string operator_discrim = (discrim=="qgMLPJet") ? "<" : ">";
  char gluon_text[500];
  sprintf( gluon_text, "Gluon Eff. (%s %s %.1f)", discrim_text.c_str(), operator_discrim.c_str(), thresh1 );
  char gluon_syst_text[500];
  sprintf( gluon_syst_text, "Gluon Eff. (%s %s %.1f) After Syst.", discrim_text.c_str(), operator_discrim.c_str(), thresh1 );
  char quark_text[500];
  sprintf( quark_text, "Quark Eff. (%s %s %.1f)", discrim_text.c_str(), operator_discrim.c_str(), thresh3 );
  char quark_syst_text[500];
  sprintf( quark_syst_text, "Quark Eff. (%s %s %.1f) After Syst.", discrim_text.c_str(), operator_discrim.c_str(), thresh3 );


  char legendTitle[300];
  sprintf( legendTitle, "Quark-Gluon %s, %s", discrim_text.c_str(), etaText );

  TLegend* legend = new TLegend( 0.18, 0.18, 0.5, 0.43, legendTitle );
  legend->SetTextSize(0.035);
  legend->SetFillColor(kWhite);
  legend->AddEntry( h1_eff_quark_thresh3, quark_text, "P");
  legend->AddEntry( h1_eff_quark_syst_thresh3, quark_syst_text, "P");
  legend->AddEntry( h1_eff_gluon_thresh1, gluon_text, "P");
  legend->AddEntry( h1_eff_gluon_syst_thresh1, gluon_syst_text, "P");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  c1->SetLogx();
  h2_axes->Draw();
  legend->Draw("same");
  labelTop->Draw("same");
  //labelEta->Draw("same");
  h1_eff_gluon_thresh1->Draw("p same");
  h1_eff_gluon_syst_thresh1->Draw("p same");
  h1_eff_quark_thresh3->Draw("p same");
  h1_eff_quark_syst_thresh3->Draw("p same");
  
  char canvasName[500];
  if( etaMin > 2.5 ) {
    sprintf( canvasName, "%s/syst_eff_%s_fwd.eps", db->get_outputdir().c_str(), discrim.c_str() );
  } else {
    sprintf( canvasName, "%s/syst_eff_%s_centr.eps", db->get_outputdir().c_str(), discrim.c_str() );
  }
  
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete labelEta;
  delete legend;

}



void drawMC_beforeAfter( DrawBase* db, const std::string& discrim, TH1D* h1_qglJet_quark, TH1D* h1_qglJetSyst_quark, TH1D* h1_qglJet_gluon, TH1D* h1_qglJetSyst_gluon, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax ) {


  h1_qglJet_quark->SetLineColor(38);
  h1_qglJet_quark->SetLineWidth(2);
  h1_qglJet_quark->SetFillColor(38);
  h1_qglJet_quark->SetFillStyle(3004);

  h1_qglJetSyst_quark->SetMarkerColor(38);
  h1_qglJetSyst_quark->SetMarkerStyle(20);
  h1_qglJetSyst_quark->SetMarkerSize(1.3);

  h1_qglJet_gluon->SetLineColor(46);
  h1_qglJet_gluon->SetLineWidth(2);
  h1_qglJet_gluon->SetFillColor(46);
  h1_qglJet_gluon->SetFillStyle(3005);

  h1_qglJetSyst_gluon->SetMarkerColor(46);
  h1_qglJetSyst_gluon->SetMarkerStyle(21);
  h1_qglJetSyst_gluon->SetMarkerSize(1.3);


  float histoMax = h1_qglJet_quark->GetMaximum()/h1_qglJet_quark->Integral();
  if( h1_qglJetSyst_quark->GetMaximum()/h1_qglJetSyst_quark->Integral() > histoMax ) histoMax = h1_qglJetSyst_quark->GetMaximum()/h1_qglJetSyst_quark->Integral();
  if( h1_qglJet_gluon->GetMaximum()/h1_qglJet_gluon->Integral() > histoMax ) histoMax = h1_qglJet_gluon->GetMaximum()/h1_qglJet_gluon->Integral();
  if( h1_qglJetSyst_gluon->GetMaximum()/h1_qglJetSyst_gluon->Integral() > histoMax ) histoMax = h1_qglJetSyst_gluon->GetMaximum()/h1_qglJetSyst_gluon->Integral();

  float yMax = 1.3*histoMax;


  std::string discrim_text = (discrim=="qgMLPJet") ? "Quark-Gluon MLP" : "Quark-Gluon Likelihood";

  TH2D* h2_axes = new TH2D("axes", "", 10, h1_qglJet_quark->GetXaxis()->GetXmin(), h1_qglJet_quark->GetXaxis()->GetXmax(), 10, 0., yMax);
  h2_axes->SetXTitle(discrim_text.c_str());
  h2_axes->SetYTitle("Normalized to Unity");

  char legendTitle[500];
  if( etaMin==0. )
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV,  |#eta| < %.0f", ptMin, ptMax, etaMax);
  else
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV,  %.0f < |#eta| < %.0f", ptMin, ptMax, etaMin, etaMax);

  TLegend* legend = new TLegend(0.25, 0.63, 0.7, 0.9, legendTitle);
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_qglJet_quark, "Quark MC", "F");
  legend->AddEntry( h1_qglJetSyst_quark, "Quark MC After Smearing", "P");
  legend->AddEntry( h1_qglJet_gluon, "Gluon MC", "F");
  legend->AddEntry( h1_qglJetSyst_gluon, "Gluon MC After Smearing", "P");

  TPaveText* labelTop = db->get_labelTop();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  h2_axes->Draw();
  legend->Draw("same");
  labelTop->Draw("same");

  h1_qglJet_gluon->DrawNormalized("histo same");
  h1_qglJet_quark->DrawNormalized("histo same");
  h1_qglJetSyst_gluon->DrawNormalized("p same");
  h1_qglJetSyst_quark->DrawNormalized("p same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/%s_MC_beforeAfter_pt%.0f_%.0f_eta%.0f_%.0f_rho%.0f_%.0f.eps", db->get_outputdir().c_str(), discrim.c_str(), ptMin, ptMax, etaMin, etaMax, rhoMin, rhoMax);

  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;

}
