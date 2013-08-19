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




struct FitParameters {

  float q_a;
  float q_b;
  float g_a;
  float g_b;
  float chi2;

};


std::string scanSingleBin( const std::string& selection, TTree* tree, TTree* tree_data, const std::string& discrim, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin=0., float rhoMax=100. );
FitParameters getSingleChiSquare( const std::string& selection, TTree* tree, TTree* tree_data, const std::string& discrim, float q_a, float q_b, float g_a, float g_b, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin=0., float rhoMax=100. );
float computeChiSquare( TH1D* h1_data, TH1D* h1_mc );




int main( int argc, char* argv[] ) {


  if( argc==1 ) {
    std::cout << "USAGE: ./checkSyst [selection] [discriminator=\"qglJet\"]" << std::endl;
    std::cout << "[selection] can be \"ZJet\" or \"DiJet\"" << std::endl;
    exit(1);
  }


  std::string selection = "ZJet";
  if( argc>1 ) {
    std::string selection_str(argv[1]);
    if( selection_str!="ZJet" && selection_str!="DiJet" ) {
      std::cout << "-> Unknown selection. Only \"ZJet\" and \"DiJet\" allowed. Exiting." << std::endl;
      exit(191);
    }
    selection = selection_str;
  }

  std::string discrim = "qglJet";
  if( argc>2 ) {
    std::string discrim_str(argv[2]);
    if( discrim_str!="qglJet" && discrim_str!="qgMLPJet" ) {
      std::cout << "-> Unknown discriminator. Only \"qglJet\" and \"qgMLPJet\" allowed. Exiting." << std::endl;
      exit(191);
    }
    discrim = discrim_str;
  }




  std::string mcFileName = (selection=="ZJet") ? "sunilFlat_ZJet_Zjets_12Jul.root" : "sunilFlat_DiJet_flatQCD_P6_Dijets_12Aug_ptHatWeight.root";
  std::cout << "-> Opening MC file: " << mcFileName << std::endl;
  TFile* file = TFile::Open(mcFileName.c_str());
  TTree* tree = (TTree*)file->Get("tree_passedEvents");

  TChain* tree_data = new TChain("tree_passedEvents");
  std::string dataFileName = (selection=="ZJet") ? "sunilFlat_ZJet_data2012ABCD_MuPD_12Jul.root" : "sunilFlat_DiJet_data2012ABCD_MBPD_12Jul.root";
  std::cout << "-> Opening data file: " << dataFileName << std::endl;
  tree_data->Add(dataFileName.c_str());


  std::string outfilename = "SystDB_new_"+selection+".txt";
  ofstream systdbfile(outfilename.c_str());
  systdbfile << scanSingleBin( selection, tree, tree_data, discrim, 20., 30., 3., 4.7 ) << std::endl;
  //systdbfile << scanSingleBin( selection, tree, tree_data, discrim, 30., 40., 3., 4.7 ) << std::endl;
  //systdbfile << scanSingleBin( selection, tree, tree_data, discrim, 40., 50., 3., 4.7 ) << std::endl;
  //systdbfile << scanSingleBin( selection, tree, tree_data, discrim, 50., 65., 3., 4.7 ) << std::endl;
  systdbfile.close();

  return 0;

}




std::string scanSingleBin( const std::string& selection, TTree* tree, TTree* tree_data, const std::string& discrim, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax ) {

  std::cout << std::endl;
  std::cout << "-> Beginning scan for new bin:" << std::endl;
  std::cout << "     pt: " << ptMin << "-" << ptMax << std::endl;
  std::cout << "     eta: " << etaMin << "-" << etaMax << std::endl;
  std::cout << "     rho: " << rhoMin << "-" << rhoMax << std::endl;
  std::cout << std::endl;
  std::cout << "-> First round: smear only gluons." << std::endl;


  float chi2min = 100000.;
  float g_a_min = 1.;
  float g_b_min = 0.;

  for( float g_a = 0.8; g_a<1.04; g_a+=0.02 ) {

    for( float g_b = -1.0; g_b<0.5; g_b+=0.1 ) {

      std::cout << "Scanning: g_a = " << g_a << "  g_b = " << g_b << std::endl;

      FitParameters fp = getSingleChiSquare( selection, tree, tree_data, discrim, 1., 0., g_a, g_b, ptMin, ptMax, etaMin, etaMax );
      
      if( fp.chi2 < chi2min ) {
        chi2min = fp.chi2;
        g_a_min = fp.g_a;
        g_b_min = fp.g_b;
      }

    } // for g_b

  } // for g_a

  std::cout << "-> Found chi2 minimum with g_a = " << g_a_min << "  g_b = " << g_b_min << std::endl;

  std::cout << std::endl;
  std::cout << "-> Second round: smear only quarks." << std::endl;

  chi2min = 100000.;
  float q_a_min = 1.;
  float q_b_min = 0.;

  for( float q_a = 0.965; q_a<1.015; q_a+=0.005 ) {

    for( float q_b = -0.2; q_b<0.2; q_b+=0.1 ) {

      std::cout << "Scanning: q_a = " << q_a << "  q_b = " << q_b << std::endl;

      FitParameters fp = getSingleChiSquare( selection, tree, tree_data, discrim, q_a, q_b, g_a_min, g_b_min, ptMin, ptMax, etaMin, etaMax );
 
      if( fp.chi2 < chi2min ) {
        chi2min = fp.chi2;
        q_a_min = fp.q_a;
        q_b_min = fp.q_b;
      }

    } // for g_b

  } // for g_a

  std::cout << "-> Found chi2 minimum with q_a = " << q_a_min << "  q_b = " << q_b_min << std::endl;


  std::cout << std::endl;
  std::cout << "-> Last round: last correction to gluons" << std::endl;
  chi2min = 100000.;
  float g_a_min2 = 1.;
  float g_b_min2 = 0.;

  for( float corr_g_a = 0.98; corr_g_a<1.02; corr_g_a+=0.01 ) {

    for( float corr_g_b = 0.98; corr_g_b<1.02; corr_g_b+=0.01 ) {

      float g_a = g_a_min*corr_g_a;
      float g_b = g_b_min*corr_g_b;

      std::cout << "Scanning: g_a = " << g_a << "  g_b = " << g_b << std::endl;

      FitParameters fp = getSingleChiSquare( selection, tree, tree_data, discrim, q_a_min, q_b_min, g_a, g_b, ptMin, ptMax, etaMin, etaMax );
      
      if( fp.chi2 < chi2min ) {
        chi2min = fp.chi2;
        g_a_min2 = fp.g_a;
        g_b_min2 = fp.g_b;
      }

    } // for g_b

  } // for g_a

  std::cout << "-> Done. Total minimum found in q_a = " << q_a_min << "  q_b = " << q_b_min << "  g_a = " << g_a_min << "  g_b = " << g_b_min << std::endl;
  std::cout << "QGLHisto " << ptMin << " " << ptMax << " " << rhoMin << " " << rhoMax << " " << etaMin << " " << etaMax << " " << q_a_min << " " << q_b_min << " " << g_a_min << " " << g_b_min << " 0. 1." << std::endl;
  std::cout << std::endl;
  char line[1023];
  sprintf( line, "QGLHisto %.0f %.0f %.0f %.0f %.1f %.1f %.3f %.3f %.3f %.3f 0. 1.", ptMin,  ptMax,  rhoMin, rhoMax, etaMin, etaMax, q_a_min, q_b_min, g_a_min, g_b_min );

  std::string line_str(line);
  return line_str;

}




FitParameters getSingleChiSquare( const std::string& selection, TTree* tree, TTree* tree_data, const std::string& discrim, float q_a, float q_b, float g_a, float g_b, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax ) {

  float eventWeight;
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

  tree->SetBranchAddress("eventWeight", &eventWeight);
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



  std::string syst_tmp_name = "systdb_tmp_" + selection + ".txt";
  ofstream systdb(syst_tmp_name.c_str());
  float lmin = 0.;
  float lmax = 1.;
  systdb << "QGLHisto " << ptMin << " " << ptMax << " " << rhoMin << " " << rhoMax << " " << etaMin << " " << etaMax << " " << q_a << " " << q_b << " " << g_a << " " << g_b << " " << lmin << " " << lmax << std::endl;
  systdb.close();

  QGSyst qgsyst;
  qgsyst.ReadDatabaseDoubleMin("systdb_tmp.txt");

  std::string nameForSyst = (discrim=="qglJet") ? "QGLHisto" : "QGLMLP";
  qgsyst.SetTagger(nameForSyst);


  int nBins = 30;

                                                             
  TH1D* h1_mc = new TH1D("mc", "", nBins, 0., 1.0001);
  h1_mc->Sumw2();
  TH1D* h1_mcSyst = new TH1D("mcSyst", "", nBins, 0., 1.0001);
  h1_mcSyst->Sumw2();



  int nentries = tree->GetEntries();

  for( unsigned iEntry = 0; iEntry < nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( njet<1 ) continue;
    if( selection=="DiJet" && njet<2 ) continue;

    if( rho<rhoMin || rho>rhoMax ) continue;


    bool smear1=true;
    std::string type = "all";
    if( fabs(pdgId[0])>0 && fabs(pdgId[0])<6 ) type = "quark";
    else if( pdgId[0]==21 ) type = "gluon";
    else { // both 0 and -999
      //smearJet = false;
      smear1=false;
      //type = "gluon";
      //std::cout << "Unknown jet PDG ID (" << pdgIdPartJet << "). Will use gluon instead." << std::endl;
    }


    int tag_index = (selection=="DiJet") ? 1 : 0;

    if( pt[tag_index]>ptMin && pt[tag_index]<ptMax ) {
      if( fabs(eta[0])>=etaMin && fabs(eta[0])<etaMax ) {
        h1_mc->Fill( qglJet[0], eventWeight );
        if( smear1 )
          h1_mcSyst->Fill( qgsyst.Smear(pt[0], eta[0], rho, qglJet[0], type), eventWeight );
        else
          h1_mcSyst->Fill( qglJet[0], eventWeight );
      }
    }

    if( selection=="DiJet" ) {

      bool smear2=true;
      std::string type2 = "all";
      if( fabs(pdgId[1])>0 && fabs(pdgId[1])<6 ) type2 = "quark";
      else if( pdgId[1]==21 ) type2 = "gluon";
      else { // both 0 and -999
        //smearJet = false;
        smear2=false;
        //type = "gluon";
        //std::cout << "Unknown jet PDG ID (" << pdgIdPartJet << "). Will use gluon instead." << std::endl;
      }

      if( pt[0]>ptMin && pt[0]<ptMax ) {
        if( fabs(eta[1])>=etaMin && fabs(eta[1])<etaMax ) {
          h1_mc->Fill( qglJet[1], eventWeight );
          if( smear2 )
            h1_mcSyst->Fill( qgsyst.Smear(pt[1], eta[1], rho, qglJet[1], type2), eventWeight );
          else
            h1_mcSyst->Fill( qglJet[1], eventWeight );
        }
      }

    } // if dijet

  } //for mc


  // now switch to data:
  int event;

  tree_data->SetBranchAddress("event", &event );
  tree_data->SetBranchAddress("nJet", &njet);
  tree_data->SetBranchAddress("ptJet", pt);
  tree_data->SetBranchAddress("etaJet", eta);
  tree_data->SetBranchAddress("rhoPF", &rho);
  tree_data->SetBranchAddress("ptD_QCJet", ptD_QC);
  tree_data->SetBranchAddress(discrim.c_str(), qglJet);


  TH1D* h1_data = new TH1D("data", "", nBins, 0., 1.0001);


  int nentries_data = tree_data->GetEntries();

  for( unsigned iEntry = 0; iEntry < nentries_data; ++iEntry ) {

    tree_data->GetEntry(iEntry);


    if( njet<1 ) continue;
    if( selection=="DiJet" && njet<2 ) continue;

    if( rho<rhoMin || rho>rhoMax ) continue;


    int tag_index = (selection=="DiJet") ? 1 : 0;

    if( pt[tag_index]>ptMin && pt[tag_index]<ptMax ) {
      if( fabs(eta[0])>=etaMin && fabs(eta[0])<etaMax ) {
        h1_data->Fill( qglJet[0] );
      }
    }

    if( selection=="DiJet" ) {

      if( pt[0]>ptMin && pt[0]<ptMax ) {
        if( fabs(eta[1])>=etaMin && fabs(eta[1])<etaMax ) {
          h1_data->Fill( qglJet[1] );
        }
      }

    } // if dijet

  } //for data


  float thisChiSquare = computeChiSquare( h1_data, h1_mcSyst );

  FitParameters fp;
  fp.q_a = q_a;
  fp.q_b = q_b;
  fp.g_a = g_a;
  fp.g_b = g_b;
  fp.chi2 = thisChiSquare;


  delete h1_data;
  delete h1_mc;
  delete h1_mcSyst;

  return fp;

}



float computeChiSquare( TH1D* h1_data, TH1D* h1_mc ) {

  float mc_integral = h1_mc->Integral();
  float data_integral = h1_data->Integral();
  h1_mc->Scale( data_integral/mc_integral );

  float chi2 = 0.;

  for( unsigned ibin=1; ibin<h1_data->GetNbinsX(); ++ibin ) {

    float data = h1_data->GetBinContent(ibin);

    if( data>5 ) { //consider only bins with 5 or more entries

      float mc = h1_mc->GetBinContent(ibin);
      float diffSquared = (data-mc)*(data-mc);
      chi2 += diffSquared/data; // error on data is sqrt(data) so error squared is data

    }

  } // for bins

  
  return chi2;

}
