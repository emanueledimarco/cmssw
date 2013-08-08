#include "DrawBase.h"
#include <cstdlib>
#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TChain.h"




void drawHistoWithQuarkGluonComponents( const std::string& selectionType, DrawBase* db, const std::string& treeName, const std::string& additionalCuts, const std::string& varName, const std::string& canvasSaveName, std::string axisName, const std::string& units="", float ptMin=0., float ptMax=10000., float etaMin=0., float etaMax =10., float rhoMin=0., float rhoMax=100., int nBins=30, float xMin=0., float xMax=1.0001, bool legendQuadrant=1, bool log_aussi=false );



int main(int argc, char* argv[]) {


  std::string selectionType;
  if( argc>1 ) {
    std::string selectionType_str(argv[1]);
    selectionType = selectionType_str;
  }

  if( selectionType!="ZJets" && selectionType!="DiJets" ) {
    std::cout << "Supported selections are only \"ZJets\" and \"DiJets\". Exiting." << std::endl;
    exit(11);
  }


  TFile* file_data = TFile::Open("sunilFlat_ZJet_data2012ABCD_MuPD_12Jul.root");
  TFile* file_mc = TFile::Open("sunilFlat_ZJet_Zjets_12Jul.root");

  DrawBase* db = new DrawBase("qgdatamc");

  std::string outputdir = "QGLDataMCPlots_" + selectionType;
  db->set_outputdir(outputdir);

  db->add_dataFile(file_data, "data");
  db->add_mcFile(file_mc, "mc", "mc_process");

  drawHistoWithQuarkGluonComponents( selectionType, db, "tree_passedEvents", "", "qglJet[0]", "qglJet", "Quark-Gluon Likelihood Discriminator", "", 40., 50., 3., 4.7);
  drawHistoWithQuarkGluonComponents( selectionType, db, "tree_passedEvents", "axis1_QCJet[0]*axis1_QCJet[0]+axis2_QCJet[0]*axis2_QCJet[0]<0.06", "qglJet[0]", "qglJet_PUID", "Quark-Gluon Likelihood Discriminator", "", 40., 50., 3., 4.7);

  return 0;

}












void drawHistoWithQuarkGluonComponents( const std::string& selectionType, DrawBase* db, const std::string& treeName, const std::string& additionalCuts, const std::string& varName, const std::string& canvasSaveName, std::string axisName, const std::string& units, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax, int nBins, float xMin, float xMax, bool legendQuadrant, bool log_aussi ) {


  TString varName_tstr(varName);



  TH1D* h1_data = new TH1D( "data", "", nBins, xMin, xMax );

  TH1D* h1_all = new TH1D( "all", "", nBins, xMin, xMax );
  TH1D* h1_quark = new TH1D( "quark", "", nBins, xMin, xMax );
  TH1D* h1_gluon = new TH1D( "gluon", "", nBins, xMin, xMax );
  TH1D* h1_pu = new TH1D( "pu", "", nBins, xMin, xMax );
  TH1D* h1_b = new TH1D( "b", "", nBins, xMin, xMax );

  //// these ones for all processes (to get the fractions right):
  //TH1D* h1_all_all = new TH1D( "all_all", "", nBins, xMin, xMax );
  //TH1D* h1_quark_all = new TH1D( "quark_all", "", nBins, xMin, xMax );
  //TH1D* h1_gluon_all = new TH1D( "gluon_all", "", nBins, xMin, xMax );
  //TH1D* h1_pu_all = new TH1D( "pu_all", "", nBins, xMin, xMax );
  //TH1D* h1_b_all = new TH1D( "b_all", "", nBins, xMin, xMax );



  char commonCondition[500];
  if( additionalCuts!="" )
    sprintf( commonCondition, "%s && ptJet[0]>%f && ptJet[0]<%f && abs(etaJet[0])>=%f && abs(etaJet[0])<%f && rhoPF>%f && rhoPF<%f", additionalCuts.c_str(), ptMin, ptMax, etaMin, etaMax, rhoMin, rhoMax ); 
    //sprintf( commonCondition, "%s && ptJet0>%f && ptJet0<%f && QGLikelihoodJet0>0. && QGLikelihoodJet0<1. && rhoPF>%f && rhoPF<%f", additionalCuts.c_str(), ptMin, ptMax, rhoMin, rhoMax ); 
  else
    sprintf( commonCondition, "ptJet[0]>%f && ptJet[0]<%f && abs(etaJet[0])>=%f && abs(etaJet[0])<%f && rhoPF>%f && rhoPF<%f", ptMin, ptMax, etaMin, etaMax, rhoMin, rhoMax ); 
    //sprintf( commonCondition, "ptJet0>%f && ptJet0<%f && QGLikelihoodJet0>0. && QGLikelihoodJet0<1. && rhoPF>%f && rhoPF<%f", ptMin, ptMax, rhoMin, rhoMax ); 
  //sprintf( commonCondition, "ptJet0>%f && ptJet0<%f", ptMin, ptMax ); 


  char allCondition[800];
  sprintf( allCondition,   "eventWeight*(%s)", commonCondition );
  char quarkCondition[800];
  sprintf( quarkCondition, "eventWeight*(abs(pdgIdJet[0])<5  && abs(pdgIdJet[0])!=0 && %s)", commonCondition );
  char gluonCondition[800];
  sprintf( gluonCondition, "eventWeight*(pdgIdJet[0]==21     && abs(pdgIdJet[0])!=0 && %s)", commonCondition );
  char bCondition[800];
  sprintf( bCondition,     "eventWeight*(abs(pdgIdJet[0])==5 && abs(pdgIdJet[0])!=0 && %s)", commonCondition );
  char puCondition[800];
  sprintf( puCondition,     "eventWeight*(pdgIdJet[0]==0 && %s)", commonCondition );

  TTree* treeDATA = (TTree*)(db->get_dataFile(0).file->Get(treeName.c_str()));
  treeDATA->Project( "data", varName.c_str(), commonCondition );


  TChain* treeMC= new TChain(treeName.c_str());
  //// this one to get the shapes (avoid huge QCD weights for gamma+jet):
  //TChain* treeMC_signal = new TChain(treeName.c_str());
  for( unsigned iFile=0; iFile<db->get_mcFiles().size(); ++iFile ) {
    std::string fileName(db->get_mcFile(iFile).file->GetName());
    std::string treeFullName = fileName + "/" + treeName;
    treeMC->Add(treeFullName.c_str());
    //if( iFile==0 ) //signal only
    //  treeMC_signal->Add(treeFullName.c_str());
  }

  treeMC->Project( "all",   varName.c_str(), allCondition );
  treeMC->Project( "quark", varName.c_str(), quarkCondition );
  treeMC->Project( "gluon", varName.c_str(), gluonCondition );
  treeMC->Project( "pu", varName.c_str(), puCondition );
  treeMC->Project( "b", varName.c_str(), bCondition );

  //treeMC_all->Project( "all_all",   varName.c_str(), allCondition );
  //treeMC_all->Project( "quark_all", varName.c_str(), quarkCondition );
  //treeMC_all->Project( "gluon_all", varName.c_str(), gluonCondition );
  //treeMC_all->Project( "pu_all", varName.c_str(), puCondition );
  //treeMC_all->Project( "b_all", varName.c_str(), bCondition );

  float data_int = h1_data->Integral();
  float mc_int = h1_all->Integral();
  //float mc_int_all = h1_all_all->Integral();
  float scaleFactor = data_int/mc_int;

  float quark_fraction = h1_quark->Integral()/mc_int;
  float gluon_fraction = h1_gluon->Integral()/mc_int;
  float pu_fraction = h1_pu->Integral()/mc_int;
  float b_fraction = h1_b->Integral()/mc_int;
  float other_fraction = 1.-quark_fraction-gluon_fraction-b_fraction;


  //h1_all->Scale( h1_all_all->Integral()/h1_all->Integral() );
  //h1_gluon->Scale( h1_gluon_all->Integral()/h1_gluon->Integral() );
  //h1_pu->Scale( h1_pu_all->Integral()/h1_pu->Integral() );
  //h1_quark->Scale( h1_quark_all->Integral()/h1_quark->Integral() );
  //h1_b->Scale( h1_b_all->Integral()/h1_b->Integral() );
  

  char quarkText[300];
  sprintf( quarkText, "udsc");
  char gluonText[300];
  sprintf( gluonText, "Gluon");
  char bText[300];
  sprintf( bText, "b");
  char puText[300];
  sprintf( puText, "Pile Up");
  char otherText[300];
  sprintf( otherText, "Undefined");

  //char quarkText[300];
  //sprintf( quarkText, "udsc (%.1f%%)", 100.*quark_fraction );
  //char gluonText[300];
  //sprintf( gluonText, "Gluons (%.1f%%)", 100.*gluon_fraction );
  //char bText[300];
  //sprintf( bText, "b (%.1f%%)", 100.*b_fraction );
  //char puText[300];
  //sprintf( puText, "Pile Up (%.1f%%)", 100.*pu_fraction );
  //char otherText[300];
  //sprintf( otherText, "Undefined (%.1f%%)", 100.*other_fraction );


  float xMin_leg = 0.32;
  float xMax_leg = 0.8;

  if( legendQuadrant==1 ) {
    xMin_leg = 0.6;
    xMax_leg = 0.93;
  }
  
  //TLegend* legend;
  //if( (ptMin !=0. || ptMax != 10000.) && (rhoMin!=0. || rhoMax !=30.) ) {
  //  char legendTitle[250];
  //  if( varName=="QGLikelihoodJet0" ) {
  //    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV,  %.0f < #rho < %.0f GeV", ptMin, ptMax, rhoMin, rhoMax );
  //    legend = new TLegend( 0.32, 0.55, 0.8, 0.9, legendTitle );
  //  } else {
  //    sprintf( legendTitle, "#splitline{%.0f < p_{T} < %.0f GeV}{%.0f < #rho < %.0f GeV}", ptMin, ptMax, rhoMin, rhoMax );
  //    legend = new TLegend( xMin_leg, 0.5, xMax_leg, 0.9, legendTitle );
  //  }
  //} else if( ptMin !=0. && ptMax != 10000. ) {
  //  char legendTitle[150];
  //  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax);
  //  legend = new TLegend( xMin_leg, 0.55, xMax_leg, 0.9, legendTitle );
  //} else {
  //  legend = new TLegend( xMin_leg, 0.6, xMax_leg, 0.9 );
  //}

  std::string selectionType_text;
  if( selectionType=="ZJets" ) selectionType_text = "Z+Jets";
  if( selectionType=="DiJets" ) selectionType_text = "DiJets";

  TLegend* legend = new TLegend( xMin_leg, 0.55, xMax_leg, 0.91, selectionType_text.c_str() );
  legend->SetFillColor( kWhite );
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_data, "Data", "p" );
  legend->AddEntry( h1_quark, quarkText, "F" );
  legend->AddEntry( h1_gluon, gluonText, "F" );
  legend->AddEntry( h1_pu, puText, "F" );
  legend->AddEntry( h1_b, bText, "F" );
  legend->AddEntry( h1_all, otherText, "F" );

  TPaveText* cutsText = new TPaveText(0.22, 0.77, 0.45, 0.93, "brNDC");
  cutsText->SetTextSize(0.038);
  cutsText->SetTextFont(42);
  cutsText->SetFillColor(0);
  char cutsText_char[500];
  sprintf( cutsText_char, "#splitline{%.0f < p_{T} < %.0f GeV}{%.0f < |#eta| < %.0f}", ptMin, ptMax, etaMin, etaMax );
  cutsText->AddText(cutsText_char);

  h1_all->Rebin( db->get_rebin() );
  h1_gluon->Rebin( db->get_rebin() );
  h1_pu->Rebin( db->get_rebin() );
  h1_quark->Rebin( db->get_rebin() );
  h1_b->Rebin( db->get_rebin() );
  h1_data->Rebin( db->get_rebin() );
  
  h1_all->Scale( scaleFactor );
  h1_gluon->Scale( scaleFactor );
  h1_pu->Scale( scaleFactor );
  h1_quark->Scale( scaleFactor );
  h1_b->Scale( scaleFactor );
  
//  h1_all_all->Rebin( db->get_rebin() );
//  h1_all_gluon->Rebin( db->get_rebin() );
//  h1_all_pu->Rebin( db->get_rebin() );
//  h1_all_quark->Rebin( db->get_rebin() );
//  h1_all_b->Rebin( db->get_rebin() );
//  
//  h1_all_all->Scale( scaleFactor );
//  h1_all_gluon->Scale( scaleFactor );
//  h1_all_pu->Scale( scaleFactor );
//  h1_all_quark->Scale( scaleFactor );
//  h1_all_b->Scale( scaleFactor );
//  
  h1_data->SetMarkerStyle( 20 );
  h1_data->SetMarkerSize( 1. );
  h1_all->SetFillColor( kGray );
  h1_gluon->SetFillColor( 46 );
  h1_pu->SetFillColor( 30 );
  h1_quark->SetFillColor( 38 );
  h1_b->SetFillColor( kYellow );

//  h1_all_all->SetFillColor( kGray );
//  h1_all_gluon->SetFillColor( 46 );
//  h1_all_pu->SetFillColor( 30 );
//  h1_all_quark->SetFillColor( 38 );
//  h1_all_b->SetFillColor( kYellow );
//
  THStack* stack = new THStack();
  stack->Add(h1_pu);
  stack->Add(h1_b);
  stack->Add(h1_gluon );
  stack->Add(h1_quark);

//  THStack* stack_all = new THStack();
//  stack_all->Add(h1_all_gluon );
//  stack_all->Add(h1_all_quark);
//  stack_all->Add(h1_all_pu);
//  stack_all->Add(h1_all_b);
//
  float dataMax = h1_data->GetMaximum();
  float mcMax = h1_all->GetMaximum();
  float yMax = (dataMax>mcMax) ? dataMax : mcMax;
  yMax *= db->get_yAxisMaxScale();


  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();

  char yAxisTitle[200];
  std::string units_text = (units!="") ? (" "+units) : "";
  if( (h1_data->GetBinWidth(1)) < 0.1 )
    sprintf( yAxisTitle, "Events / (%.2f%s)", h1_data->GetBinWidth(1), units_text.c_str() );
  else if( ((int)(10.*h1_data->GetBinWidth(1)) % 10) == 0 )
    sprintf( yAxisTitle, "Events / (%.0f%s)", h1_data->GetBinWidth(1), units_text.c_str() );
  else
    sprintf( yAxisTitle, "Events / (%.1f%s)", h1_data->GetBinWidth(1), units_text.c_str() );


  if( units!="" ) axisName = axisName + " [" + units + "]";


  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 220. );
  //TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( yAxisTitle );
  if( yMax>1000. )
    h2_axes->GetYaxis()->SetTitleOffset(1.55); 

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();
  legend->Draw("same");
  cutsText->Draw("same");
  h1_all->Draw("same");
  stack->Draw("histo same");
  h1_data->Draw("e same");
  sqrtLabel->Draw("Same");

  gPad->RedrawAxis();

  //std::string canvasName = db->get_outputdir() + "/" + varName + "_components.eps";

  char canvasNameChar[400];
  if( etaMax>4. )
    sprintf( canvasNameChar, "%s/%s_pt%d%d_fwd.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax );
  else if( etaMax>2. )
    sprintf( canvasNameChar, "%s/%s_pt%d%d_trans.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax );
  else 
    sprintf( canvasNameChar, "%s/%s_pt%d%d_centr.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax );

  c1->SaveAs(canvasNameChar);
  std::string epstopdf_command(canvasNameChar);
  epstopdf_command = "epstopdf " + epstopdf_command;
  system(epstopdf_command.c_str());

  if( log_aussi ) {

    c1->Clear();
    c1->SetLogy();

    float ymin_log = 0.1;
    if( varName=="betaStarJet0" ) ymin_log = 0.01;

    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, ymin_log, yMax*10 );
    h2_axes_log->SetXTitle( axisName.c_str() );
    h2_axes_log->SetYTitle( yAxisTitle );


    h2_axes_log->Draw();
    legend->Draw("same");
    h1_all->Draw("same");
    stack->Draw("histo same");
    h1_data->Draw("e same");
    sqrtLabel->Draw("Same");

    gPad->RedrawAxis();

    //std::string canvasName = db->get_outputdir() + "/" + varName + "_components.eps";

    char canvasNameChar_log[400];
    if( rhoMin==0. && rhoMax==30. )
      sprintf( canvasNameChar_log, "%s/%s_pt%d%d.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax );
    else
      sprintf( canvasNameChar_log, "%s/%s_pt%d%d_rho%d%d.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax, (int)rhoMin, (int)rhoMax );

    c1->SaveAs(canvasNameChar_log);

    delete h2_axes_log;

  }
  
  delete c1;
  delete h2_axes;

  delete h1_data;
  delete h1_all;
  delete h1_quark;
  delete h1_gluon;
  delete h1_pu;
  delete h1_b;
  
  //delete h1_all_all;
  //delete h1_quark_all;
  //delete h1_gluon_all;
  //delete h1_pu_all;
  //delete h1_b_all;
  
}

