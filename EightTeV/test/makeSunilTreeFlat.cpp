#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "QuarkGluonTagger/EightTeV/interface/QGLikelihoodCalculator.h"
#include "/afs/cern.ch/work/p/pandolf/CMSSW_5_3_9_patch1_QGSyst/src/QuarkGluonTagger/EightTeV/interface/QGMLPCalculator.h"



float get_xSecWeight( const std::string& dataset );
float correctNCharged(int nChg_QCJet, float eta);


int main( int argc, char* argv[] ) {


  if( argc<3 ) {
    std::cout << "USAGE: ./makeSunilTreeFlat [DATASET] [SELECTIONTYPE=\"ZJet\" or \"DiJet\"] [DIRECTORY=\"/afs/cern.ch/work/s/sunil/public/forTom/\"]" << std::endl;
    exit(11);
  }


  std::string dataset;
  if( argc>1 ) {
    std::string dataset_str(argv[1]);
    dataset = dataset_str;
  }
  TString dataset_tstr(dataset);


  std::string selectionType;
  if( argc>2 ) {
    std::string selectionType_str(argv[2]);
    selectionType = selectionType_str;
  }

  if( selectionType!="ZJet" && selectionType!="DiJet" ) {
    std::cout << "Invalid selectionType: Only 'ZJet' and 'DiJet' currently supported." << std::endl;
    exit(113);
  }


  std::string dir="/afs/cern.ch/work/s/sunil/public/forTom/";
  if( argc>3 ) {
    std::string dir_str(argv[3]);
    dir = dir_str;
  }
    

  std::cout << "dataset: " << dataset << std::endl;
  std::cout << "selectionType: " << selectionType << std::endl;
  std::cout << "dir: " << dir << std::endl;

  std::string sunilTreeName = dir + "/analysis_" + dataset + ".root/Hbb/events";

  TChain* sunilTree = new TChain("events");
  sunilTree->Add(sunilTreeName.c_str());
  
  std::cout << "-> Opening: " <<  sunilTreeName << std::endl;
  int nentries = sunilTree->GetEntries();
  std::cout << "-> Tree has: " << nentries << std::endl;


  std::string flatFileName = "sunilFlat_" + selectionType + "_" + dataset + ".root";
  TFile* newFile = TFile::Open(flatFileName.c_str(), "recreate");
  newFile->cd();

  TTree* flatTree = new TTree("tree_passedEvents", "");

  int run;
  int event;
  float rho;
  float rhoMLP;
  int nvtx;
  std::vector<bool> *triggerResult = 0;
  float jetPt[20];
  float jetEta[20];
  float jetPhi[20];
  float jetEnergy[20];
  float jetBeta[20];
  float jetBtag[20];
  int nCharged[20];
  int nNeutral[20];
  float jetPtD[20];
  float jetPtD_QC[20];
  Float_t jetAxis_QC[2][4];
  Float_t jetAxis[2][4];
  float jetPull_QC[4];
  float jetLead[4];
  int nCharged_QC[20];
  int nCharged_ptCut[20];
  int nNeutral_ptCut[20];
  std::vector<int> *partonId=0;
  std::vector<int> *partonSt=0;
  std::vector<float> *partonPt=0;
  std::vector<float> *partonEta=0;
  std::vector<float> *partonPhi=0;
  std::vector<float> *partonE=0;
  float MuEta[2];
  float MuPhi[2];
  float MuPt[2];
  float MuEnergy[2];


  sunilTree->SetBranchAddress("runNo", &run);
  sunilTree->SetBranchAddress("evtNo", &event);
  sunilTree->SetBranchAddress("nvtx", &nvtx);
  sunilTree->SetBranchAddress("rhoIso", &rho);
  sunilTree->SetBranchAddress("rho", &rhoMLP);
  sunilTree->SetBranchAddress("triggerResult", &triggerResult);

  sunilTree->SetBranchAddress("jetPt", jetPt);
  sunilTree->SetBranchAddress("jetEta", jetEta);
  sunilTree->SetBranchAddress("jetPhi", jetPhi);
  sunilTree->SetBranchAddress("jetEnergy", jetEnergy);
  sunilTree->SetBranchAddress("jetBeta", jetBeta);
  sunilTree->SetBranchAddress("jetBtag", jetBtag);
  sunilTree->SetBranchAddress("jetChgPart", nCharged);
  sunilTree->SetBranchAddress("jetNeutralPart", nNeutral);
  sunilTree->SetBranchAddress("jetPtD", jetPtD);
  sunilTree->SetBranchAddress("jetPtD_QC", jetPtD_QC);
  sunilTree->SetBranchAddress("jetAxis",&jetAxis);
  sunilTree->SetBranchAddress("jetAxis_QC",&jetAxis_QC);
  sunilTree->SetBranchAddress("jetPull_QC",&jetPull_QC);
  sunilTree->SetBranchAddress("jetLead",&jetLead);
  sunilTree->SetBranchAddress("jetChgPart_QC", nCharged_QC);
  sunilTree->SetBranchAddress("jetChgPart_ptcut", nCharged_ptCut);
  sunilTree->SetBranchAddress("jetNeutralPart_ptcut", nNeutral_ptCut);
  sunilTree->SetBranchAddress("partonId",&partonId);
  sunilTree->SetBranchAddress("partonSt",&partonSt);
  sunilTree->SetBranchAddress("partonPt",&partonPt);
  sunilTree->SetBranchAddress("partonEta",&partonEta);
  sunilTree->SetBranchAddress("partonPhi",&partonPhi);
  sunilTree->SetBranchAddress("partonE",&partonE);
  if( selectionType=="ZJet" ) {
    sunilTree->SetBranchAddress("MuPt",MuPt);
    sunilTree->SetBranchAddress("MuEta",MuEta);
    sunilTree->SetBranchAddress("MuPhi",MuPhi);
    sunilTree->SetBranchAddress("MuEnergy",MuEnergy);
  }


  int nJet;
  float eventWeight;
  float wt_pu;
  float wt_pteta;
  float wt_xsec;
  int jetPdgId[20];
  float axis2_QC[20];
  float axis1_QC[20];
  float rmsCand_QC[20];
  float qglMLPJet[20];
  float qglJet[20];
  float qglCorrJet[20];
  float ptZ;

  flatTree->Branch("event", &event, "event/I");
  flatTree->Branch("eventWeight", &eventWeight, "eventWeight/F");
  flatTree->Branch("wt_pu", &wt_pu, "wt_pu/F");
  flatTree->Branch("wt_pteta", &wt_pteta, "wt_pteta/F");
  flatTree->Branch("wt_xsec", &wt_xsec, "wt_xsec/F");
  flatTree->Branch("nvertex", &nvtx, "nvertex/I");
  flatTree->Branch("rhoPF", &rho, "rhoPF/F");
  flatTree->Branch("rhoPF_allJet", &rhoMLP, "rhoPF_allJet/F");
  flatTree->Branch("ptZ", &ptZ, "ptZ/F");
  flatTree->Branch("nJet", &nJet, "nJet/I");
  flatTree->Branch("ptJet", jetPt, "ptJet[nJet]/F");
  flatTree->Branch("etaJet", jetEta, "etaJet[nJet]/F");
  flatTree->Branch("pdgIdJet", jetPdgId, "pdgIdJet[nJet]/I");
  flatTree->Branch("nChargedJet", nCharged, "nChargedJet[nJet]/I");
  flatTree->Branch("nNeutralJet", nNeutral, "nNeutralJet[nJet]/I");
  flatTree->Branch("ptDJet", jetPtD, "ptDJet[nJet]/F");
  flatTree->Branch("ptD_QCJet", jetPtD_QC, "ptDJet_QC[nJet]/F");
  flatTree->Branch("pull_QCJet", jetPull_QC, "pull_QC[nJet]/F");
  flatTree->Branch("RJet", jetLead, "RJet[nJet]/F");
  flatTree->Branch("axis1_QCJet", axis1_QC, "axis1_QCJet[nJet]/F");
  flatTree->Branch("axis2_QCJet", axis2_QC, "axis2_QCJet[nJet]/F");
  flatTree->Branch("rmsCand_QCJet", rmsCand_QC, "rmsCand_QCJet[nJet]/F");
  flatTree->Branch("nChg_QCJet", nCharged_QC, "nChg_QCJet[nJet]/I");
  flatTree->Branch("nNeutral_ptCutJet", nNeutral_ptCut, "nNeutral_ptCutJet[nJet]/I");
  flatTree->Branch("qgMLPJet", qglMLPJet, "qgMLPJet[nJet]/F");
  flatTree->Branch("qglJet", qglJet, "qglJet[nJet]/F");
  flatTree->Branch("qglCorrJet", qglCorrJet, "qglCorrJet[nJet]/F");



  std::cout << "-> Booking QGLikelihoodCalculator..." << std::endl;
  QGLikelihoodCalculator* qglc = new QGLikelihoodCalculator("QuarkGluonTagger/EightTeV/data/");
  std::cout << "-> Booking QGMLPCalculator..." << std::endl;
  QGMLPCalculator* qgmlp = new QGMLPCalculator("MLP","QuarkGluonTagger/EightTeV/data/", true);


  //std::string puFileName = (selectionType=="ZJet") ? "PU_rewt_ZJets.root" : "PU_rewt_flatP6.root";
  //TFile *fPU = TFile::Open(puFileName.c_str());
  //std::string puHistName = (selectionType=="ZJet") ? "hist_puWT" : "hist_WT"; //absurd
  //TH1F *hPU = (TH1F*)fPU->Get(puHistName.c_str());

  //TFile* filePtEtaWeights = TFile::Open("Jetpteta_rewt2D_flatP6.root");
  //TH1F* hPtEta_wt = (TH1F*)filePtEtaWeights->Get("hist_WT");

  std::string rhoWeightFileName = "rhoWeights_" + dataset + ".root";
  TFile* fileRhoWeights = TFile::Open(rhoWeightFileName.c_str());
  TH1F* hPU;
  if( fileRhoWeights!=0 )
    hPU = (TH1F*)fileRhoWeights->Get("rho_weights");


  std::string ptWeightFileName = "ptWeights_" + dataset + ".root";
  TFile* filePtWeights = TFile::Open(ptWeightFileName.c_str());
  TH1F* hPt_wt;
  if( filePtWeights!=0 )
    hPt_wt = (TH1F*)filePtWeights->Get("ptAve_weights");



  std::cout << "-> Begin loop." << std::endl;

  for( unsigned int ientry=0; ientry<nentries; ++ientry ) {

    sunilTree->GetEntry(ientry);

    if( ientry % 100000 == 0 ) std::cout << "entry: " << ientry << "/" << nentries << std::endl;

    bool isMC = (run < 10000);

    // trigger
    if( !isMC ) {
      if(selectionType=="ZJet" ) {
        if(!triggerResult->at(0)) continue;
      } else {
        if( dataset_tstr.Contains("_MBPD_") && !triggerResult->at(19)) continue;
        if( dataset_tstr.Contains("_JetPD_") && !triggerResult->at(1)) continue;
      }
    }


    // event selection:

    if( selectionType=="DiJet" ) {

      if( jetPt[0] < 20.) continue; 
      if( jetPt[1] < 20.) continue;

      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiE( jetPt[0], jetEta[0], jetPhi[0], jetEnergy[0]);
      jet2.SetPtEtaPhiE( jetPt[1], jetEta[1], jetPhi[1], jetEnergy[1]);
      if( fabs(jet1.DeltaPhi(jet2)) < 2.5 ) continue;

      if( jetPt[2] > 0.3*(jetPt[0]+jetPt[1])/2. )continue;

      ptZ = 0.;

    } else if( selectionType=="ZJet" ) {

      if( MuPt[0]<20. ) continue;
      if( MuPt[1]<20. ) continue;

      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiE( MuPt[0], MuEta[0], MuPhi[0], MuEnergy[0] );
      mu2.SetPtEtaPhiE( MuPt[1], MuEta[1], MuPhi[1], MuEnergy[1] );

      TLorentzVector Zmm = mu1 + mu2;

      if( Zmm.M()<70. || Zmm.M()>110. ) continue;
      if( jetPt[0]<20. ) continue; 
      TLorentzVector jet;
      jet.SetPtEtaPhiE( jetPt[0], jetEta[0], jetPhi[0], jetEnergy[0] );
      if( fabs(Zmm.DeltaPhi(jet)) < 2.5 ) continue;
      if( jetPt[1]>0.3*Zmm.Pt() ) continue;

      ptZ = Zmm.Pt();

    }

    // common jet ID:
    if( fabs(jetEta[0]) < 2.5 && jetBeta[0] < ( 1.0 -  0.2 * TMath::Log( nvtx - 0.67))) continue; 
    if( jetBtag[0]>0.244 ) continue;
    if( jetPtD_QC[0]>0.9 ) continue;


    // set event weights:
    wt_pu = 1.;
    wt_pteta = 1.;
    wt_xsec = 1.;

    if( isMC ) {

      // kinematic reweighting only for dijets:
      if( selectionType=="DiJet" ) {
        int ptAveBin = hPt_wt->FindBin( 0.5*(jetPt[0]+jetPt[1]) );
        wt_pteta = hPt_wt->GetBinContent(ptAveBin);

        //int ptetabin = hPtEta_wt->FindBin(jetPt[0],fabs(jetEta[0]));
        //wt_pteta = hPtEta_wt->GetBinContent(ptetabin);
      }

      wt_xsec =  get_xSecWeight(dataset);

      // pu reweighting:
      int bin = hPU->FindBin(rho);
      wt_pu = hPU->GetBinContent(bin);

    }


    eventWeight = wt_xsec*wt_pu*wt_pteta;




    // and now set variables:

    // **** FIRST FOR FIRST JET

    axis1_QC[0] = jetAxis_QC[0][0];
    axis2_QC[0] = jetAxis_QC[1][0];
    rmsCand_QC[0] = (axis1_QC[0]>0. && axis2_QC[0]>0.) ? sqrt( axis1_QC[0]*axis1_QC[0] + axis2_QC[0]*axis2_QC[0] ) : -1.;


    nJet=0;
    for( unsigned ijet=0; ijet<4; ++ijet )
      if( jetPt[ijet]>20. ) nJet++;


    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE( jetPt[0], jetEta[0], jetPhi[0], jetEnergy[0] );

    // match to parton:
    float deltaR_min = 999.;
    int foundPart = -1;

    float deltaR_min_charm = 999.;
    int foundPart_charm = -1;

    float deltaR_min_bottom = 999.;
    int foundPart_bottom = -1;


    for(int iPart=0;iPart<int(partonPt->size());iPart++) {

      if( partonSt->at(iPart) != 3 ) continue;
      if( !( fabs(partonId->at(iPart))<6 || fabs(partonId->at(iPart))>0 || partonId->at(iPart)==21) ) continue;
      TLorentzVector thisPart;
      thisPart.SetPtEtaPhiE( partonPt->at(iPart), partonEta->at(iPart), partonPhi->at(iPart), partonE->at(iPart) );
      float deltaR_part = thisPart.DeltaR(thisJet);
      if(deltaR_part< deltaR_min) {
        deltaR_min = deltaR_part;
        foundPart = iPart;
      }
      if( fabs(partonId->at(iPart))==4 && deltaR_part< deltaR_min_charm) {
        deltaR_min_charm = deltaR_part;
        foundPart_charm = iPart;
      }
      if( fabs(partonId->at(iPart))==5 && deltaR_part< deltaR_min_bottom) {
        deltaR_min_bottom = deltaR_part;
        foundPart_bottom = iPart;
      }

    }


    if( deltaR_min_charm<0.3 && foundPart_charm>=0 ) { // priority to charm
      jetPdgId[0] = partonId->at(foundPart_charm);
    } else if( deltaR_min_bottom<0.3 && foundPart_bottom>=0 ) { // then to bottom
      jetPdgId[0] = partonId->at(foundPart_bottom);
    } else if(deltaR_min < 0.3 && foundPart>=0) {
      jetPdgId[0] = partonId->at(foundPart);
    } else {
      jetPdgId[0] = 0;
    }


    if( !isMC && fabs(jetEta[0])>2.5 ) {
      qglJet[0] = qglc->computeQGLikelihood2012( jetPt[0], jetEta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0]-1, jetPtD_QC[0], axis2_QC[0]);
      qglCorrJet[0] = qglJet[0];
    } else {
      qglJet[0] = qglc->computeQGLikelihood2012( jetPt[0], jetEta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], jetPtD_QC[0], axis2_QC[0]);
      qglCorrJet[0] = (fabs(jetEta[0])<1.9) ? qglJet[0] : qglc->computeQGLikelihood2012( jetPt[0], jetEta[0], rho, correctNCharged(nCharged_QC[0],jetEta[0])+nNeutral_ptCut[0], jetPtD_QC[0], axis2_QC[0]);
    }

    std::map<TString,float> variables_MLP;
    variables_MLP["axis1"]=axis1_QC[0];
    variables_MLP["axis2"]=axis2_QC[0];
    variables_MLP["ptD"]=jetPtD_QC[0];
    if( fabs(jetEta[0])<2.5 ) {
      variables_MLP["mult"]=nCharged_QC[0];
    } else {
      if( isMC ) {
        variables_MLP["mult"]=nCharged_ptCut[0]+nNeutral_ptCut[0];
      } else {
        variables_MLP["mult"]=nCharged_ptCut[0]+nNeutral_ptCut[0]-1;
      }
    }
    
    variables_MLP["pt"]=jetPt[0];
    variables_MLP["eta"]=jetEta[0];
    variables_MLP["rho"]=rhoMLP;

    qglMLPJet[0] = qgmlp->QGvalue(variables_MLP);






    if( selectionType=="DiJet" ) {

      // **** THEN FOR SECOND JET

      axis1_QC[1] = jetAxis_QC[0][1];
      axis2_QC[1] = jetAxis_QC[1][1];
      rmsCand_QC[1] = (axis1_QC[1]>0. && axis2_QC[1]>0.) ? sqrt( axis1_QC[1]*axis1_QC[1] + axis2_QC[1]*axis2_QC[1] ) : -1.;


      TLorentzVector secondJet;
      secondJet.SetPtEtaPhiE( jetPt[1], jetEta[1], jetPhi[1], jetEnergy[1] );

      // match to parton:
      float deltaR_min = 999.;
      int foundPart = -1;

      for(int iPart=0;iPart<int(partonPt->size());iPart++) {

        if( partonSt->at(iPart) != 3 ) continue;
        if( !( fabs(partonId->at(iPart))<6 || fabs(partonId->at(iPart))>0 || partonId->at(iPart)==21) ) continue;
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( partonPt->at(iPart), partonEta->at(iPart), partonPhi->at(iPart), partonE->at(iPart) );
        float deltaR_part = thisPart.DeltaR(secondJet);
        if(deltaR_part< deltaR_min) {
          deltaR_min = deltaR_part;
          foundPart = iPart;
        }

      }

      if(deltaR_min < 0.3 && foundPart>=0) {
        jetPdgId[1] = partonId->at(foundPart);
      } else {
        jetPdgId[1] = 0;
      }


      if( !isMC && fabs(jetEta[1])>2.5 ) {
        qglJet[1] = qglc->computeQGLikelihood2012( jetPt[1], jetEta[1], rho, nCharged_QC[1]+nNeutral_ptCut[1]-1, jetPtD_QC[1], axis2_QC[1]);
        qglCorrJet[1] = qglJet[1];
      } else {
        qglJet[1] = qglc->computeQGLikelihood2012( jetPt[1], jetEta[1], rho, nCharged_QC[1]+nNeutral_ptCut[1], jetPtD_QC[1], axis2_QC[1]);
        qglCorrJet[1] = (fabs(jetEta[1])<1.9) ? qglJet[1] : qglc->computeQGLikelihood2012( jetPt[1], jetEta[1], rho, correctNCharged(nCharged_QC[1],jetEta[1])+nNeutral_ptCut[1], jetPtD_QC[1], axis2_QC[1]);
      }

      std::map<TString,float> variables_MLP;
      variables_MLP["axis1"]=axis1_QC[1];
      variables_MLP["axis2"]=axis2_QC[1];
      variables_MLP["ptD"]=jetPtD_QC[1];
      if( fabs(jetEta[1])<2.5 ) {
        variables_MLP["mult"]=nCharged_QC[1];
      } else {
        if( isMC ) {
          variables_MLP["mult"]=nCharged_ptCut[1]+nNeutral_ptCut[1];
        } else {
          variables_MLP["mult"]=nCharged_ptCut[1]+nNeutral_ptCut[1]-1;
        }
      }
      
      variables_MLP["pt"]=jetPt[1];
      variables_MLP["eta"]=jetEta[1];
      variables_MLP["rho"]=rhoMLP;

      qglMLPJet[1] = qgmlp->QGvalue(variables_MLP);

    } //second jet in dijets


    flatTree->Fill();

  } // for entries


  newFile->cd();
  flatTree->Write();
  newFile->Close();

  return 0;

}




float get_xSecWeight( const std::string& dataset ) {

  return 1.;

}


float correctNCharged(int nChg_QCJet, float eta) {

  float scale = 1.;

  if( fabs(eta)<2.4 && fabs(eta)>1.9 ) {
    // compute the area of the cone outside of the tracker:
    // using the notation found here: http://upload.wikimedia.org/wikipedia/commons/f/fb/Circularsegment.svg
    float R = 0.5;
    float h = fabs(eta)+0.5-2.4;
    float d = R-h;
    float theta = 2.*TMath::ACos(d/R);
    float area_tot = R*R*TMath::Pi();
    float area_outside = R*R*(theta - sin(theta))/2.;
    float fraction_area = area_outside/area_tot;
    scale = 1./(1.-fraction_area);
  }

  return scale*(float)nChg_QCJet;

}

