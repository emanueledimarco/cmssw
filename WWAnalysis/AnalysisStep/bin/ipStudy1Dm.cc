#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include <DataFormats/FWLite/interface/ChainEvent.h>
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include <WWAnalysis/DataFormats/interface/SkimEvent.h>


#include <TSystem.h>
#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/VectorUtil.h"


#include <sstream>
#include <vector>
//using namespace reco;
 #include "DataFormats/Common/interface/View.h"

//typedef edm::View<reco::RecoCandidate> coll;
//typedef std::vector<pat::Electron> coll;
typedef std::vector<pat::Muon> coll;

void setInput(std::string inputFolder,std::vector<std::string> samples,
	      std::vector<std::string>& fileNames);

void fillHistos(coll::const_iterator it,int type,
		double ip2Dmax, double ip3Dmax, double ip2DSmax, double ip3DSmax, double vtxMax, 
		TH1* h2d, TH1* h3d,TH1* h2ds, TH1* h3ds, TH1* vtx);


//int histoProducer(int sample=0) {
int main(int argc,char* argv[]) {
  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 


  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  // Get the python configuration
  PythonProcessDesc builder(argv[1]);
  edm::ParameterSet const& fwliteParameters = 
    builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Params");
  
  // Now get each parameter
  using namespace std;  
  string inputFolder( fwliteParameters.getParameter<string>("inputFolder") );
  vector<string> sampleStrings( fwliteParameters.getParameter<vector<string> >("input") );
  int maxEvents(  fwliteParameters.getParameter<int>("maxEvents") );

  int typeFill(  fwliteParameters.getParameter<int>("typeFill") );

  string lep1Label(fwliteParameters.getParameter<string>("lep1Label") );
  string lep2Label(fwliteParameters.getParameter<string>("lep2Label") );
  string lep3Label(fwliteParameters.getParameter<string>("lep3Label") );
  string lep4Label(fwliteParameters.getParameter<string>("lep4Label") );

  //int assType(  fwliteParameters.getParameter<int>("assType") );
  //double ptCut(  fwliteParameters.getParameter<double>("ptCut") );

  string fileOutName(  fwliteParameters.getParameter<string>("fileOutName") );


  using namespace ROOT::Math::VectorUtil;

  //TCanvas* myC = new TCanvas("myCanvas","myCanvas",600,600); myC->cd(1)
  TH1D *ip2D_1,*ip3D_1,*ip2DS_1,*ip3DS_1,*vtx_1;
  TH1D *ip2D_2,*ip3D_2,*ip2DS_2,*ip3DS_2,*vtx_2;
  TH1D *ip2D_3,*ip3D_3,*ip2DS_3,*ip3DS_3,*vtx_3;
  TH1D *ip2D_4,*ip3D_4,*ip2DS_4,*ip3DS_4,*vtx_4;


  int ip2Dnbins(200);  double ip2Dmin(0.),ip2Dmax(1000.);
  int ip3Dnbins(200);  double ip3Dmin(0.),ip3Dmax(1000.);
  int ip2DSnbins(100);  double ip2DSmin(0.),ip2DSmax(10.);
  int ip3DSnbins(100);  double ip3DSmin(0.),ip3DSmax(10.);
  int vtxNbins(100);  double vtxMin(0.),vtxMax(40.);




  ip2D_1 = new TH1D("ip2D_1","ip2D_1",ip2Dnbins,ip2Dmin,ip2Dmax);
  ip2D_2 = new TH1D("ip2D_2","ip2D_2",ip2Dnbins,ip2Dmin,ip2Dmax);
  ip2D_3 = new TH1D("ip2D_3","ip2D_3",ip2Dnbins,ip2Dmin,ip2Dmax);
  ip2D_4 = new TH1D("ip2D_4","ip2D_4",ip2Dnbins,ip2Dmin,ip2Dmax);

  ip3D_1 = new TH1D("ip3D_1","ip3D_1",ip3Dnbins,ip3Dmin,ip3Dmax);
  ip3D_2 = new TH1D("ip3D_2","ip3D_2",ip3Dnbins,ip3Dmin,ip3Dmax);
  ip3D_3 = new TH1D("ip3D_3","ip3D_3",ip3Dnbins,ip3Dmin,ip3Dmax);
  ip3D_4 = new TH1D("ip3D_4","ip3D_4",ip3Dnbins,ip3Dmin,ip3Dmax);

  ip2DS_1 = new TH1D("ip2DS_1","ip2DS_1",ip2DSnbins,ip2DSmin,ip2DSmax);
  ip2DS_2 = new TH1D("ip2DS_2","ip2DS_2",ip2DSnbins,ip2DSmin,ip2DSmax);
  ip2DS_3 = new TH1D("ip2DS_3","ip2DS_3",ip2DSnbins,ip2DSmin,ip2DSmax);
  ip2DS_4 = new TH1D("ip2DS_4","ip2DS_4",ip2DSnbins,ip2DSmin,ip2DSmax);

  ip3DS_1 = new TH1D("ip3DS_1","ip3DS_1",ip3DSnbins,ip3DSmin,ip3DSmax);
  ip3DS_2 = new TH1D("ip3DS_2","ip3DS_2",ip3DSnbins,ip3DSmin,ip3DSmax);
  ip3DS_3 = new TH1D("ip3DS_3","ip3DS_3",ip3DSnbins,ip3DSmin,ip3DSmax);
  ip3DS_4 = new TH1D("ip3DS_4","ip3DS_4",ip3DSnbins,ip3DSmin,ip3DSmax);

  vtx_1 = new TH1D("vtx_1","vtx_1",vtxNbins,vtxMin,vtxMax);
  vtx_2 = new TH1D("vtx_2","vtx_2",vtxNbins,vtxMin,vtxMax);
  vtx_3 = new TH1D("vtx_3","vtx_3",vtxNbins,vtxMin,vtxMax);
  vtx_4 = new TH1D("vtx_4","vtx_4",vtxNbins,vtxMin,vtxMax);


  vector<string> fileNames;
  setInput(inputFolder,sampleStrings,fileNames);
  fwlite::ChainEvent ev(fileNames);

  //######### End of initialization ###############
  //###############################################


  //########## Main event LOOP ####################
  int counter(0);
  for( ev.toBegin(); ! ev.atEnd(); ++ev) {
    counter++; 
    if(counter == maxEvents) break;
    //if(counter == 20) break;
    if(counter % 1000 == 1) {   
      cout << "counter: " << counter << endl;
      //cout << "file name: " << ev.getTFile()->GetName() << endl;
    }
    

    fwlite::Handle<coll > leptons1;
    leptons1.getByLabel(ev,lep1Label.c_str());

    fwlite::Handle<coll > leptons2;
    leptons2.getByLabel(ev,lep2Label.c_str());

    fwlite::Handle<coll > leptons3;
    leptons3.getByLabel(ev,lep3Label.c_str());

    fwlite::Handle<coll > leptons4;
    leptons4.getByLabel(ev,lep4Label.c_str());



    //if(skimEventH->size()==0) continue;    
    if(! leptons1.isValid()) continue;
    for(coll::const_iterator it1 = leptons1.ptr()->begin();
	it1 != leptons1.ptr()->end(); it1++){
	if(it1->pt()>10.0)
	  fillHistos(it1,typeFill,ip2Dmax,ip3Dmax,ip2DSmax,ip3DSmax,vtxMax,
		     ip2D_1,ip3D_1,ip2DS_1,ip3DS_1,vtx_1);
    }

    if(! leptons2.isValid()) continue;
    for(coll::const_iterator it2 = leptons2.ptr()->begin();
	it2 != leptons2.ptr()->end(); it2++){
      fillHistos(it2,typeFill,ip2Dmax,ip3Dmax,ip2DSmax,ip3DSmax,vtxMax,
		 ip2D_2,ip3D_2,ip2DS_2,ip3DS_2,vtx_2);
    }

    if(! leptons3.isValid()) continue;
    for(coll::const_iterator it3 = leptons3.ptr()->begin();
	it3 != leptons3.ptr()->end(); it3++){
      fillHistos(it3,typeFill,ip2Dmax,ip3Dmax,ip2DSmax,ip3DSmax,vtxMax,
		 ip2D_3,ip3D_3,ip2DS_3,ip3DS_3,vtx_3);
    }

    /*
    for(coll::const_iterator it4 = leptons4.ptr()->begin();
	it4 != leptons4.ptr()->end(); it4++){
      fillHistos(it4,typeFill,ip2Dmax,ip3Dmax,ip2DSmax,ip3DSmax,vtxMax,
		 ip2D_4,ip3D_4,ip2DS_4,ip3DS_4,vtx_4);
    }
    */


    
  }// end loop over edm::events 
  //###############################################


  // write histograms on disk
  TFile output(fileOutName.c_str(),"RECREATE");  
  ip2D_1->Write();
  ip2D_2->Write();
  ip2D_3->Write();
  ip2D_4->Write();

  ip3D_1->Write();
  ip3D_2->Write();
  ip3D_3->Write();
  ip3D_4->Write();

  ip2DS_1->Write();
  ip2DS_2->Write();
  ip2DS_3->Write();
  ip2DS_4->Write();

  ip3DS_1->Write();
  ip3DS_2->Write();
  ip3DS_3->Write();
  ip3DS_4->Write();

  vtx_1->Write();
  vtx_2->Write();
  vtx_3->Write();
  vtx_4->Write();

  //output.Write();   

  return 0;
}


 

void setInput(std::string inputFolder,std::vector<std::string> samples,
	      std::vector<std::string>& fileNames){  
  using namespace std;
  for(unsigned int i=0;i<samples.size(); ++i){ 
    std::string tmpString=inputFolder+samples[i];
    //tmpString+=".root";
    cout << "input: " << tmpString << endl;    
    fileNames.push_back(tmpString.c_str());
  }
}



void fillHistos(coll::const_iterator it,int type,
		double ip2Dmax, double ip3Dmax, double ip2DSmax, double ip3DSmax, double vtxMax, 
		TH1* h2d, TH1* h3d,TH1* h2ds, TH1* h3ds, TH1* vtx){
	       
  using namespace std;


  if(type==1){
    h2d->Fill(std::min((double) it->userFloat("tip")*10000.,ip2Dmax-0.00001));
    h3d->Fill(std::min((double) it->userFloat("ip")*10000.,ip3Dmax-0.00001));

    h2ds->Fill(std::min((double) it->userFloat("tip")/it->userFloat("tipErr"),ip2DSmax-0.00001));
    h3ds->Fill(std::min((double) it->userFloat("ip")/it->userFloat("ipErr"),ip3DSmax-0.00001));

    vtx->Fill(std::min((double) it->userFloat("vtxComp"),vtxMax-0.00001));
  }

  if(type==2){
    h2d->Fill(std::min((double) it->userFloat("tip2")*10000.,ip2Dmax-0.00001));
    h3d->Fill(std::min((double) it->userFloat("ip2")*10000.,ip3Dmax-0.00001));

    h2ds->Fill(std::min((double) it->userFloat("tip2")/it->userFloat("tipErr2"),ip2DSmax-0.00001));
    h3ds->Fill(std::min((double) it->userFloat("ip2")/it->userFloat("ipErr2"),ip3DSmax-0.00001));

    vtx->Fill(std::min((double) it->userFloat("vtxComp2"),vtxMax-0.00001));
  }

}		 
