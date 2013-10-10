// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" 
#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include <cmath>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

//
// class decleration
//

using namespace edm;
using namespace reco;
using namespace std;
class MuonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MuonAnalyzer(const edm::ParameterSet&);
      ~MuonAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  
  ParameterSet conf_;
	
	MuonMVAEstimator                *fMuonIDMVA;
  MuonMVAEstimator                *fMuonIsoMVA;
	MuonMVAEstimator                *fMuonIsoRingsRadMVA;
  MuonMVAEstimator                *fMuonIDIsoCombinedMVA;

  
  TMVA::Reader             *myTMVAReader;

  unsigned int ev;
      // ----------member data ---------------------------

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig):
  conf_(iConfig)

{
  fMuonIsoMVA = new MuonMVAEstimator();
  vector<string> muoniso_weightfiles;

  muoniso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");                 
  muoniso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
  fMuonIsoMVA->initialize("MuonIso_BDTG_IsoRings",
                   MuonMVAEstimator::kIsoRings,
                   kTRUE,
                   muoniso_weightfiles);
  fMuonIsoMVA->SetPrintMVADebug(kTRUE);
	
	fMuonIDMVA = new MuonMVAEstimator();
  vector<string> muonid_weightfiles;

  muonid_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Tracker_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Global_V0_BDTG.weights.xml");
  fMuonIDMVA->initialize("MuonID_BDTG",
                   MuonMVAEstimator::kID,
                   kTRUE,
                   muonid_weightfiles);
  fMuonIDMVA->SetPrintMVADebug(kTRUE);


 	fMuonIsoRingsRadMVA = new MuonMVAEstimator();
  vector<string> muonisoRingsRad_weightfiles;

  muonisoRingsRad_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V1_BDTG.weights.xml");
  muonisoRingsRad_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V1_BDTG.weights.xml");
  muonisoRingsRad_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V1_BDTG.weights.xml");
  muonisoRingsRad_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V1_BDTG.weights.xml");                 
  muonisoRingsRad_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V1_BDTG.weights.xml");
  muonisoRingsRad_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V1_BDTG.weights.xml");
  fMuonIsoRingsRadMVA->initialize("MuonIso_BDTG_IsoRingsRad",
                   MuonMVAEstimator::kIsoRingsRadial,
                   kTRUE,
                   muonisoRingsRad_weightfiles);
  fMuonIsoRingsRadMVA->SetPrintMVADebug(kTRUE);

 	fMuonIDIsoCombinedMVA = new MuonMVAEstimator();
  vector<string> muonidiso_weightfiles;
  muonidiso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_barrel_lowpt.weights.xml");
  muonidiso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_endcap_lowpt.weights.xml");
  muonidiso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_barrel_highpt.weights.xml");
  muonidiso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_endcap_highpt.weights.xml");
  muonidiso_weightfiles.push_back("Muon/MuonAnalysisTools/data/MuonIDIsoCombinedMVA_V0_tracker.weights.xml");
  fMuonIDIsoCombinedMVA->initialize("MuonIso_BDTG_IsoRings",
                       MuonMVAEstimator::kIDIsoRingsCombined,
                       kTRUE,
                       muonidiso_weightfiles);
  fMuonIDIsoCombinedMVA->SetPrintMVADebug(kTRUE);


}


MuonAnalyzer::~MuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{																	 
																		 
  edm::Handle<reco::VertexCollection> hVertex;
  iEvent.getByLabel("offlinePrimaryVertices", hVertex);
  const reco::VertexCollection *pvCol = hVertex.product();

  Handle<double> hRho;
  edm::InputTag tag("kt6PFJets","rho");
  iEvent.getByLabel(tag,hRho);
  double Rho = *hRho;

  Handle<reco::PFCandidateCollection> hPfCandProduct;
	iEvent.getByLabel("particleFlow", hPfCandProduct);
  const reco::PFCandidateCollection &inPfCands = *(hPfCandProduct.product());

  InputTag gsfEleLabel(string("gsfElectrons"));
  Handle<GsfElectronCollection> theEGammaCollection;
  iEvent.getByLabel(gsfEleLabel,theEGammaCollection);
  const GsfElectronCollection inElectrons = *(theEGammaCollection.product());

  Handle<reco::MuonCollection> hMuonProduct;
  iEvent.getByLabel("muons", hMuonProduct);  
  const reco::MuonCollection inMuons = *(hMuonProduct.product());  

	reco::MuonCollection IdentifiedMuons; //this collection is empty because we already veto pf candidates of type muon
  reco::GsfElectronCollection IdentifiedElectrons;

  for (reco::GsfElectronCollection::const_iterator iE = inElectrons.begin(); 
       iE != inElectrons.end(); ++iE) {

    double electronTrackZ = 0;
    if (iE->gsfTrack().isNonnull()) {
      electronTrackZ = iE->gsfTrack()->dz(pvCol->at(0).position());
    } else if (iE->closestCtfTrackRef().isNonnull()) {
      electronTrackZ = iE->closestCtfTrackRef()->dz(pvCol->at(0).position());
    }    
    if(fabs(electronTrackZ) > 0.2)  continue;

    
    if(fabs(iE->superCluster()->eta())<1.479) {     
      if(iE->pt() > 20) {
        if(iE->sigmaIetaIeta()       > 0.01)  continue;
        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
        if(iE->hadronicOverEm()       > 0.15)  continue;    
      } else {
        if(iE->sigmaIetaIeta()       > 0.012)  continue;
        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
        if(iE->hadronicOverEm()       > 0.15) continue;    
      } 
    } else {     
      if(iE->pt() > 20) {
        if(iE->sigmaIetaIeta()       > 0.03)  continue;
        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
      } else {
        if(iE->sigmaIetaIeta()       > 0.032)  continue;
        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
      }
    }
    IdentifiedElectrons.push_back(*iE);
  }

  for (reco::MuonCollection::const_iterator iM = inMuons.begin(); 
       iM != inMuons.end(); ++iM) {



  	double isomva = fMuonIsoMVA->mvaValue( *iM, pvCol->at(0), 
                                     inPfCands, Rho, 
                                     MuonEffectiveArea::kMuEAFall11MC, 
                                     IdentifiedElectrons, IdentifiedMuons);

																		 
   	double idmva = fMuonIDMVA->mvaValue_ID( *iM, pvCol->at(0));
		 
  	double isoringsradmva = fMuonIsoRingsRadMVA->mvaValue( *iM, pvCol->at(0), 
                                     inPfCands, Rho, 
                                     MuonEffectiveArea::kMuEAFall11MC, 
                                     IdentifiedElectrons, IdentifiedMuons);



  	double idisomva = fMuonIDIsoCombinedMVA->mvaValue( *iM, pvCol->at(0), 
                                     inPfCands, Rho, 
                                     MuonEffectiveArea::kMuEAFall11MC, 
                                     IdentifiedElectrons, IdentifiedMuons);

  	cout << "isomva: " << isomva << endl;

   	cout << "idmva: " << idmva << endl;

  	cout << "isoringsradmva: " << isoringsradmva << endl;
    
		cout << "idisomva: " << idisomva << endl;

	}
}


void 
MuonAnalyzer::beginJob(const edm::EventSetup&)
{

  ev = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonAnalyzer::endJob() {
  cout << " endJob:: #events " << ev << endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyzer);
