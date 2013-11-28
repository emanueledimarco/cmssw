// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//BDT MuonID
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"


#include <vector>
#include <string>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatMuonBoosterBDTID : public edm::EDProducer {
    public:
        explicit PatMuonBoosterBDTID(const edm::ParameterSet&);
        ~PatMuonBoosterBDTID();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        edm::InputTag muonTag_;
        edm::InputTag vertexsTag_;
        edm::InputTag pfCandsTag_;
        edm::InputTag rhoTag_;
        double dzCut_;
        std::string outputName_;
        MuonMVAEstimator* muMVANonTrig;
        std::vector<std::string> manualCatNonTrigWeigths;
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
PatMuonBoosterBDTID::PatMuonBoosterBDTID(const edm::ParameterSet& iConfig) :
  muonTag_(iConfig.getParameter<edm::InputTag>("src")),
  vertexsTag_(iConfig.getParameter<edm::InputTag>("vertexs")),
  pfCandsTag_(iConfig.getParameter<edm::InputTag>("pfCands")),
  rhoTag_(iConfig.getParameter<edm::InputTag>("rho")),
  dzCut_(iConfig.getParameter<double>("dzCut")),
  outputName_(iConfig.getParameter<std::string>("outputName"))
{
  produces<pat::MuonCollection>();  

  // ---- Here I initialize the BDT for the mva MuonID

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/Muon/MuonAnalysisTools/data/";

  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-Tracker_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-Global_V0_BDTG.weights.xml");

  muMVANonTrig  = new MuonMVAEstimator();
  muMVANonTrig->initialize("MuonID_BDTG",MuonMVAEstimator::kID,true,manualCatNonTrigWeigths);
  muMVANonTrig->SetPrintMVADebug(kFALSE);
  //muMVANonTrig->SetPrintMVADebug(kTRUE);

  // ---------
}


PatMuonBoosterBDTID::~PatMuonBoosterBDTID() {
  delete muMVANonTrig;
}


void PatMuonBoosterBDTID::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<reco::Candidate> > muons;
    iEvent.getByLabel(muonTag_,muons);

    edm::Handle<reco::VertexCollection> vertexs;
    iEvent.getByLabel(vertexsTag_,vertexs);

    Handle<reco::PFCandidateCollection> pfCands;
    iEvent.getByLabel(pfCandsTag_,pfCands);

    edm::Handle<double> hRho;
    iEvent.getByLabel(rhoTag_,hRho);


    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);

    // ----- here is the real loop over the muons ----
    for(edm::View<reco::Candidate>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){    
      const pat::MuonRef musRef = edm::RefToBase<reco::Candidate>(muons,mu-muons->begin()).castTo<pat::MuonRef>();
      pat::Muon clone = *edm::RefToBase<reco::Candidate>(muons,mu-muons->begin()).castTo<pat::MuonRef>();
      
      const reco::GsfElectronCollection dummyIdentifiedEleCollection;
      const reco::MuonCollection dummyIdentifiedMuCollection;

      double mvaValueNonTrig = muMVANonTrig->mvaValue(clone,
						      vertexs->front(),
						      *pfCands,
						      *hRho,MuonEffectiveArea::kMuEAFall11MC,
						      dummyIdentifiedEleCollection,
						      dummyIdentifiedMuCollection,
						      dzCut_);     
     
      clone.addUserFloat(outputName_,mvaValueNonTrig);

      // -----------------------------
      

      pOut->push_back(clone);
	
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatMuonBoosterBDTID::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatMuonBoosterBDTID::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(PatMuonBoosterBDTID);
