// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//BDT ElectronID
//#include "WWAnalysis/AnalysisStep/interface/ElectronIDMVAHZZ.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"


#include <vector>
#include <string>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatElectronBoosterBDTID : public edm::EDProducer {
    public:
        explicit PatElectronBoosterBDTID(const edm::ParameterSet&);
        ~PatElectronBoosterBDTID();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        edm::InputTag electronTag_;
        EGammaMvaEleEstimator* eleMVATrig;
        EGammaMvaEleEstimator* eleMVANonTrig;
        std::vector<std::string> manualCatTrigWeigths;
        std::vector<std::string> manualCatNonTrigWeigths;
        std::string postfix_;
        bool doTrig_, doNonTrig_;
        bool debug_;
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
PatElectronBoosterBDTID::PatElectronBoosterBDTID(const edm::ParameterSet& iConfig) :
        electronTag_(iConfig.getParameter<edm::InputTag>("src")),
        eleMVATrig(0),eleMVANonTrig(0),
        postfix_(iConfig.existsAs<std::string>("postfix") ? iConfig.getParameter<std::string>("postfix") : ""),
        doTrig_(iConfig.existsAs<bool>("doTrigMVA") ? iConfig.getParameter<bool>("doTrigMVA") : true),
        doNonTrig_(iConfig.existsAs<bool>("doNonTrigMVA") ? iConfig.getParameter<bool>("doNonTrigMVA") : true),
        debug_(iConfig.getUntrackedParameter<bool>("verbose",false))
{
  produces<pat::ElectronCollection>();  

  // ---- Here I initialize the BDT for the mva ElectronID

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/WWAnalysis/AnalysisStep/data/ElectronMVAWeights/";
  manualCatTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_TrigV0_Cat1.weights.xml");
  manualCatTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_TrigV0_Cat2.weights.xml");
  manualCatTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_TrigV0_Cat3.weights.xml");
  manualCatTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_TrigV0_Cat4.weights.xml");
  manualCatTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_TrigV0_Cat5.weights.xml");
  manualCatTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_TrigV0_Cat6.weights.xml");  

  manualCatNonTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_NonTrigV0_Cat1.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_NonTrigV0_Cat2.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_NonTrigV0_Cat3.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_NonTrigV0_Cat4.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_NonTrigV0_Cat5.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/Electrons_BDTG_NonTrigV0_Cat6.weights.xml");  

  if (doTrig_) {
      eleMVATrig  = new EGammaMvaEleEstimator();
      eleMVATrig->initialize("BDT",EGammaMvaEleEstimator::kTrig,true,manualCatTrigWeigths);
  }
  if (doNonTrig_) {
      eleMVANonTrig  = new EGammaMvaEleEstimator();
      eleMVANonTrig->initialize("BDT",EGammaMvaEleEstimator::kNonTrig,true,manualCatNonTrigWeigths);
  }
  
  //  eleMVATrig->SetPrintMVADebug(kTRUE);
//   eleMVANonTrig->SetPrintMVADebug(kTRUE);


  // ---------
}


PatElectronBoosterBDTID::~PatElectronBoosterBDTID() {
  delete eleMVATrig;
  delete eleMVANonTrig;
}


void PatElectronBoosterBDTID::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<reco::Candidate> > electrons;
    iEvent.getByLabel(electronTag_,electrons);

    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    // ----- here is the real loop over the electrons ----
    for(edm::View<reco::Candidate>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){    
      const pat::ElectronRef elecsRef = edm::RefToBase<reco::Candidate>(electrons,ele-electrons->begin()).castTo<pat::ElectronRef>();
      pat::Electron clone = *edm::RefToBase<reco::Candidate>(electrons,ele-electrons->begin()).castTo<pat::ElectronRef>();

      // ------ HERE I ADD THE BDT ELE ID VALUE TO THE ELECTRONS
      double xieSign  = ( (-clone.userFloat("dxyPV")) >=0 )  ? 1: -1;
      float ctfChi2 = clone.userFloat("ctfChi2"); if (ctfChi2 == -1) ctfChi2 = 0; // fix a posteriori problem in booster

      if (doTrig_) {
      double mvaValueTrig = eleMVATrig->mvaValue(
                     clone.fbrem(),
                     ctfChi2,
                     clone.userInt("ctfHits"),   
                     clone.gsfTrack()->normalizedChi2(),
					 clone.deltaEtaSuperClusterTrackAtVtx(),
					 clone.deltaPhiSuperClusterTrackAtVtx(),
					 clone.deltaEtaSeedClusterTrackAtCalo(),
					 clone.sigmaIetaIeta(), 
					 clone.userFloat("sigmaIphiIphi"),
                     clone.userFloat("etaWidth"), 
                     clone.userFloat("phiWidth"), 
					 (clone.e5x5()) !=0. ? 1.-(clone.e1x5()/clone.e5x5()) : -1.,
                     clone.userFloat("e3x3")/(clone.superCluster()->rawEnergy()),
					 clone.hcalOverEcal(),
					 clone.eSuperClusterOverP(), 
                     (1.0 / clone.ecalEnergy()) - (1.0 / clone.p()),   
                     clone.eEleClusterOverPout(),
                     clone.superCluster()->preshowerEnergy()/clone.superCluster()->rawEnergy(),          
					 -clone.userFloat("dxyPV"),
					 xieSign*clone.userFloat("ip"),
                     clone.superCluster()->eta(), 
                     clone.pt(), 
                     debug_);	
      clone.addUserFloat(std::string("bdttrig")+postfix_,mvaValueTrig);
      }

      if (doNonTrig_) {
      double mvaValueNonTrig = eleMVANonTrig->mvaValue(
                     clone.fbrem(),
                     ctfChi2,
                     clone.userInt("ctfHits"),   
                     clone.gsfTrack()->normalizedChi2(),
                     clone.deltaEtaSuperClusterTrackAtVtx(),
                     clone.deltaPhiSuperClusterTrackAtVtx(),
                     clone.deltaEtaSeedClusterTrackAtCalo(),
                     clone.sigmaIetaIeta(),
                     clone.userFloat("sigmaIphiIphi"),
                     clone.userFloat("etaWidth"),
                     clone.userFloat("phiWidth"),
					 (clone.e5x5()) !=0. ? 1.-(clone.e1x5()/clone.e5x5()) : -1.,
                     clone.userFloat("e3x3")/(clone.superCluster()->rawEnergy()),
                     clone.hcalOverEcal(),
                     clone.eSuperClusterOverP(), 
                     (1.0 / clone.ecalEnergy()) - (1.0 / clone.p()),
                     clone.eEleClusterOverPout(),
                     clone.superCluster()->preshowerEnergy()/clone.superCluster()->rawEnergy(),
                     clone.superCluster()->eta(),
                     clone.pt(),
                     debug_);
      clone.addUserFloat(std::string("bdtnontrig"+postfix_),mvaValueNonTrig);
      }

      // -----------------------------

      pOut->push_back(clone);
	
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatElectronBoosterBDTID::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatElectronBoosterBDTID::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronBoosterBDTID);
