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

class PatElectronBoosterBDTIso : public edm::EDProducer {
    public:
        explicit PatElectronBoosterBDTIso(const edm::ParameterSet&);
        ~PatElectronBoosterBDTIso();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        edm::InputTag electronTag_;
        std::string   rho_;
        std::string   chargedOption_, neutralsOption_, outputName_;  
        ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_;
        EGammaMvaEleEstimator* eleMVANonTrig;
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
PatElectronBoosterBDTIso::PatElectronBoosterBDTIso(const edm::ParameterSet& iConfig) :
        electronTag_(iConfig.getParameter<edm::InputTag>("src")),
        rho_(iConfig.getParameter<std::string>("rho")),
        chargedOption_(iConfig.existsAs<std::string>("chargedOption") ? iConfig.getParameter<std::string>("chargedOption") : ""),
        neutralsOption_(iConfig.existsAs<std::string>("neutralsOption") ? iConfig.getParameter<std::string>("neutralsOption") : ""),
        outputName_(iConfig.existsAs<std::string>("outputName") ? iConfig.getParameter<std::string>("outputName") : "bdtisonontrig")
{
  produces<pat::ElectronCollection>();  

  // ---- Here I initialize the BDT for the mva ElectronIso

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/WWAnalysis/AnalysisStep/data/ElectronMVAWeights/";

  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_BarrelPt5To10.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_EndcapPt5To10.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_BarrelPt10ToInf.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_EndcapPt10ToInf.weights.xml");  

  eleMVANonTrig  = new EGammaMvaEleEstimator();
  eleMVANonTrig->initialize("EleIso_BDTG_IsoRings",EGammaMvaEleEstimator::kIsoRings,true,manualCatNonTrigWeigths);
  //eleMVANonTrig->SetPrintMVADebug(kFALSE);

  // ---------

  std::string eaTarget = iConfig.getParameter<std::string>("effectiveAreaTarget");
  if      (eaTarget == "NoCorr") effAreaTarget_ = ElectronEffectiveArea::kEleEANoCorr;
  else if (eaTarget == "Data2011") effAreaTarget_ = ElectronEffectiveArea::kEleEAData2011;
  else if (eaTarget == "Data2012") effAreaTarget_ = ElectronEffectiveArea::kEleEAData2012;
  else if (eaTarget == "Summer11MC") effAreaTarget_ = ElectronEffectiveArea::kEleEASummer11MC;
  else if (eaTarget == "Fall11MC") effAreaTarget_ = ElectronEffectiveArea::kEleEAFall11MC;
  else throw cms::Exception("Configuration") << "Unknown effective area " << eaTarget << std::endl;
}


PatElectronBoosterBDTIso::~PatElectronBoosterBDTIso() {
  delete eleMVANonTrig;
}


void PatElectronBoosterBDTIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByLabel(electronTag_,electrons);

    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    // ----- here is the real loop over the electrons ----
    for(edm::View<pat::Electron>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){    
      pat::Electron clone = *ele;
      
      
      // ------ HERE I ADD THE BDT ELE ID VALUE TO THE ELECTRONS
      double mvaValueNonTrig = eleMVANonTrig->isoMvaValue(clone.pt(),
                          clone.superCluster()->eta(),
                          clone.userFloat(rho_),
                          effAreaTarget_,
                          clone.userFloat("electronPFIsoChHad01"+chargedOption_),
                          clone.userFloat("electronPFIsoChHad02"+chargedOption_) - clone.userFloat("electronPFIsoChHad01"+chargedOption_),
                          clone.userFloat("electronPFIsoChHad03"+chargedOption_) - clone.userFloat("electronPFIsoChHad02"+chargedOption_),
                          clone.userFloat("electronPFIsoChHad04"+chargedOption_) - clone.userFloat("electronPFIsoChHad03"+chargedOption_),
                          clone.userFloat("electronPFIsoChHad05"+chargedOption_) - clone.userFloat("electronPFIsoChHad04"+chargedOption_),
                          clone.userFloat("electronPFIsoPhoton01"+neutralsOption_),
                          clone.userFloat("electronPFIsoPhoton02"+neutralsOption_) - clone.userFloat("electronPFIsoPhoton01"+neutralsOption_),
                          clone.userFloat("electronPFIsoPhoton03"+neutralsOption_) - clone.userFloat("electronPFIsoPhoton02"+neutralsOption_),
                          clone.userFloat("electronPFIsoPhoton04"+neutralsOption_) - clone.userFloat("electronPFIsoPhoton03"+neutralsOption_),
                          clone.userFloat("electronPFIsoPhoton05"+neutralsOption_) - clone.userFloat("electronPFIsoPhoton04"+neutralsOption_),
                          clone.userFloat("electronPFIsoNHad01"+neutralsOption_),
                          clone.userFloat("electronPFIsoNHad02"+neutralsOption_) - clone.userFloat("electronPFIsoNHad01"+neutralsOption_),
                          clone.userFloat("electronPFIsoNHad03"+neutralsOption_) - clone.userFloat("electronPFIsoNHad02"+neutralsOption_),
                          clone.userFloat("electronPFIsoNHad04"+neutralsOption_) - clone.userFloat("electronPFIsoNHad03"+neutralsOption_),
                          clone.userFloat("electronPFIsoNHad05"+neutralsOption_) - clone.userFloat("electronPFIsoNHad04"+neutralsOption_),
                          false);

      clone.addUserFloat(outputName_, mvaValueNonTrig);

      // -----------------------------
      

      pOut->push_back(clone);
	
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatElectronBoosterBDTIso::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatElectronBoosterBDTIso::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronBoosterBDTIso);
