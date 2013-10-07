// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"

#include <vector>
#include <string>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatElectronEffAreaIso : public edm::EDProducer {
    public:
        explicit PatElectronEffAreaIso(const edm::ParameterSet&);
        ~PatElectronEffAreaIso() {}
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
        // ----------member data ---------------------------
        edm::InputTag electronTag_;
        ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_;
        ElectronEffectiveArea::ElectronEffectiveAreaType   effAreaGamma_, effAreaNeutralHad_, effAreaTotal_;
        std::string rho_; 
        std::string deltaR_, chargedOption_, neutralsOption_; 
        std::string label_; 
        enum Truncation { Both, Sum, None } truncate_;
        bool separateEAs_;
};

PatElectronEffAreaIso::PatElectronEffAreaIso(const edm::ParameterSet& iConfig) :
    electronTag_(iConfig.getParameter<edm::InputTag>("src")),
    rho_(iConfig.getParameter<std::string>("rho")),
    deltaR_(iConfig.getParameter<std::string>("deltaR")),
    chargedOption_(iConfig.existsAs<std::string>("chargedOption") ? iConfig.getParameter<std::string>("chargedOption") : ""),
    neutralsOption_(iConfig.existsAs<std::string>("neutralsOption") ? iConfig.getParameter<std::string>("neutralsOption") : ""),
    label_(iConfig.getParameter<std::string>("label")),
    separateEAs_(iConfig.getParameter<bool>("separatePhotonAndHadronEAs"))
{
  std::string truncate = iConfig.getParameter<std::string>("truncateAtZero");
  if (separateEAs_) {
      if (truncate == "both") truncate_ = Both;
      else if (truncate == "sum") truncate_ = Sum;
      else if (truncate == "no")  truncate_ = None;
      else throw cms::Exception("Configuration") << "With separatePhotonAndHadronEAs=true, truncateAtZero must be 'both', 'sum' or 'no'\n";
  } else {
      if (truncate == "yes") truncate_ = Sum;
      else if (truncate == "no")  truncate_ = None;
      else throw cms::Exception("Configuration") << "With separatePhotonAndHadronEAs=false, truncateAtZero must be 'yes' or 'no'\n";
  }
  produces<pat::ElectronCollection>();  
  std::string eaTarget = iConfig.getParameter<std::string>("effectiveAreaTarget");
  if      (eaTarget == "NoCorr") effAreaTarget_ = ElectronEffectiveArea::kEleEANoCorr;
  else if (eaTarget == "Data2011") effAreaTarget_ = ElectronEffectiveArea::kEleEAData2011;
  else if (eaTarget == "Data2012") effAreaTarget_ = ElectronEffectiveArea::kEleEAData2012;
  else if (eaTarget == "Summer11MC") effAreaTarget_ = ElectronEffectiveArea::kEleEASummer11MC;
  else if (eaTarget == "Fall11MC") effAreaTarget_ = ElectronEffectiveArea::kEleEAFall11MC;
  else throw cms::Exception("Configuration") << "Unknown effective area " << eaTarget << "\n";
  if (deltaR_ == "03") {
      effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso03;
      effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso03;
      effAreaTotal_ = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;
  } else if (deltaR_ == "04") {
      effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso04;
      effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso04;
      effAreaTotal_ = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04;
  } else throw cms::Exception("Configuration") << "Unsupported deltaR " << deltaR_ << "\n";
}


void PatElectronEffAreaIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByLabel(electronTag_,electrons);

    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    // ----- here is the real loop over the electrons ----
    for(edm::View<pat::Electron>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){    
      pat::Electron clone = *ele;
      float rho = clone.userFloat(rho_), abseta = abs(clone.superCluster()->eta());
      float eff_area_ga = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
      float eff_area_nh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
      float eff_area_tot = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaTotal_, abseta, effAreaTarget_);
      clone.addUserFloat(label_+"EApho",    eff_area_ga);
      clone.addUserFloat(label_+"EAneuHad", eff_area_nh);
      clone.addUserFloat(label_+"EAtot",    eff_area_tot);
      clone.addUserFloat(label_+"Rho",      rho);
      float iso = 0;
      float nhiso = clone.userFloat("electronPFIsoNHad"+deltaR_+neutralsOption_), phiso = clone.userFloat("electronPFIsoPhoton"+deltaR_+neutralsOption_);
      clone.addUserFloat(label_+"NHad", nhiso);
      clone.addUserFloat(label_+"Pho",  phiso);
      if (separateEAs_) {
        if (truncate_ == Both) {
            iso += max<float>(0.f, nhiso - eff_area_nh*rho);
            iso += max<float>(0.f, phiso - eff_area_ga*rho);
        } else {
            iso += nhiso + phiso - (eff_area_nh+eff_area_ga)*rho;
        }
      } else {
        iso += nhiso + phiso - eff_area_tot*rho;
      }
      if (truncate_ == Sum && iso < 0) iso = 0;
      float chiso = clone.userFloat("electronPFIsoChHad"+deltaR_+chargedOption_);
      float chpu  = clone.userFloat("electronPFIsoChHadPU"+deltaR_);
      iso += chiso;
      clone.addUserFloat(label_+"ChPU", chpu);
      clone.addUserFloat(label_+"ChHad", chiso);
      clone.addUserFloat(label_, iso);
      pOut->push_back(clone);
    }
    iEvent.put(pOut);
}

class PatMuonEffAreaIso : public edm::EDProducer {
    public:
        explicit PatMuonEffAreaIso(const edm::ParameterSet&);
        ~PatMuonEffAreaIso() {}
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
        // ----------member data ---------------------------
        edm::InputTag muonTag_;
        MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_;
        MuonEffectiveArea::MuonEffectiveAreaType   effAreaGamma_, effAreaNeutralHad_, effAreaTotal_;
        std::string rho_; 
        std::string deltaR_, chargedOption_, neutralsOption_; 
        std::string label_; 
        enum Truncation { Both, Sum, None } truncate_;
        bool separateEAs_;
};

PatMuonEffAreaIso::PatMuonEffAreaIso(const edm::ParameterSet& iConfig) :
    muonTag_(iConfig.getParameter<edm::InputTag>("src")),
    rho_(iConfig.getParameter<std::string>("rho")),
    deltaR_(iConfig.getParameter<std::string>("deltaR")),
    chargedOption_(iConfig.existsAs<std::string>("chargedOption") ? iConfig.getParameter<std::string>("chargedOption") : ""),
    neutralsOption_(iConfig.existsAs<std::string>("neutralsOption") ? iConfig.getParameter<std::string>("neutralsOption") : ""),
    label_(iConfig.getParameter<std::string>("label")),
    separateEAs_(iConfig.getParameter<bool>("separatePhotonAndHadronEAs"))
{
  std::string truncate = iConfig.getParameter<std::string>("truncateAtZero");
  if (separateEAs_) {
      if (truncate == "both") truncate_ = Both;
      else if (truncate == "sum") truncate_ = Sum;
      else if (truncate == "no")  truncate_ = None;
      else throw cms::Exception("Configuration") << "With separatePhotonAndHadronEAs=true, truncateAtZero must be 'both', 'sum' or 'no'\n";
  } else {
      if (truncate == "yes") truncate_ = Sum;
      else if (truncate == "no")  truncate_ = None;
      else throw cms::Exception("Configuration") << "With separatePhotonAndHadronEAs=false, truncateAtZero must be 'yes' or 'no'\n";
  }
  produces<pat::MuonCollection>();  
  std::string eaTarget = iConfig.getParameter<std::string>("effectiveAreaTarget");
  if      (eaTarget == "NoCorr") effAreaTarget_ = MuonEffectiveArea::kMuEANoCorr;
  else if (eaTarget == "Data2011") effAreaTarget_ = MuonEffectiveArea::kMuEAData2011;
  else if (eaTarget == "Data2012") effAreaTarget_ = MuonEffectiveArea::kMuEAData2012;
  else if (eaTarget == "Summer11MC") effAreaTarget_ = MuonEffectiveArea::kMuEASummer11MC;
  else if (eaTarget == "Fall11MC") effAreaTarget_ = MuonEffectiveArea::kMuEAFall11MC;
  else throw cms::Exception("Configuration") << "Unknown effective area " << eaTarget << "\n";
  if (deltaR_ == "03") {
      effAreaGamma_   = MuonEffectiveArea::kMuGammaIso03;
      effAreaNeutralHad_ = MuonEffectiveArea::kMuNeutralHadronIso03;
      effAreaTotal_ = MuonEffectiveArea::kMuGammaAndNeutralHadronIso03;
  } else if (deltaR_ == "04") {
      effAreaGamma_   = MuonEffectiveArea::kMuGammaIso04;
      effAreaNeutralHad_ = MuonEffectiveArea::kMuNeutralHadronIso04;
      effAreaTotal_ = MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;
  } else throw cms::Exception("Configuration") << "Unsupported deltaR " << deltaR_ << "\n";
}


void PatMuonEffAreaIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByLabel(muonTag_,muons);

    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);

    // ----- here is the real loop over the muons ----
    for(edm::View<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){    
      pat::Muon clone = *mu;
      float rho = clone.userFloat(rho_), abseta = abs(clone.eta());
      float eff_area_ga  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
      float eff_area_nh  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
      float eff_area_tot = MuonEffectiveArea::GetMuonEffectiveArea(effAreaTotal_, abseta, effAreaTarget_);
      clone.addUserFloat(label_+"EApho",    eff_area_ga);
      clone.addUserFloat(label_+"EAneuHad", eff_area_nh);
      clone.addUserFloat(label_+"EAtot",    eff_area_tot);
      clone.addUserFloat(label_+"Rho",      rho);
      float iso = 0;
      float nhiso = clone.userFloat("muonPFIsoNHad"+deltaR_+neutralsOption_), phiso = clone.userFloat("muonPFIsoPhoton"+deltaR_+neutralsOption_);
      clone.addUserFloat(label_+"NHad", nhiso);
      clone.addUserFloat(label_+"Pho",  phiso);
      if (separateEAs_) {
        if (truncate_ == Both) {
            iso += max<float>(0.f, nhiso - eff_area_nh*rho);
            iso += max<float>(0.f, phiso - eff_area_ga*rho);
        } else {
            iso += nhiso + phiso - (eff_area_nh+eff_area_ga)*rho;
        }
      } else {
        iso += nhiso + phiso - eff_area_tot*rho;
      }
      if (truncate_ == Sum && iso < 0) iso = 0;
      float chiso = clone.userFloat("muonPFIsoChHad"+deltaR_+chargedOption_);
      float chpu  = clone.userFloat("muonPFIsoChHadPU"+deltaR_);
      iso += chiso;
      clone.addUserFloat(label_+"ChPU", chpu);
      clone.addUserFloat(label_+"ChHad", chiso);
      clone.addUserFloat(label_, iso);
      pOut->push_back(clone);
    }
    iEvent.put(pOut);
}


//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronEffAreaIso);
DEFINE_FWK_MODULE(PatMuonEffAreaIso);
