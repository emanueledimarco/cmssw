#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <cmath>
#include <algorithm>

namespace { const double MZ = 91.188; }

class FourLeptonBlinder : public edm::EDFilter {
    public:
        explicit FourLeptonBlinder(const edm::ParameterSet&);
        ~FourLeptonBlinder();
        
        virtual bool filter(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag electrons_;
        edm::InputTag muons_;
        StringCutObjectSelector<reco::Candidate> electronSelection_;
        StringCutObjectSelector<reco::Candidate> muonSelection_;
        bool blindMC_, blindData_;
        uint32_t firstBlindedDataRun_;
};




FourLeptonBlinder::FourLeptonBlinder(const edm::ParameterSet& iConfig) :
    electrons_(iConfig.getParameter<edm::InputTag>("electrons")),
    muons_(iConfig.getParameter<edm::InputTag>("muons")),
    electronSelection_(iConfig.getParameter<std::string>("electronSelection")),
    muonSelection_(iConfig.getParameter<std::string>("muonSelection")),
    blindMC_(iConfig.existsAs<bool>("blindMC") ? iConfig.getParameter<bool>("blindMC") : false),
    blindData_(iConfig.existsAs<bool>("blindData") ? iConfig.getParameter<bool>("blindData") : true),
    firstBlindedDataRun_(iConfig.existsAs<uint32_t>("firstBlindedDataRun") ? iConfig.getParameter<uint32_t>("firstBlindedDataRun") : 190000)
{
}

FourLeptonBlinder::~FourLeptonBlinder() {}


bool FourLeptonBlinder::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace std;

    if (iEvent.isRealData()) {
        if (!blindData_) return true;
        if (iEvent.id().run() < firstBlindedDataRun_) return true; 
    } else {
        if (!blindMC_) return true;
    }

    // Read Candidates from Event
    edm::Handle<edm::View<reco::Candidate> > muons, electrons;
    iEvent.getByLabel(muons_, muons);
    iEvent.getByLabel(electrons_, electrons);

    std::vector<const reco::Candidate*> leptons;
    for (edm::View<reco::Candidate>::const_iterator it = muons->begin(), ed = muons->end(); it != ed; ++it) {
        if (muonSelection_(*it)) leptons.push_back(&*it);
    }
    for (edm::View<reco::Candidate>::const_iterator it = electrons->begin(), ed = electrons->end(); it != ed; ++it) {
        if (electronSelection_(*it)) leptons.push_back(&*it);
    }
    if (leptons.size() < 4) return true;

    for (unsigned int i1 = 0, n = leptons.size(); i1 < n; ++i1) {
        for (unsigned int i2 = 0; i2 < n; ++i2) {
            if (i2 == i1) continue;
            if (leptons[i2]->pdgId() != -leptons[i1]->pdgId()) continue; // req. same flavour and opposite sign
            for (unsigned int i3 = 0; i3 < n; ++i3) {
                if (i3 == i1 || i3 == i2) continue;
                for (unsigned int i4 = 0; i4 < n; ++i4) {
                    if (i4 == i1 || i4 == i2 || i4 == i3) continue;
                    if (leptons[i3]->pdgId() != -leptons[i4]->pdgId()) continue; // req. same flavour and opposite sign
                    double m4l = (leptons[i1]->p4() + leptons[i2]->p4() + leptons[i3]->p4() + leptons[i4]->p4()).mass();
                    if ((110 < m4l && m4l < 140) || m4l > 300) return false;
                }
            }
        }
    }

    return true;
}


DEFINE_FWK_MODULE(FourLeptonBlinder);
