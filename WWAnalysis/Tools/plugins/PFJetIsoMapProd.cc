#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

class PFJetIsoDepMapProd : public edm::EDProducer {
    public:
        explicit PFJetIsoDepMapProd(const edm::ParameterSet&);
        ~PFJetIsoDepMapProd();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag vtxLabel_;
        edm::InputTag candLabel_;
        edm::InputTag pfLabel_;
        edm::InputTag pfjetLabel_;

        enum ConeMode { JetCone, FixedCone, MaxCone, SumCone };
        ConeMode coneMode_;
        double coneSize_;

        double dzCut_;
};



PFJetIsoDepMapProd::PFJetIsoDepMapProd(const edm::ParameterSet& iConfig) :
    vtxLabel_(iConfig.getParameter<edm::InputTag>("vtxLabel")),
    candLabel_(iConfig.getParameter<edm::InputTag>("candLabel")),
    pfjetLabel_(iConfig.getParameter<edm::InputTag>("pfjetLabel"))
{
    std::string coneMode(iConfig.getParameter<std::string>("coneMode"));
    if      (coneMode == "jet")   coneMode_ = JetCone;
    else if (coneMode == "fixed") coneMode_ = FixedCone;
    else if (coneMode == "max")   coneMode_ = MaxCone;
    else if (coneMode == "sum")   coneMode_ = SumCone;
    else throw cms::Exception("Configuration") << "coneMode must be one of 'jet', 'fixed', 'max', 'sum'.\n";
    if (coneMode_ != JetCone) {
        pfLabel_ =  iConfig.getParameter<edm::InputTag>("pfLabel");
        coneSize_ = iConfig.getParameter<double>("coneSize");
        dzCut_    = iConfig.getParameter<double>("dzCut");
    }
    /* 
    * Mode 'jet':   use jet constituents
    * Mode 'fixed': use fixed cone
    * Mode 'max':   use dynamic cone R=max(dR(jl), <cone-size>) 
    * Mode 'sum':   use dynamic cone R=dR(jl) + <cone-size>
    */

    produces<edm::ValueMap<reco::IsoDeposit> >("charged");
    produces<edm::ValueMap<reco::IsoDeposit> >("neutralHad");
    produces<edm::ValueMap<reco::IsoDeposit> >("photon");
}

void PFJetIsoDepMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::RecoCandidate> > candH;
    iEvent.getByLabel(candLabel_,candH);

    edm::Handle<reco::PFCandidateCollection> pfH;
    if (coneMode_ != JetCone) iEvent.getByLabel(pfLabel_,pfH);

    edm::Handle<reco::PFJetCollection> pfJH;
    iEvent.getByLabel(pfjetLabel_,pfJH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxLabel_,vtxH);

    size_t size = candH->size();
    std::vector<reco::IsoDeposit> chargedV; chargedV.reserve(size);
    std::vector<reco::IsoDeposit> neutralV; neutralV.reserve(size);
    std::vector<reco::IsoDeposit> photonV; photonV.reserve(size);
    std::auto_ptr<edm::ValueMap<reco::IsoDeposit> > chargedM(new edm::ValueMap<reco::IsoDeposit> ());
    std::auto_ptr<edm::ValueMap<reco::IsoDeposit> > neutralM(new edm::ValueMap<reco::IsoDeposit> ());
    std::auto_ptr<edm::ValueMap<reco::IsoDeposit> > photonM(new edm::ValueMap<reco::IsoDeposit> ());
    edm::ValueMap<reco::IsoDeposit>::Filler chargedF(*chargedM);
    edm::ValueMap<reco::IsoDeposit>::Filler neutralF(*neutralM);
    edm::ValueMap<reco::IsoDeposit>::Filler photonF(*photonM);

    for(size_t i=0; i<size;++i) {
        const reco::RecoCandidate &cand = (*candH)[i];
        float leta = cand.eta(), lphi = cand.phi();

        double zLepton = 0.0;
        if(cand.track().isNonnull()) zLepton = cand.track()->dz(vtxH->at(0).position());

        reco::TrackRef    candTrackRef;
        reco::GsfTrackRef candGsfTrackRef;
        const reco::Muon        *mu = dynamic_cast<const reco::Muon *>(&cand);
        const reco::GsfElectron *el = dynamic_cast<const reco::GsfElectron *>(&cand);
        if (mu != 0) {
            // get track, bypassing PAT embedding
            candTrackRef = mu->reco::Muon::innerTrack();
        } else if (el != 0) {
            // get tracks, bypassing PAT embedding
            candTrackRef = el->reco::GsfElectron::closestCtfTrackRef();
            candGsfTrackRef = el->reco::GsfElectron::gsfTrack();
        }

        const reco::PFJet *nearest = 0; double dr2Nearest = 9999;
        for (size_t j=0, nj = pfJH->size(); j < nj; j++) { 
            double thisdr2 = reco::deltaR2(cand, (*pfJH)[j]);
            if (thisdr2 < dr2Nearest) {
                dr2Nearest = thisdr2;
                nearest = & (*pfJH)[j];
            }
        }
        if (nearest == 0) {
            chargedV.push_back(reco::IsoDeposit());
            neutralV.push_back(reco::IsoDeposit());
            photonV.push_back(reco::IsoDeposit());
            continue;
        }

        float jeta = nearest->eta(), jphi = nearest->phi();
        neutralV.push_back(reco::IsoDeposit(jeta, jphi));
        chargedV.push_back(reco::IsoDeposit(jeta, jphi));
        photonV.push_back(reco::IsoDeposit(jeta, jphi));
        reco::IsoDeposit & chargedD = chargedV.back(), & neutralD = neutralV.back(), & photonD = photonV.back(); 

        double dr2max = 9999;
        if (coneMode_ != JetCone) {
            float drjl = reco::deltaR(jeta,jphi,leta,lphi);
            dr2max = (coneMode_ == FixedCone ? coneSize_ : (coneMode_ == MaxCone ? std::max<double>(coneSize_, drjl) : drjl + coneSize_));
            dr2max = dr2max *dr2max;
        }
        size_t nk = (coneMode_ == JetCone ? nearest->numberOfDaughters() : pfH->size());
        for (size_t k = 0; k < nk; ++k) {
            const reco::PFCandidate & pf = (coneMode_ == JetCone ? 
                                             dynamic_cast<const reco::PFCandidate &>(*nearest->daughter(k)) :
                                             pfH->at(k));

            bool charged = pf.trackRef().isNonnull();

            // apply cone and dz cuts when not using PF jet constituents
            if (coneMode_ != JetCone) {
                if (reco::deltaR2<float>(jeta, jphi, pf.eta(), pf.phi()) > dr2max) continue;
                if (charged && std::abs<double>(zLepton - pf.trackRef()->dz(vtxH->at(0).position())) > dzCut_) continue;
            }

            // remove self 
            if (pf.trackRef().isNonnull()    && pf.trackRef()    == candTrackRef)   continue;
            if (pf.gsfTrackRef().isNonnull() && pf.gsfTrackRef() == candGsfTrackRef) continue;

            reco::IsoDeposit &where = (charged ? chargedD : (pf.pdgId() == 22 ? photonD : neutralD));
            where.addDeposit(reco::IsoDeposit::Direction(pf.eta(), pf.phi()), pf.pt());
        }
    }

    neutralF.insert(candH, neutralV.begin(), neutralV.end());
    chargedF.insert(candH, chargedV.begin(), chargedV.end());
    photonF.insert(candH, photonV.begin(), photonV.end());
    neutralF.fill();
    photonF.fill();
    chargedF.fill();
    iEvent.put(neutralM, "neutralHad");
    iEvent.put(chargedM, "charged");
    iEvent.put(photonM,  "photon");
}

PFJetIsoDepMapProd::~PFJetIsoDepMapProd() { }
void PFJetIsoDepMapProd::beginJob() { }
void PFJetIsoDepMapProd::endJob() { }
DEFINE_FWK_MODULE(PFJetIsoDepMapProd);
