#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

class MuonPFIsoMapProd : public edm::EDProducer {
    public:
        explicit MuonPFIsoMapProd(const edm::ParameterSet&);
        ~MuonPFIsoMapProd();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag vtxLabel_;
        edm::InputTag muonLabel_;
        edm::InputTag pfLabel_;
        double deltaR_;
};



MuonPFIsoMapProd::MuonPFIsoMapProd(const edm::ParameterSet& iConfig) :
    vtxLabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxLabel")),
    muonLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonLabel")),
    pfLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfLabel")),
    deltaR_(iConfig.getUntrackedParameter<double>("deltaR")) {

    produces<edm::ValueMap<float> >().setBranchAlias("pfMuIso");
}

void MuonPFIsoMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::MuonCollection> muH;
    iEvent.getByLabel(muonLabel_,muH);

    edm::Handle<reco::PFCandidateCollection> pfH;
    iEvent.getByLabel(pfLabel_,pfH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxLabel_,vtxH);

    std::vector<float> isoV;
    std::auto_ptr<edm::ValueMap<float> > isoM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler isoF(*isoM);

    for(size_t i=0; i<muH->size();++i) {
        const reco::Muon &mu = muH->at(i);

        Double_t zLepton = 0.0;
        if(mu.track().isNonnull()) zLepton = mu.track()->dz(vtxH->at(0).position());

        Double_t ptSum =0.;  
        for (size_t j=0; j<pfH->size();j++) {   
            const reco::PFCandidate &pf = pfH->at(j);

            // exclude muon
            if(pf.trackRef().isNonnull() && mu.track().isNonnull() &&
                    pf.trackRef() == mu.track()) continue;


            // pt cut applied to neutrals
            if(!pf.trackRef().isNonnull() && pf.pt() <= 1.0) continue;

            // ignore the pf candidate if it is too far away in Z
            Double_t deltaZ = 0.0;
            if(pf.trackRef().isNonnull()) {
                deltaZ = fabs(pf.trackRef()->dz(vtxH->at(0).position()) - zLepton);
            }
            if (deltaZ >= 0.1) continue;

            // inner cone veto for gammas
            Double_t dr = ROOT::Math::VectorUtil::DeltaR(mu.momentum(), pf.momentum());
            if (pf.particleId() == reco::PFCandidate::gamma && dr < 0.0) continue;

            // add the pf pt if it is inside the extRadius 
            if ( dr < deltaR_ ) ptSum += pf.pt();

        }
        isoV.push_back(ptSum);
    }

    isoF.insert(muH,isoV.begin(),isoV.end());

    isoF.fill();
    iEvent.put(isoM);

}

MuonPFIsoMapProd::~MuonPFIsoMapProd() { }
void MuonPFIsoMapProd::beginJob() { }
void MuonPFIsoMapProd::endJob() { }
DEFINE_FWK_MODULE(MuonPFIsoMapProd);
