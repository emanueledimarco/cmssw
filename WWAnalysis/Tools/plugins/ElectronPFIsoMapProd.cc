#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

class ElectronPFIsoMapProd : public edm::EDProducer {
    public:
        explicit ElectronPFIsoMapProd(const edm::ParameterSet&);
        ~ElectronPFIsoMapProd();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag vtxLabel_;
        edm::InputTag eleLabel_;
        edm::InputTag pfLabel_;
        double deltaR_;
};



ElectronPFIsoMapProd::ElectronPFIsoMapProd(const edm::ParameterSet& iConfig) :
    vtxLabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxLabel")),
    eleLabel_(iConfig.getUntrackedParameter<edm::InputTag>("eleLabel")),
    pfLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfLabel")),
    deltaR_(iConfig.getUntrackedParameter<double>("deltaR")) {

    produces<edm::ValueMap<float> >().setBranchAlias("pfElIso");
}

void ElectronPFIsoMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::GsfElectronCollection> eleH;
    iEvent.getByLabel(eleLabel_,eleH);

    edm::Handle<reco::PFCandidateCollection> pfH;
    iEvent.getByLabel(pfLabel_,pfH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxLabel_,vtxH);

    std::vector<float> isoV;
    std::auto_ptr<edm::ValueMap<float> > isoM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler isoF(*isoM);

    for(size_t i=0; i<eleH->size();++i) {
        const reco::GsfElectron &ele = eleH->at(i);

        Double_t zLepton = 0.0;
        if(ele.track().isNonnull()) zLepton = ele.track()->dz(vtxH->at(0).position());

        Double_t ptSum =0.;  
        for (size_t j=0; j<pfH->size();j++) {   
            const reco::PFCandidate &pf = pfH->at(j);

            // supercluster veto on electrons failing missing hits cut
            if (ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits()>0 && pf.mva_nothing_gamma() > 0.99 && 
                ele.superCluster().isNonnull() && pf.superClusterRef().isNonnull() && 
                ele.superCluster() == pf.superClusterRef()) continue;

            // exclude electron
            if(pf.trackRef().isNonnull() && ele.closestCtfTrackRef().isNonnull() &&
                    pf.trackRef() == ele.closestCtfTrackRef()) continue;

            // exclude electron
            if(pf.gsfTrackRef().isNonnull() && ele.gsfTrack().isNonnull() &&
                    pf.gsfTrackRef() == ele.gsfTrack()) continue;


            // pt cut applied to neutrals
            if(!pf.trackRef().isNonnull() && pf.pt() <= 1.0) continue;

            // ignore the pf candidate if it is too far away in Z
            Double_t deltaZ = 0.0;
            if(pf.trackRef().isNonnull()) {
                deltaZ = fabs(pf.trackRef()->dz(vtxH->at(0).position()) - zLepton);
            }
            if (deltaZ >= 0.1) continue;

//             Double_t dr = ROOT::Math::VectorUtil::DeltaR(ele.momentum(), pf.momentum());
//             // add the pf pt if it is inside the extRadius 
//             if ( dr < 0.4 && dr > 0.0) ptSum += pf.pt();

//             // inner cone veto for gammas
//             if (pf.particleId() == reco::PFCandidate::gamma && dr < 0.0) continue;

            Double_t dr = ROOT::Math::VectorUtil::DeltaR(ele.momentum(), pf.momentum());
            // add the pf pt if it is inside the extRadius and outside the intRadius
            if ( dr < deltaR_ && dr >= 0.0 ) {

                //EtaStrip Veto for Gamma 
                if (pf.particleId() == reco::PFCandidate::gamma && fabs(ele.eta() - pf.eta()) < 0.025) continue;

                //InnerCone (One Tower = dR < 0.07) Veto for non-gamma neutrals
                if (!pf.trackRef().isNonnull() && pf.particleId() == reco::PFCandidate::h0
                        && ROOT::Math::VectorUtil::DeltaR(ele.momentum(), pf.momentum()) < 0.07 ) continue; 

                ptSum += pf.pt();            

            }



        }
        isoV.push_back(ptSum);
    }

    isoF.insert(eleH,isoV.begin(),isoV.end());

    isoF.fill();
    iEvent.put(isoM);

}

ElectronPFIsoMapProd::~ElectronPFIsoMapProd() { }
void ElectronPFIsoMapProd::beginJob() { }
void ElectronPFIsoMapProd::endJob() { }
DEFINE_FWK_MODULE(ElectronPFIsoMapProd);
