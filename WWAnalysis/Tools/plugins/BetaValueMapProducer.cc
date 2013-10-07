// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"



#include <Math/VectorUtil.h>


class BetaValueMapProducer : public edm::EDProducer {
    public:
        explicit BetaValueMapProducer(const edm::ParameterSet&);
        ~BetaValueMapProducer();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        const edm::InputTag leptonTag_;    
        const edm::InputTag trackTag_;
        double etaValue_;
        double dRValue_;
        double outerConeValue_;
        double ptCut_;

};

BetaValueMapProducer::BetaValueMapProducer(const edm::ParameterSet& iConfig):
        leptonTag_(iConfig.getUntrackedParameter<edm::InputTag>("leptonTag")),
        trackTag_(iConfig.getUntrackedParameter<edm::InputTag>("trackTag")),
        etaValue_(iConfig.getUntrackedParameter<double>("dEtaVeto")),
        dRValue_(iConfig.getUntrackedParameter<double>("dRVeto")),
        outerConeValue_(iConfig.getUntrackedParameter<double>("outerCone")),
        ptCut_(iConfig.getUntrackedParameter<double>("ptCut")) {
        
    produces<edm::ValueMap<float> >();
}

void BetaValueMapProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::RecoCandidate> > leptonH;
    iEvent.getByLabel(leptonTag_,  leptonH);

    edm::Handle<std::vector<reco::Track> > tkH;
    iEvent.getByLabel(trackTag_ ,tkH);

    std::vector<float> betas;

    // fill
    edm::View<reco::RecoCandidate>::const_iterator itLep, endprobes = leptonH->end();
    for (itLep = leptonH->begin(); itLep != endprobes; ++itLep) {

        if( !itLep->track().isNonnull() && !itLep->gsfTrack().isNonnull() ) {
            betas.push_back(1);
            continue;
        }

        math::XYZVector mom =  itLep->track().isNonnull() ? itLep->track()->momentum() : itLep->gsfTrack()->momentum() ;
        float num = 0;
        float den = 0;
        for(size_t k=0;k<tkH->size();++k) {
            reco::TrackRef tk(tkH,k);
            

            if( tk->pt() < ptCut_) continue;
            if( fabs(tk->eta()-mom.eta()) < etaValue_) continue;
            if( ROOT::Math::VectorUtil::DeltaR(tk->momentum(),mom) < dRValue_ ) continue;

            if( ROOT::Math::VectorUtil::DeltaR(tk->momentum(),mom) < outerConeValue_ ) {
                den += tk->pt();
                float dZ =  itLep->track().isNonnull() ? itLep->track()->dz(tk->vertex()) : itLep->gsfTrack()->dz(tk->vertex()) ;
                if( fabs( dZ ) < 0.2 )
                    num += tk->pt();
            }
        }
        betas.push_back(((den==0)?1:num/den));

    }

    // convert into ValueMap and store -- betaMap
    std::auto_ptr<edm::ValueMap<float> > betaMap(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*betaMap);
    filler.insert(leptonH, betas.begin(), betas.end());
    filler.fill();
    iEvent.put(betaMap);


}

BetaValueMapProducer::~BetaValueMapProducer() { }
void BetaValueMapProducer::beginJob() { }
void BetaValueMapProducer::endJob() { } 
DEFINE_FWK_MODULE(BetaValueMapProducer);
