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

#include <Math/VectorUtil.h>


class RhoValueMapProducer : public edm::EDProducer {
    public:
        explicit RhoValueMapProducer(const edm::ParameterSet&);
        ~RhoValueMapProducer();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        const edm::InputTag leptonTag_;    
        const edm::InputTag rhoTag_;

};

RhoValueMapProducer::RhoValueMapProducer(const edm::ParameterSet& iConfig):
        leptonTag_(iConfig.getUntrackedParameter<edm::InputTag>("leptonTag")),
        rhoTag_(iConfig.getUntrackedParameter<edm::InputTag>("rhoTag")) {
        
    produces<edm::ValueMap<float> >();
}

void RhoValueMapProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::Candidate> > leptonH;
    iEvent.getByLabel(leptonTag_,  leptonH);

    edm::Handle<double> rhoH;
    iEvent.getByLabel(rhoTag_,rhoH);

    std::auto_ptr<edm::ValueMap<float> > rhoMap(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*rhoMap);
    std::vector<float> myRhos(leptonH->size(),*rhoH);
    filler.insert(leptonH, myRhos.begin(), myRhos.end());
    filler.fill();
    iEvent.put(rhoMap);

}

RhoValueMapProducer::~RhoValueMapProducer() { }
void RhoValueMapProducer::beginJob() { }
void RhoValueMapProducer::endJob() { } 
DEFINE_FWK_MODULE(RhoValueMapProducer);
