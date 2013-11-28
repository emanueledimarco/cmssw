// system include files
#include <memory>
#include <cmath>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"



class InterestingVertexRefProducer : public edm::EDProducer {
    public:
        explicit InterestingVertexRefProducer(const edm::ParameterSet&);
        ~InterestingVertexRefProducer();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        std::vector<edm::InputTag> leptonTags_;    
        edm::InputTag vertexTag_;

};

InterestingVertexRefProducer::InterestingVertexRefProducer(const edm::ParameterSet& iConfig):
        leptonTags_(iConfig.getParameter<std::vector<edm::InputTag> >("leptonTags")),
        vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")) {
        
    produces<edm::RefVector<reco::VertexCollection> >();
}

void InterestingVertexRefProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::VertexCollection > vtxH;
    iEvent.getByLabel(vertexTag_ ,vtxH);

    std::vector<size_t> vtxIds(1,0);

    edm::Handle<edm::View<reco::Candidate> > leptonH;
    for(size_t i=0;i<leptonTags_.size();++i) {
        iEvent.getByLabel(leptonTags_[i],  leptonH);

        edm::View<reco::Candidate>::const_iterator itLep, endprobes = leptonH->end();
        for (itLep = leptonH->begin(); itLep != endprobes; ++itLep) {

            const reco::RecoCandidate *temp = dynamic_cast<const reco::RecoCandidate*>(&*itLep);

            float minDZ = temp->track().isNonnull() ? temp->track()->dz(vtxH->at(0).position()) : temp->gsfTrack()->dz(vtxH->at(0).position()) ;
            size_t minIndex = 0;
            for(size_t j=1;j<vtxH->size();++j) {
                float thisDZ = temp->track().isNonnull() ? temp->track()->dz(vtxH->at(j).position()) : temp->gsfTrack()->dz(vtxH->at(j).position()) ;
                if(  fabs(thisDZ) < fabs(minDZ) ) {
                    minIndex = j;
                    minDZ = thisDZ;
                }
            }
            vtxIds.push_back(minIndex);

        }
    }

    std::sort(vtxIds.begin(), vtxIds.end());
    std::vector<size_t>::iterator end = std::unique(vtxIds.begin(), vtxIds.end());

    std::auto_ptr<edm::RefVector<reco::VertexCollection> > myRefVector(new edm::RefVector<reco::VertexCollection>);
    for(std::vector<size_t>::const_iterator it=vtxIds.begin();it!=end;++it) {
        myRefVector->push_back(reco::VertexRef(vtxH,it-vtxIds.begin()));
    }

    iEvent.put(myRefVector);


}

InterestingVertexRefProducer::~InterestingVertexRefProducer() { }
void InterestingVertexRefProducer::beginJob() { }
void InterestingVertexRefProducer::endJob() { } 
DEFINE_FWK_MODULE(InterestingVertexRefProducer);
