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

#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Utilities/interface/Exception.h"


class ChargedMetProducer : public edm::EDProducer {
    public:
        explicit ChargedMetProducer(const edm::ParameterSet&);
        ~ChargedMetProducer();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag collectionTag_;    
        edm::InputTag vertexTag_;
        double dzCut_;
        double minNeutralPt_;
        double maxNeutralEta_;

};

ChargedMetProducer::ChargedMetProducer(const edm::ParameterSet& iConfig):
        collectionTag_(iConfig.getParameter<edm::InputTag>("collectionTag")),
        vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")),
        dzCut_(iConfig.getParameter<double>("dzCut")),
        minNeutralPt_(iConfig.getParameter<double>("minNeutralPt")),
        maxNeutralEta_(iConfig.getParameter<double>("maxNeutralEta")) {
        
    produces<edm::ValueMap<reco::PFMET> >();
}

void ChargedMetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::Vertex> > vtxH;
    iEvent.getByLabel(vertexTag_ ,vtxH);

    if( vtxH->size() == 0 ) {
        throw cms::Exception("CorruptData") << "Weird that the vertex collection is empty, how is this possible?\n";
    }


    edm::Handle<edm::View<reco::Candidate> > leptonH;
    iEvent.getByLabel(collectionTag_,  leptonH);

    std::vector<math::XYZTLorentzVector> p4s(vtxH->size(),math::XYZTLorentzVector(0,0,0,0));
    std::vector<float> sumets(vtxH->size(),0);

    edm::View<reco::Candidate>::const_iterator itLep, endprobes = leptonH->end();
    float tempDz;
    for (itLep = leptonH->begin(); itLep != endprobes; ++itLep) {

        //Neutral case:
        if( itLep->charge() == 0 ) {
            for(size_t j=0;j<vtxH->size();++j) {
                if(  itLep->pt() > minNeutralPt_ && fabs(itLep->eta()) < maxNeutralEta_ ) {
                    p4s[j] += itLep->p4();
                    sumets[j] += itLep->pt();
                }
            }
        //Charged case:
        } else {
            if( dynamic_cast<const reco::RecoCandidate*>(&*itLep) ) {
                tempDz =  dynamic_cast<const reco::RecoCandidate*>(&*itLep)->track().isNonnull() ?
                          dynamic_cast<const reco::RecoCandidate*>(&*itLep)->track()->vz() :
                          dynamic_cast<const reco::RecoCandidate*>(&*itLep)->gsfTrack()->vz() ;
            } else if( dynamic_cast<const reco::PFCandidate*>(&*itLep) ) {
                tempDz =  dynamic_cast<const reco::PFCandidate*>(&*itLep)->vz();
            } else throw cms::Exception("InvalidData") << "Data pointed to by collectionTag must inherit from reco::RecoCandidate or reco::PFCandidate\n";

            float closeDz = fabs(tempDz - vtxH->at(0).z());
            size_t iClose = 0;
            for(size_t j=1;j<vtxH->size();++j) {
                if( fabs(tempDz - vtxH->at(j).z()) < closeDz ) {
                    closeDz = fabs(tempDz - vtxH->at(j).z());
                    iClose = j;
                }
            }
            if( fabs(tempDz - vtxH->at(iClose).z()) < dzCut_ ) {
                p4s[iClose] += itLep->p4();
                sumets[iClose] += itLep->pt();
            }
        }

    }

    std::vector<reco::PFMET> mets;
    CommonMETData output;
    for(size_t i=0;i<p4s.size();++i) {
        output.mex = -p4s[i].px();
        output.mey = -p4s[i].py();
        output.mez = -p4s[i].pz();
        output.met = p4s[i].pt();
        output.sumet = sumets[i];
        output.phi = atan2(-p4s[i].py(),-p4s[i].px());

        PFSpecificAlgo pf;
        mets.push_back(pf.addInfo(leptonH,output));
    }


    //Get original vertex collection
    edm::Handle<edm::View<reco::Vertex> > originalVtxH;
    iEvent.get( vtxH->refAt(0).id() , originalVtxH );

    // set up the ValueMap for adding to the event
    std::auto_ptr<edm::ValueMap<reco::PFMET> > metMap(new edm::ValueMap<reco::PFMET>);

    // blank set of pfmets (one for each original vertex
    std::vector<reco::PFMET> emptymets(originalVtxH->size(),reco::PFMET());

    //fill the value map
    edm::ValueMap<reco::PFMET>::Filler filler(*metMap);
    filler.insert(originalVtxH, emptymets.begin(), emptymets.end());
    filler.fill();

    // add the ones I want to the ValueMap
    for(size_t i=0;i<vtxH->size();++i) {
        for(size_t j=0;j<originalVtxH->size();++j) {
            if( vtxH->refAt(i) == originalVtxH->refAt(j) ) {
                (*metMap)[originalVtxH->refAt(j)] = mets[i];
            }
        }
    }

    // and finally put it in the event
    iEvent.put(metMap);

}

ChargedMetProducer::~ChargedMetProducer() { }
void ChargedMetProducer::beginJob() { }
void ChargedMetProducer::endJob() { } 
DEFINE_FWK_MODULE(ChargedMetProducer);
