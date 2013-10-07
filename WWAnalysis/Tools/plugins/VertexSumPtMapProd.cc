#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/Common/interface/ValueMap.h>

#include <vector>

class VertexSumPtMapProd : public edm::EDProducer {
    public:
        explicit VertexSumPtMapProd(const edm::ParameterSet&);
        ~VertexSumPtMapProd();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag vertexTag_;
};

VertexSumPtMapProd::VertexSumPtMapProd(const edm::ParameterSet& iConfig) :
    vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")) {
    produces<edm::ValueMap<float> >("sumPt").setBranchAlias("pvSumPt");
    produces<edm::ValueMap<float> >("sumPt2").setBranchAlias("pvSumPt");
}




void VertexSumPtMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vertexTag_,vtxH);

    std::vector<float> sumPtsV;
    std::auto_ptr<edm::ValueMap<float> > sumPtsM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler sumPtsF(*sumPtsM);

    std::vector<float> sumPts2V;
    std::auto_ptr<edm::ValueMap<float> > sumPts2M(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler sumPts2F(*sumPts2M);

    float tempSum = 0., tempSum2 = 0.;
    reco::VertexCollection::const_iterator itVtx;
    reco::Vertex::trackRef_iterator itTk;
    for(itVtx=vtxH->begin();itVtx!=vtxH->end();++itVtx) {

        tempSum = 0.;
        for ( itTk = itVtx->tracks_begin(); itTk != itVtx->tracks_end(); itTk++) {
            tempSum += (**itTk).pt();
            tempSum2 += (**itTk).pt()*(**itTk).pt();
        }
        sumPtsV.push_back(tempSum);
        sumPts2V.push_back(tempSum2);


    }

    sumPtsF.insert(vtxH,sumPtsV.begin(),sumPtsV.end());
    sumPtsF.fill();
    iEvent.put(sumPtsM,"sumPt");

    sumPts2F.insert(vtxH,sumPts2V.begin(),sumPts2V.end());
    sumPts2F.fill();
    iEvent.put(sumPts2M,"sumPt2");

}

VertexSumPtMapProd::~VertexSumPtMapProd() { }
void VertexSumPtMapProd::beginJob() { }
void VertexSumPtMapProd::endJob() { }
DEFINE_FWK_MODULE(VertexSumPtMapProd);
