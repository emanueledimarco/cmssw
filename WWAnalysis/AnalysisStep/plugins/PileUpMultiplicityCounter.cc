#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class PileUpMultiplicityCounter : public edm::EDProducer {
    public:
        explicit PileUpMultiplicityCounter(const edm::ParameterSet&);
        ~PileUpMultiplicityCounter();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
    
        edm::InputTag puTag_;
        bool hasSrc_;
        edm::InputTag src_;
};

PileUpMultiplicityCounter::PileUpMultiplicityCounter(const edm::ParameterSet& iConfig) :
    puTag_(iConfig.getParameter<edm::InputTag>("puLabel")),
    hasSrc_(iConfig.existsAs<edm::InputTag>("src")),
    src_(hasSrc_ ? iConfig.getParameter<edm::InputTag>("src") : edm::InputTag())
{ 
    if (hasSrc_) {
        produces<edm::ValueMap<float> > ("it");
        produces<edm::ValueMap<float> > ("m1");
        produces<edm::ValueMap<float> > ("p1");
    } else {
        produces<float> ("it");
        produces<float> ("m1");
        produces<float> ("p1");
    }
} 


void PileUpMultiplicityCounter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    std::auto_ptr<float> it(new float(0));
    std::auto_ptr<float> m1(new float(0));
    std::auto_ptr<float> p1(new float(0));

    edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
    iEvent.getByLabel(puTag_,puInfoH);
    for(size_t i=0;i<puInfoH->size();++i) {
        if( puInfoH->at(i).getBunchCrossing()==0 ) *it = puInfoH->at(i).getPU_NumInteractions();
        if( puInfoH->at(i).getBunchCrossing()==-1) *m1 = puInfoH->at(i).getPU_NumInteractions();
        if( puInfoH->at(i).getBunchCrossing()==1 ) *p1 = puInfoH->at(i).getPU_NumInteractions();
    }

    if(hasSrc_) {
        edm::Handle<edm::View<reco::Candidate> > srcH;
        iEvent.getByLabel(src_ ,srcH);

        std::auto_ptr<edm::ValueMap<float> > itMap(new edm::ValueMap<float>);
        std::auto_ptr<edm::ValueMap<float> > m1Map(new edm::ValueMap<float>);
        std::auto_ptr<edm::ValueMap<float> > p1Map(new edm::ValueMap<float>);

        std::vector<float> itVec(srcH->size(),*it);
        std::vector<float> m1Vec(srcH->size(),*m1);
        std::vector<float> p1Vec(srcH->size(),*p1);

        edm::ValueMap<float>::Filler itFill(*itMap);
        edm::ValueMap<float>::Filler m1Fill(*m1Map);
        edm::ValueMap<float>::Filler p1Fill(*p1Map);

        itFill.insert(srcH, itVec.begin(), itVec.end());
        m1Fill.insert(srcH, m1Vec.begin(), m1Vec.end());
        p1Fill.insert(srcH, p1Vec.begin(), p1Vec.end());

        itFill.fill();
        m1Fill.fill();
        p1Fill.fill();

        iEvent.put(itMap,"it");
        iEvent.put(p1Map,"p1");
        iEvent.put(m1Map,"m1");
    } else {
        iEvent.put(it,"it");
        iEvent.put(p1,"p1");
        iEvent.put(m1,"m1");
    }
}

PileUpMultiplicityCounter::~PileUpMultiplicityCounter() { }
DEFINE_FWK_MODULE(PileUpMultiplicityCounter);


