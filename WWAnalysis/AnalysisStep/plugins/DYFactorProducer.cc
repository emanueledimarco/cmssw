#include <memory>
#include <iostream>

#include "WWAnalysis/AnalysisStep/interface/DYFactorProducer.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"


DYFactorProducer::DYFactorProducer(const edm::ParameterSet& iConfig): 
    src_       (iConfig.getParameter<edm::InputTag>("src")),
    weightFile_(edm::FileInPath(iConfig.getParameter<std::string>("weightFile")).fullPath().c_str() ),
    dyWeighter_(new DYWeighter(weightFile_.c_str())) {

    produces<edm::ValueMap<float> > ();
}



void DYFactorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<std::vector<reco::SkimEvent> > skimH;
    iEvent.getByLabel(src_, skimH);

    const reco::GenParticle *gen1;
    const reco::GenParticle *gen2;

    std::vector<float> floatVec(skimH->size(),1);
    for(size_t i=0;i<skimH->size();++i) {
        gen1 = skimH->at(i).genParticle(0);
        gen2 = skimH->at(i).genParticle(1);
        if( gen1 && gen2 ) {
            const math::XYZTLorentzVector &p41 = gen1->p4();
            const math::XYZTLorentzVector &p42 = gen2->p4();
            
            floatVec[i] =  dyWeighter_->GetWeight((p41+p42).mass(),(p41+p42).pt(),(p41+p42).Rapidity());
        }
    }

    std::auto_ptr<edm::ValueMap<float> > floatMap(new edm::ValueMap<float>);
    edm::ValueMap<float>::Filler filler(*floatMap);
    filler.insert(skimH, floatVec.begin(), floatVec.end());
    filler.fill();
    iEvent.put(floatMap);

}

DYFactorProducer::~DYFactorProducer() {}
void DYFactorProducer::beginJob(const edm::EventSetup&) {}
void DYFactorProducer::endJob() { delete dyWeighter_; }
DEFINE_FWK_MODULE(DYFactorProducer);
