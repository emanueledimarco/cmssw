#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

class WgammaFilter : public edm::EDFilter {
    public:
        explicit WgammaFilter(const edm::ParameterSet&);
        ~WgammaFilter();

    private:
        virtual void beginJob() ;
        virtual bool filter(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

};

bool WgammaFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace std;
    using namespace edm;
    using namespace reco;

    Handle<GenParticleCollection> particles;
    iEvent.getByLabel("prunedGen",particles);

    bool low = false;
    bool high = false;
    for(GenParticleCollection::const_iterator it=particles->begin(); it!=particles->end(); ++it){
        if(it->status()!=3) continue;
        if(fabs(it->pdgId())!=11 && fabs(it->pdgId())!=13) continue;
        for(GenParticleCollection::const_iterator it2=it+1; it2!=particles->end(); ++it2){
            if(it2->status()!=3) continue;
            if(fabs(it2->pdgId())!=11 && fabs(it2->pdgId())!=13) continue;

            if( it->pdgId() * it2->pdgId() > 0 ) continue;

            if( (it->p4() + it2->p4()).mass() > 12.) high = true;
            else                                     low  = true;
        }
    }
    bool skip = high && !low ;


    return (!skip);
}


WgammaFilter::WgammaFilter(const edm::ParameterSet& iConfig) { } 
WgammaFilter::~WgammaFilter() { }
void WgammaFilter::beginJob() { }
void WgammaFilter::endJob() { }

DEFINE_FWK_MODULE(WgammaFilter);
