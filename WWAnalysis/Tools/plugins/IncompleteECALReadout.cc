#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


class IncompleteECALReadout : public edm::EDFilter {

    public:
        explicit IncompleteECALReadout (const edm::ParameterSet&);
        ~IncompleteECALReadout();

    private:
        void beginJob();
        void endJob();
        bool filter(edm::Event&, const edm::EventSetup&);

};

IncompleteECALReadout::IncompleteECALReadout(const edm::ParameterSet& iConfig) {}
IncompleteECALReadout::~IncompleteECALReadout() {} 
void IncompleteECALReadout::beginJob() {}
void IncompleteECALReadout::endJob() {}

bool IncompleteECALReadout::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
    edm::Handle<EcalRecHitCollection> pBarrelEcalRecHits ;
    iEvent.getByLabel ("ecalRecHit","EcalRecHitsEB", pBarrelEcalRecHits) ;
    const EcalRecHitCollection* theBarrelEcalRecHits = pBarrelEcalRecHits.product () ;

    edm::Handle<EcalRecHitCollection> pEndcapEcalRecHits ;
    iEvent.getByLabel ("ecalRecHit","EcalRecHitsEE", pEndcapEcalRecHits) ;
    const EcalRecHitCollection* theEndcapEcalRecHits = pEndcapEcalRecHits.product () ;


    //if((theBarrelEcalRecHits->size() == 0 && theEndcapEcalRecHits->size() != 0) ||
    //   (theBarrelEcalRecHits->size() != 0 && theEndcapEcalRecHits->size() == 0))
    if((theBarrelEcalRecHits->size() == 0 || theEndcapEcalRecHits->size() == 0))
    {
        return false;
    }

    else
    {
        return true;  
    }
}


DEFINE_FWK_MODULE(IncompleteECALReadout);

