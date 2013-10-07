#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"


class GioTriggerFilter : public edm::EDFilter {
    public:
        explicit GioTriggerFilter(const edm::ParameterSet&);
        ~GioTriggerFilter();

    private:
        virtual void beginJob() ;
        virtual bool filter(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag triggerTag_;
        TriggerBitChecker checker_;

};

GioTriggerFilter::GioTriggerFilter(const edm::ParameterSet& iConfig) :
    triggerTag_(iConfig.getParameter<edm::InputTag>("triggerTag")),
    checker_(iConfig.getParameter<std::vector<std::string> > ("pathsToKeep"))
{ } 


bool GioTriggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerTag_,triggerResults);

    return checker_.check(iEvent,*triggerResults);

}

GioTriggerFilter::~GioTriggerFilter() { }
void GioTriggerFilter::beginJob() { }
void GioTriggerFilter::endJob() { }
DEFINE_FWK_MODULE(GioTriggerFilter);


