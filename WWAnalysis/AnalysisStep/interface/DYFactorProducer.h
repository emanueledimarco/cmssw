#ifndef HWWKFACTORPRODUCER
#define HWWKFACTORPRODUCER

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "WWAnalysis/AnalysisStep/interface/DYWeighter.h"

class DYFactorProducer : public edm::EDProducer {
    public:
        explicit DYFactorProducer(const edm::ParameterSet&);
        ~DYFactorProducer();


    private:
        virtual void beginJob(const edm::EventSetup&) ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------

        edm::InputTag    src_;
        std::string weightFile_;
        DYWeighter *dyWeighter_;
};

#endif
