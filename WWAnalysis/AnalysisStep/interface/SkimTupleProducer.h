#ifndef MYPRODUCERS_SKIMTUPLEPRODUCER_SKIMTUPLEPRODUCER_H
#define MYPRODUCERS_SKIMTUPLEPRODUCER_SKIMTUPLEPRODUCER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TTree.h>
#include <vector>
#include <string>


class SkimTupleProducer : public edm::EDAnalyzer {
    public:
        explicit SkimTupleProducer(const edm::ParameterSet&);
        ~SkimTupleProducer();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        edm::InputTag skimTag_;
        std::vector<edm::ParameterSet> plotsInfo_;
        float weight_;
        std::vector< StringObjectFunction< reco::SkimEvent > > objFuncs_;
        std::vector<float> placeHolder_;
        TTree *tree_;
        TFileService *fs_;

};

#endif
