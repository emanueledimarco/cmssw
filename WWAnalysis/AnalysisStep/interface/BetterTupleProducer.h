#ifndef MYPRODUCERS_SKIMTUPLEPRODUCER_SKIMTUPLEPRODUCER_H
#define MYPRODUCERS_SKIMTUPLEPRODUCER_SKIMTUPLEPRODUCER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TTree.h>
#include <vector>
#include <string>


template <class T> class BetterTupleProducer : public edm::EDAnalyzer {
    public:
        explicit BetterTupleProducer(const edm::ParameterSet&);
        ~BetterTupleProducer();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        edm::InputTag skimTag_;
        std::vector<edm::ParameterSet> plotsInfo_;

        TFileService *fs_;
        TTree *tree_;

        float weight_;

        edm::InputTag vtxTag_;
        int nVtx_;

        std::vector<edm::InputTag> countTags_;
        std::vector<int> nCounts_;

        std::vector<edm::InputTag> rhoTags_;
        std::vector<float> nRhos_;

        edm::InputTag puTag_;
        int nPU_;

        edm::ParameterSet vtxWeightInfo_;
        std::vector<std::vector<double> > vtxWeights_;
        std::vector<float> specificVtxWeight_;

        std::vector< StringObjectFunction< T > > objFuncs_;
        std::vector<float> placeHolder_;

};

#endif

