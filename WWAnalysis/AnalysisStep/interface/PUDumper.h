#ifndef MYPRODUCERS_PUDumper_PUDumper_H
#define MYPRODUCERS_PUDumper_PUDumper_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "WWAnalysis/AnalysisStep/interface/LatinoReWeighting.h"

#include <TTree.h>
#include <vector>
#include <string>


class PUDumper : public edm::EDAnalyzer {
    public:
        explicit PUDumper(const edm::ParameterSet&);
        ~PUDumper();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        TFileService *fs_;
        TTree *tree_;
        float nvtx_;
        float itpu_;
        float trpu_ ;
        float ootpu1m_;
        float ootpu2m_;
        float ootpu3m_;
        float ootpu4m_;
        float ootpu1p_;
        float ootpu2p_;
        float ootpu3p_;
        float ootpu4p_;
};

#endif
