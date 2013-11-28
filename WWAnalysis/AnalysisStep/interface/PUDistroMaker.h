#ifndef MYPRODUCERS_PUDistroMaker_PUDistroMaker_H
#define MYPRODUCERS_PUDistroMaker_PUDistroMaker_H

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


class PUDistroMaker : public edm::EDAnalyzer {
    public:
        explicit PUDistroMaker(const edm::ParameterSet&);
        ~PUDistroMaker();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        TFileService *fs_;
        TTree *tree_;
        float nvtx_;
        float pun1_;
        float pu0_ ;
        float pup1_;
        float puavg_;
        float weight_;
        float weight3BX_;
        float weightOOT_;
        float weightLatino_;
        std::vector<double> s3_;
        std::vector<double> s4_;
        std::vector<double> lumi_;
        edm::LatinoReWeighting *lrw_;

};

#endif
