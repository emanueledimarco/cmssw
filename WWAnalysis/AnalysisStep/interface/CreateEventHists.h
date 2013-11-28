#ifndef WWANALYSIS_ANALYSISSTEP_CREATEEVENTHISTS_H
#define WWANALYSIS_ANALYSISSTEP_CREATEEVENTHISTS_H

#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "PhysicsTools/SelectorUtils/interface/RunLumiSelector.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "WWAnalysis/AnalysisStep/interface/EventBitsFiller.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

#include<vector>
#include<map>

class CreateEventHists : public edm::BasicAnalyzer {

    public:
        CreateEventHists(const edm::ParameterSet& cfg, TFileDirectory& fs);
        virtual ~CreateEventHists(){};

        void beginJob(){};
        void endJob();
        void analyze(const edm::EventBase& event);

    private:
        EventBitsFiller eventFiller_;
        std::vector<double> puWeights_;
        edm::InputTag puTag_;
        edm::InputTag ptWeightTag_;
        bool doNMinus1_;
        bool doByCuts_;
        bool printEverything_;
        RunLumiSelector myLumiSel_;
        std::vector<std::string> hypoNames_;
        std::vector<edm::InputTag> branchTags_;
        std::map<std::string,std::vector<StringCutObjectSelector<reco::SkimEvent> > > stringSelectors_;

};



#endif


