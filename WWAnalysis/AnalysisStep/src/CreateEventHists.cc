#include "WWAnalysis/AnalysisStep/interface/CreateEventHists.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "DataFormats/Common/interface/Handle.h"


CreateEventHists::CreateEventHists(const edm::ParameterSet& cfg, TFileDirectory& fs): 
        edm::BasicAnalyzer::BasicAnalyzer(cfg, fs),
        eventFiller_(cfg.getParameter<edm::ParameterSet>("histParams"), fs, cfg.getParameter<std::string>("sampleName")),
        puWeights_(cfg.getParameter<std::vector<double> > ("puWeights")),
        doNMinus1_(cfg.getParameter<bool>("doNMinus1")),
        doByCuts_ (cfg.getParameter<bool>("doByCuts")),
        printEverything_ (cfg.getParameter<bool>("printSummary")),
        myLumiSel_(cfg) {

    // Get a list of all hypotheses
    edm::ParameterSet hypotheses = cfg.getParameter<edm::ParameterSet>("hypotheses");

    // For each PSet in hypotheses, a new hypothesis is made
    hypoNames_ = hypotheses.getParameterNamesForType<edm::ParameterSet>();
    for(size_t hn=0;hn<hypoNames_.size();++hn) {

        // Get the information about this particular hypothesis
        edm::ParameterSet thisPset = hypotheses.getParameter<edm::ParameterSet>(hypoNames_[hn]);

        // Grab the input tag for this hypothesis
        branchTags_.push_back( thisPset.getParameter<edm::InputTag>("src") );

        // Fill the cut labels  and grab the stringselectors for each cut
        std::vector<std::string> cutLabels;
        std::vector<edm::ParameterSet> cuts = thisPset.getParameter<std::vector<edm::ParameterSet> >("cuts");
        for(size_t i=0;i<cuts.size();++i) {
            cutLabels.push_back( cuts[i].getParameter<std::string>("label") );
            stringSelectors_[hypoNames_[hn]].push_back( StringCutObjectSelector<reco::SkimEvent>( cuts[i].getParameter<std::string>("cut") ) );
        }
        eventFiller_.setCutLabels(hypoNames_[hn], cutLabels);
        eventFiller_.setTotalNumberOfCuts(hypoNames_[hn], cutLabels.size());
    }

    if( cfg.existsAs<edm::InputTag>("puLabel") ) {
        puTag_  = cfg.getParameter<edm::InputTag>("puLabel");
    }

    if( cfg.existsAs<edm::InputTag>("ptWeight") ) {
        ptWeightTag_ = cfg.getParameter<edm::InputTag> ("ptWeight");
    }

}


void CreateEventHists::analyze(const edm::EventBase& evt) {

    //setup weight file
    std::vector<float> weights;

    //add PU infomationn
    if( !(puTag_ == edm::InputTag()) ) {
        edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
        evt.getByLabel(puTag_,puInfoH);
        int nPU = 0;
        for(size_t i=0;i<puInfoH->size();++i) {
            if( puInfoH->at(i).getBunchCrossing()==0 ) {
                nPU = puInfoH->at(i).getPU_NumInteractions();
            }
        }
        weights.push_back( puWeights_[std::min(nPU,(int)(puWeights_.size()-1))] );
    }

    //add HqT infomationn
    if( !(ptWeightTag_ == edm::InputTag()) ) {
        edm::Handle<double> ptWeightH;
        evt.getByLabel(ptWeightTag_,ptWeightH);
        weights.push_back((float)*ptWeightH);
    }

    edm::Handle<std::vector<reco::SkimEvent> > skimH;
    if( myLumiSel_(evt) ) {
            for(size_t hn=0;hn<hypoNames_.size();++hn) {

            // get the skim evt from the file
            evt.getByLabel(branchTags_[hn],skimH);

            //loop on no of hypos in evt
            for(std::vector<reco::SkimEvent>::const_iterator mySkimEvent = skimH->begin(); mySkimEvent != skimH->end(); mySkimEvent++){

                size_t instance = mySkimEvent - skimH->begin();


                // loop over cuts for this hypothesis
                for(size_t i=0;i<stringSelectors_[hypoNames_[hn]].size();++i) {
                    if( stringSelectors_[hypoNames_[hn]][i]( (*mySkimEvent) ) ) {
                        eventFiller_(&evt,hypoNames_[hn],instance,i,*mySkimEvent,weights);
                    }
                }

            } //end of loop over SkimEvent

        } // end loop over hypothesis
    } 
//     } else { std::cout << "Skipping lumi : " << evt.luminosityBlock() << std::endl; }

}

void CreateEventHists::endJob() {
    eventFiller_.writeAllYieldHists();
    if(doNMinus1_) eventFiller_.writeAllNMinus1Plots();
    if(doByCuts_) eventFiller_.writeAllByCutPlots();
    if(printEverything_) eventFiller_.printFuckingEverything();
}
