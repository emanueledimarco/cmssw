#include "WWAnalysis/AnalysisStep/interface/SkimTupleProducer.h"


SkimTupleProducer::SkimTupleProducer(const edm::ParameterSet& cfg) : 
        skimTag_(cfg.getUntrackedParameter<edm::InputTag>("src")), 
        plotsInfo_(cfg.getUntrackedParameter<std::vector<edm::ParameterSet> >("variables")),
        weight_(cfg.getUntrackedParameter<double>("weight")) {

    fs_ = edm::Service<TFileService>().operator->();
    tree_ = fs_->make<TTree>(cfg.getParameter<std::string>("@module_label").c_str(),cfg.getParameter<std::string>("@module_label").c_str());
    tree_->Branch("weight",&weight_,"weight/F");

}

void SkimTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

//     std::cout << "Beginning analyze()" << std::endl;
    edm::Handle<std::vector< reco::SkimEvent> > skimH;
    iEvent.getByLabel(skimTag_,skimH);
//     std::cout << "Received SkimEvent" << std::endl;
    

    //Loop Over Each Plot Type
    for(size_t i=0;i<skimH->size();++i) {
        for(size_t j=0;j<plotsInfo_.size();++j) {
            placeHolder_[j] = objFuncs_[j](skimH->at(i));
        }            
//         std::cout << "Looped" << std::endl;
        tree_->Fill();
    }

//     std::cout << "End of analyze()" << std::endl;
}

void SkimTupleProducer::beginJob() { 

    placeHolder_.clear();
    objFuncs_.clear();

//     std::cout << "Before Resize" << std::endl;
    placeHolder_.resize(plotsInfo_.size());
//     std::cout << "After Resize" << std::endl;

    for(size_t j=0;j<plotsInfo_.size();++j) {
        objFuncs_.push_back( StringObjectFunction<reco::SkimEvent>(
                             plotsInfo_[j].getUntrackedParameter<std::string>("quantity")));
        tree_->Branch( plotsInfo_[j].getUntrackedParameter<std::string>("tag").c_str() , 
                       &placeHolder_[j], 
                       (plotsInfo_[j].getUntrackedParameter<std::string>("tag")+"/F").c_str() );
    }

//     std::cout << "Successfully finished beginJob()" << std::endl;

}

SkimTupleProducer::~SkimTupleProducer() { } 
void SkimTupleProducer::endJob() { }
DEFINE_FWK_MODULE(SkimTupleProducer);

