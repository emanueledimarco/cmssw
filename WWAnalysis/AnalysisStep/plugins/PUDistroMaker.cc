#include "WWAnalysis/AnalysisStep/interface/PUDistroMaker.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"



PUDistroMaker::PUDistroMaker(const edm::ParameterSet& cfg) :
        s3_(cfg.getParameter<std::vector<double> >("s3")) ,
        s4_(cfg.getParameter<std::vector<double> >("s4")) ,
        lumi_(cfg.getParameter<std::vector<double> >("lumi")) ,
        lrw_(new edm::LatinoReWeighting(s3_,s4_,lumi_)) {

    fs_ = edm::Service<TFileService>().operator->();
    tree_ = fs_->make<TTree>(cfg.getParameter<std::string>("@module_label").c_str(),cfg.getParameter<std::string>("@module_label").c_str());
    tree_->Branch("nvtx_" ,&nvtx_ ,"nvtx_/F" );
    tree_->Branch("pun1_" ,&pun1_ ,"pun1_/F" );
    tree_->Branch("pu0_"  ,&pu0_  ,"pu0_/F"  );
    tree_->Branch("pup1_" ,&pup1_ ,"pup1_/F" );
    tree_->Branch("puavg_",&puavg_,"puavg_/F");
    tree_->Branch("weight_",&weight_,"weight_/F");
    tree_->Branch("weight3BX_",&weight3BX_,"weight3BX_/F");
    tree_->Branch("weightOOT_",&weightOOT_,"weightOOT_/F");
    tree_->Branch("weightLatino_",&weightLatino_,"weightLatino_/F");
}

void PUDistroMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<std::vector< PileupSummaryInfo > >  puH;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), puH);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float tot = 0.;
    for(PVI = puH->begin(); PVI != puH->end(); ++PVI) {
        int BX = PVI->getBunchCrossing();
        tot += PVI->getPU_NumInteractions();
        if(BX == -1) pun1_ = PVI->getPU_NumInteractions();
        if(BX ==  0) pu0_  = PVI->getPU_NumInteractions();
        if(BX ==  1) pup1_ = PVI->getPU_NumInteractions();
    }
    puavg_ = tot/3.;

    edm::Handle<reco::VertexCollection>  vtxH;
    iEvent.getByLabel(edm::InputTag("goodPrimaryVertices"), vtxH);
    nvtx_ = vtxH->size();

//     weight_ = lrw_->weight(iEvent);
//     weight3BX_ = lrw_->weight3BX(iEvent);
//     weightOOT_ = lrw_->weightOOT(iEvent);
    weightLatino_ = lrw_->fullWeight(iEvent);

    tree_->Fill();

}

void PUDistroMaker::beginJob() { }
PUDistroMaker::~PUDistroMaker() { } 
void PUDistroMaker::endJob() { }
DEFINE_FWK_MODULE(PUDistroMaker);

