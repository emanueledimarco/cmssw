#include "WWAnalysis/AnalysisStep/interface/PUDumper.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"



PUDumper::PUDumper(const edm::ParameterSet& cfg) {

    fs_ = edm::Service<TFileService>().operator->();
    tree_ = fs_->make<TTree>(cfg.getParameter<std::string>("@module_label").c_str(),cfg.getParameter<std::string>("@module_label").c_str());
    tree_->Branch("nvtx" ,&nvtx_ ,"nvtx/F" );
    tree_->Branch("itpu" ,&itpu_ ,"itpu/F" );
    tree_->Branch("trpu" ,&trpu_ ,"trpu/F" );
    tree_->Branch("ootpu1m" ,&ootpu1m_ ,"ootpu1m/F" );
    tree_->Branch("ootpu2m" ,&ootpu2m_ ,"ootpu2m/F" );
    tree_->Branch("ootpu3m" ,&ootpu3m_ ,"ootpu3m/F" );
    tree_->Branch("ootpu4m" ,&ootpu4m_ ,"ootpu4m/F" );
    tree_->Branch("ootpu1p" ,&ootpu1p_ ,"ootpu1p/F" );
    tree_->Branch("ootpu2p" ,&ootpu2p_ ,"ootpu2p/F" );
    tree_->Branch("ootpu3p" ,&ootpu3p_ ,"ootpu3p/F" );
    tree_->Branch("ootpu4p" ,&ootpu4p_ ,"ootpu4p/F" );
}

void PUDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<std::vector< PileupSummaryInfo > >  puH;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), puH);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float tot = 0.;
    for(PVI = puH->begin(); PVI != puH->end(); ++PVI) {
        int BX = PVI->getBunchCrossing();
        tot += PVI->getPU_NumInteractions();
        if(BX == -4) ootpu4m_ = PVI->getPU_NumInteractions();
        if(BX == -3) ootpu3m_ = PVI->getPU_NumInteractions();
        if(BX == -2) ootpu2m_ = PVI->getPU_NumInteractions();
        if(BX == -1) ootpu1m_ = PVI->getPU_NumInteractions();
        if(BX ==  0) {
         itpu_    = PVI->getPU_NumInteractions();
         trpu_    = PVI->getTrueNumInteractions();
        }
        if(BX ==  1) ootpu1p_ = PVI->getPU_NumInteractions();
        if(BX ==  2) ootpu2p_ = PVI->getPU_NumInteractions();
        if(BX ==  3) ootpu3p_ = PVI->getPU_NumInteractions();
        if(BX ==  4) ootpu4p_ = PVI->getPU_NumInteractions();
    }

    edm::Handle<reco::VertexCollection>  vtxH;
    iEvent.getByLabel(edm::InputTag("goodPrimaryVertices"), vtxH);
    nvtx_ = vtxH->size();

    tree_->Fill();

}

void PUDumper::beginJob() { }
PUDumper::~PUDumper() { } 
void PUDumper::endJob() { }
DEFINE_FWK_MODULE(PUDumper);

