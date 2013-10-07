#include "WWAnalysis/AnalysisStep/interface/SkimEventProducer2L2N.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include<vector> 
#include "Math/VectorUtil.h"


SkimEventProducer2L2N::SkimEventProducer2L2N(const edm::ParameterSet& cfg) :
    //   branchAlias_(cfg.getParameter<std::string>("branchAlias")),
    hypoType_(reco::SkimEvent::hypoTypeByName(cfg.getParameter<std::string>("hypoType"))),
    triggerTag_(cfg.getParameter<edm::InputTag>("triggerTag")),
    singleMuData_  ( cfg.getParameter<std::vector<std::string> >("singleMuDataPaths") ),
    singleElData_  ( cfg.getParameter<std::vector<std::string> >("singleElDataPaths") ),
    doubleMuData_  ( cfg.getParameter<std::vector<std::string> >("doubleMuDataPaths") ),
    doubleElData_  ( cfg.getParameter<std::vector<std::string> >("doubleElDataPaths") ),
    muEGData_      ( cfg.getParameter<std::vector<std::string> >("muEGDataPaths") ),
    singleMuMC_    ( cfg.getParameter<std::vector<std::string> >("singleMuMCPaths") ),
    singleElMC_    ( cfg.getParameter<std::vector<std::string> >("singleElMCPaths") ),
    doubleMuMC_    ( cfg.getParameter<std::vector<std::string> >("doubleMuMCPaths") ),
    doubleElMC_    ( cfg.getParameter<std::vector<std::string> >("doubleElMCPaths") ),
    muEGMC_        ( cfg.getParameter<std::vector<std::string> >("muEGMCPaths") ), 
    AllEmbed_      ( cfg.getParameter<std::vector<std::string> >("AllEmbedPaths") )
{
    mcGenEventInfoTag_ = cfg.getParameter<edm::InputTag>("mcGenEventInfoTag"); 
    mcGenWeightTag_    = cfg.getParameter<edm::InputTag>("mcGenWeightTag"); 
    genParticlesTag_   = cfg.getParameter<edm::InputTag>("genParticlesTag"); 
    muTag_             = cfg.getParameter<edm::InputTag>("muTag"     ); 
    elTag_             = cfg.getParameter<edm::InputTag>("elTag"     ); 
    softMuTag_         = cfg.getParameter<edm::InputTag>("softMuTag" ); 
    jetTag_            = cfg.getParameter<edm::InputTag>("jetTag"    ); 
    tagJetTag_         = cfg.getParameter<edm::InputTag>("tagJetTag" ); 
    pfMetTag_          = cfg.getParameter<edm::InputTag>("pfMetTag"  ); 
    tcMetTag_          = cfg.getParameter<edm::InputTag>("tcMetTag"  ); 
    chargedMetTag_     = cfg.getParameter<edm::InputTag>("chargedMetTag" ); 
    vtxTag_            = cfg.getParameter<edm::InputTag>("vtxTag"        );
    //    allCandsTag_    = cfg.getParameter<edm::InputTag>("allCandsTag"   );  // Needed for MVAMet
    chCandsTag_        = cfg.getParameter<edm::InputTag>("chCandsTag"    ); 

    if (cfg.exists("sptTag"    )) sptTag_     = cfg.getParameter<edm::InputTag>("sptTag"    ); 
    else                          sptTag_     = edm::InputTag("","","");
    if (cfg.exists("spt2Tag"   )) spt2Tag_    = cfg.getParameter<edm::InputTag>("spt2Tag"   ); 
    else                          spt2Tag_    = edm::InputTag("","","");

    produces<std::vector<reco::SkimEvent> >().setBranchAlias(cfg.getParameter<std::string>("@module_label"));

    getDYMVA_v0 = new GetDYMVA(0);
    getDYMVA_v1 = new GetDYMVA(1);

    // Needed for MVAMet
    //    fMVAMet = new MVAMet(0.1);
    //    fMVAMet->Initialize(cfg,
    //			TString((getenv("CMSSW_BASE") + std::string("/src/pharris/MVAMet/data/gbrmet_52.root"))),
    //			TString((getenv("CMSSW_BASE") + std::string("/src/pharris/MVAMet/data/gbrmetphi_52.root"))),
    //			TString((getenv("CMSSW_BASE") + std::string("/src/pharris/MVAMet/data/gbrmetu1cov_52.root"))),
    //			TString((getenv("CMSSW_BASE") + std::string("/src/pharris/MVAMet/data/gbrmetu2cov_52.root"))));
}


void SkimEventProducer2L2N::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<reco::SkimEvent> > skimEvent(new std::vector<reco::SkimEvent> );

    //SkimEvent::hypoType type = SkimEvent::hypoTypeByName(hypoType_);

    edm::Handle<reco::GenParticleCollection> genParticles;
    if(!(genParticlesTag_==edm::InputTag(""))) {
     iEvent.getByLabel(genParticlesTag_,genParticles);
    }

    edm::Handle<pat::JetCollection> jetH;
    iEvent.getByLabel(jetTag_,jetH);

    edm::Handle<double> rhoJetIso;
    edm::InputTag rho_src ("kt6PFJetsForIso","rho");
    iEvent.getByLabel(rho_src,rhoJetIso);

    edm::Handle<pat::JetCollection> tagJetH;
    if(!(tagJetTag_==edm::InputTag(""))) iEvent.getByLabel(tagJetTag_,tagJetH);

    edm::Handle<reco::PFMETCollection> pfMetH;
    iEvent.getByLabel(pfMetTag_,pfMetH);

//    edm::Handle<reco::METCollection> tcMetH;
//    iEvent.getByLabel(tcMetTag_,tcMetH);

//    edm::Handle<edm::ValueMap<reco::PFMET> > chargedMetH;
//    iEvent.getByLabel(chargedMetTag_,chargedMetH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxTag_,vtxH);

    // Needed for MVAMet
    //    edm::Handle<reco::CandidateView> allCandsH;
    //    iEvent.getByLabel(allCandsTag_, allCandsH);

    edm::Handle<reco::CandidateView> candsH;
    iEvent.getByLabel(chCandsTag_,candsH);

    edm::Handle<edm::ValueMap<float> > sptH;
    if(!(sptTag_==edm::InputTag(""))) iEvent.getByLabel(sptTag_,sptH);

    edm::Handle<edm::ValueMap<float> > spt2H;
    if(!(spt2Tag_==edm::InputTag(""))) iEvent.getByLabel(spt2Tag_,spt2H);

    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerTag_,triggerResults);

    // May God have mercy on my soul ...
    std::vector<bool> passBits;
    passBits.push_back( singleMuData_.check(iEvent,*triggerResults) );
    passBits.push_back( singleElData_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleMuData_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleElData_.check(iEvent,*triggerResults) );
    passBits.push_back( muEGData_.check(    iEvent,*triggerResults) );
    passBits.push_back( singleMuMC_.check(iEvent,*triggerResults) );
    passBits.push_back( singleElMC_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleMuMC_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleElMC_.check(iEvent,*triggerResults) );
    passBits.push_back( muEGMC_.check(    iEvent,*triggerResults) );
    passBits.push_back( AllEmbed_.check(  iEvent,*triggerResults) );


    edm::Handle<edm::View<reco::RecoCandidate> > muons;
    iEvent.getByLabel(muTag_,muons);
    edm::Handle<edm::View<reco::RecoCandidate> > softs;
    iEvent.getByLabel(softMuTag_,softs);
    edm::Handle<edm::View<reco::RecoCandidate> > electrons;
    iEvent.getByLabel(elTag_,electrons);

//    edm::Handle<double> mcGenWeight;
    edm::Handle<GenFilterInfo> mcGenWeight;
    if (!(mcGenWeightTag_==edm::InputTag(""))) {
     iEvent.getByLabel(mcGenWeightTag_, mcGenWeight);
    }

    edm::Handle<GenEventInfoProduct> GenInfoHandle;
    if (!(mcGenEventInfoTag_==edm::InputTag(""))) {
     iEvent.getByLabel(mcGenEventInfoTag_, GenInfoHandle);
    }

    // Needed for MVAMetsetJets
    //    reco::VertexCollection lVertices = *vtxH;
    //    reco::Vertex *lPV = 0;
    //    if (lVertices.size() > 0) lPV = &lVertices[0]; 
    //    
    //    //    makeCandidates(lPFInfo, allCandsH, lPV);  // Needed for MVAMet
    //    makeCandidates(lPFInfo, candsH, lPV);           // For now, use candsH instead of allCandsH
    //    makeJets      (lJetInfo, jetH, lVertices);
    //    makeVertices  (lVtxInfo, lVertices);


    //---- choice of best candidate: highest pt lepton pair

    int type1 = -1;  // 0 = electron, 1 = muon
    int type2 = -1;  // 0 = electron, 1 = muon
    int type3 = -1;  // 0 = electron, 1 = muon
    int type4 = -1;  // 0 = electron, 1 = muon
    int index1 = -1;
    int index2 = -1;
    int index3 = -1;
    int index4 = -1;
    double ptmax1 = -1;
    double ptmax2 = -1;
    double ptmax3 = -1;
    double ptmax4 = -1;
    // loop on electrons
    for(size_t i=0;i<electrons->size();++i) {
     double temppt = electrons->at(i).pt();
     if (temppt > ptmax1) {
      // what was in 3 now is in 4
      ptmax4 = ptmax3;
      type4 = type3;
      index4 = index3;

      // what was in 2 now is in 3
      ptmax3 = ptmax2;
      type3 = type2;
      index3 = index2;

      // what was in 1 now is in 2
      ptmax2 = ptmax1;
      type2 = type1;
      index2 = index1;

      // the new should be in 1
      ptmax1 = temppt;
      type1 = 0;
      index1 = i;
     }
     else if (temppt > ptmax2) {
      // what was in 3 now is in 4
      ptmax4 = ptmax3;
      type4 = type3;
      index4 = index3;

      // what was in 2 now is in 3
      ptmax3 = ptmax2;
      type3 = type2;
      index3 = index2;

      // the new should be in 1
      ptmax2 = temppt;
      type2 = 0;
      index2 = i;
     }
     else if (temppt > ptmax3) {
      // what was in 3 now is in 4
      ptmax4 = ptmax3;
      type4 = type3;
      index4 = index3;

      // the new should be in 1
      ptmax3 = temppt;
      type3 = 0;
      index3 = i;
     }
     else if (temppt > ptmax4) {

      // the new should be in 1
      ptmax4 = temppt;
      type4 = 0;
      index4 = i;
     }
    }

    // loop on muons
    for(size_t i=0;i<muons->size();++i) {
     double temppt = muons->at(i).pt();
     if (temppt > ptmax1) {
      // what was in 3 now is in 4
      ptmax4 = ptmax3;
      type4 = type3;
      index4 = index3;

      // what was in 2 now is in 3
      ptmax3 = ptmax2;
      type3 = type2;
      index3 = index2;

      // what was in 1 now is in 2
      ptmax2 = ptmax1;
      type2 = type1;
      index2 = index1;

      // the new should be in 1
      ptmax1 = temppt;
      type1 = 1;
      index1 = i;
     }
     else if (temppt > ptmax2) {
      // what was in 3 now is in 4
      ptmax4 = ptmax3;
      type4 = type3;
      index4 = index3;

      // what was in 2 now is in 3
      ptmax3 = ptmax2;
      type3 = type2;
      index3 = index2;

      // the new should be in 1
      ptmax2 = temppt;
      type2 = 1;
      index2 = i;
     }
     else if (temppt > ptmax3) {
      // what was in 3 now is in 4
      ptmax4 = ptmax3;
      type4 = type3;
      index4 = index3;

      // the new should be in 1
      ptmax3 = temppt;
      type3 = 1;
      index3 = i;
     }
     else if (temppt > ptmax4) {

      // the new should be in 1
      ptmax4 = temppt;
      type4 = 1;
      index4 = i;
     }
    }

    hypoType_ = reco::SkimEvent::hypoTypeByName("WWNONE");

    if (type1 == 0 && type2 == 0) {
     hypoType_ = reco::SkimEvent::hypoTypeByName("WWELEL");
    }
    if (type1 == 1 && type2 == 0) {
     hypoType_ = reco::SkimEvent::hypoTypeByName("WWMUEL");
    }
    if (type1 == 0 && type2 == 1) {
     hypoType_ = reco::SkimEvent::hypoTypeByName("WWELMU");
    }
    if (type1 == 1 && type2 == 1) {
     hypoType_ = reco::SkimEvent::hypoTypeByName("WWMUMU");
    }

//     std::cout << " ho scelto :: " << hypoType_ << " type1 = " << type1 << " :: " << " type2 = " << type2 << std::endl;
//     std::cout << "   electrons->size() = " << electrons->size() << std::endl;
//     std::cout << "   muons->size() =     " << muons->size() << std::endl;
//     std::cout << "   index1 = " << index1 << std::endl;
//     std::cout << "   index2 = " << index2 << std::endl;
//     std::cout << "   index3 = " << index3 << std::endl;
//     std::cout << "   index4 = " << index4 << std::endl;
    
    if(hypoType_==reco::SkimEvent::WWELEL){//ELEL
     size_t i = index1;
     size_t j = index2;
     skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );
     skimEvent->back().setEventInfo(iEvent);
                // Leptons
     skimEvent->back().setLepton(electrons,i);
     skimEvent->back().setLepton(electrons,j);


     for(size_t k=0;k<electrons->size();++k) if(k!=i && k!=j) skimEvent->back().setExtraLepton(electrons,k);
     for(size_t k=0;k<muons->size();++k) {
      float delta1 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(i).p4());
      float delta2 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(j).p4());
      if(delta1 > 0.1 && delta2 > 0.1) skimEvent->back().setExtraLepton(muons,k);
     }
     for(size_t k=0;k<softs->size();++k) skimEvent->back().setSoftMuon(softs,k);
                // Everything else
     skimEvent->back().setTriggerBits(passBits);
     skimEvent->back().setJets(jetH);
     skimEvent->back().setJetRhoIso(rhoJetIso);
     skimEvent->back().setPFMet(pfMetH);
//                skimEvent->back().setTCMet(tcMetH);
//                skimEvent->back().setChargedMet(chargedMetH->get(0));
     skimEvent->back().setVertex(vtxH);
     if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
     if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
     if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
     else                  skimEvent->back().setTagJets(jetH);
     skimEvent->back().setChargedMetSmurf(doChMET(candsH,&electrons->at(i),&electrons->at(j)));
        //      skimEvent->back().setMvaMet(getMvaMet(&electrons->at(i), &electrons->at(j), lPV, *pfMetH));
     if(genParticles.isValid()) {
      skimEvent->back().setGenParticles(genParticles);
     }
     if (!(mcGenWeightTag_==edm::InputTag(""))) {
      skimEvent->back().setGenWeight(mcGenWeight);
     }
     if (!(mcGenEventInfoTag_==edm::InputTag(""))) {
      skimEvent->back().setGenInfo(GenInfoHandle);
     }

     if (type3 == 0 && index3 != -1) skimEvent->back().setLepton(electrons,index3);
     if (type3 == 1 && index3 != -1) skimEvent->back().setLepton(muons,index3);

     if (type4 == 0 && index4 != -1) skimEvent->back().setLepton(electrons,index4);
     if (type4 == 1 && index4 != -1) skimEvent->back().setLepton(muons,index4);


    }else if(hypoType_ == reco::SkimEvent::WWELMU){
     size_t i = index1;
     size_t j = index2;

     if( muons->at(j).pt() < electrons->at(i).pt() ) { //---- only a cross check!
      skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );
      skimEvent->back().setEventInfo(iEvent);
                // Leptons
      skimEvent->back().setLepton(electrons,i);
      skimEvent->back().setLepton(muons,j);


      for(size_t k=0;k<electrons->size();++k) {
       float delta1 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(j).p4());
       if(k!=i && delta1 > 0.1) 
        skimEvent->back().setExtraLepton(electrons,k);
      }
      for(size_t k=0;k<muons->size();++k) {
       float delta1 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(i).p4());
       if(delta1 > 0.1 && k!=j) 
        skimEvent->back().setExtraLepton(muons,k);
      }
      for(size_t k=0;k<softs->size();++k) {
       if(softs->at(k).pt() != muons->at(j).pt() || softs->at(k).eta() != muons->at(j).eta()) 
        skimEvent->back().setSoftMuon(softs,k);
      }
                // Everything else
      skimEvent->back().setTriggerBits(passBits);
      skimEvent->back().setJets(jetH);
      skimEvent->back().setJetRhoIso(rhoJetIso);
      skimEvent->back().setPFMet(pfMetH);
//                skimEvent->back().setTCMet(tcMetH);
//                skimEvent->back().setChargedMet(chargedMetH->get(0));
      skimEvent->back().setVertex(vtxH);
      if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
      if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
      if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
      else                  skimEvent->back().setTagJets(jetH);
      skimEvent->back().setChargedMetSmurf(doChMET(candsH,&electrons->at(i),&muons->at(j)));
        //      skimEvent->back().setMvaMet(getMvaMet(&electrons->at(i), &muons->at(j), lPV, *pfMetH));
      if(genParticles.isValid()) {
       skimEvent->back().setGenParticles(genParticles);
      }
      if (!(mcGenWeightTag_==edm::InputTag(""))) {
       skimEvent->back().setGenWeight(mcGenWeight);
      }
      if (!(mcGenEventInfoTag_==edm::InputTag(""))) {
       skimEvent->back().setGenInfo(GenInfoHandle);
      }

      if (type3 == 0 && index3 != -1) skimEvent->back().setLepton(electrons,index3);
      if (type3 == 1 && index3 != -1) skimEvent->back().setLepton(muons,index3);

      if (type4 == 0 && index4 != -1) skimEvent->back().setLepton(electrons,index4);
      if (type4 == 1 && index4 != -1) skimEvent->back().setLepton(muons,index4);

     }

    }else if(hypoType_ == reco::SkimEvent::WWMUEL){
     size_t i = index1;
     size_t j = index2;
     if( muons->at(j).pt() >= electrons->at(i).pt() )  { //---- only a cross check!
      skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
      skimEvent->back().setEventInfo(iEvent);
                // Leptons
      skimEvent->back().setLepton(electrons,i);
      skimEvent->back().setLepton(muons,j);

      for(size_t k=0;k<electrons->size();++k) {
       float delta1 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(j).p4());
       if(k!=i && delta1 > 0.1) 
        skimEvent->back().setExtraLepton(electrons,k);
      }
      for(size_t k=0;k<muons->size();++k) {
       float delta1 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(i).p4());
       if(delta1 > 0.1 && k!=j) 
        skimEvent->back().setExtraLepton(muons,k);
      }
      for(size_t k=0;k<softs->size();++k) {
       if(softs->at(k).pt() != muons->at(j).pt() || softs->at(k).eta() != muons->at(j).eta()) 
        skimEvent->back().setSoftMuon(softs,k);
      }
                // Everything else
      skimEvent->back().setTriggerBits(passBits);
      skimEvent->back().setJets(jetH);
      skimEvent->back().setJetRhoIso(rhoJetIso);
      skimEvent->back().setPFMet(pfMetH);
//                skimEvent->back().setTCMet(tcMetH);
//                skimEvent->back().setChargedMet(chargedMetH->get(0));
      skimEvent->back().setVertex(vtxH);
      if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
      if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
      if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
      else                  skimEvent->back().setTagJets(jetH);
      skimEvent->back().setChargedMetSmurf(doChMET(candsH,&electrons->at(i),&muons->at(j)));
        //      skimEvent->back().setMvaMet(getMvaMet(&electrons->at(i), &muons->at(j), lPV, *pfMetH));
      if(genParticles.isValid()) {
       skimEvent->back().setGenParticles(genParticles);
      }
      if (!(mcGenWeightTag_==edm::InputTag(""))) {
       skimEvent->back().setGenWeight(mcGenWeight);
      }
      if (!(mcGenEventInfoTag_==edm::InputTag(""))) {
       skimEvent->back().setGenInfo(GenInfoHandle);
      }

      if (type3 == 0 && index3 != -1) skimEvent->back().setLepton(electrons,index3);
      if (type3 == 1 && index3 != -1) skimEvent->back().setLepton(muons,index3);

      if (type4 == 0 && index4 != -1) skimEvent->back().setLepton(electrons,index4);
      if (type4 == 1 && index4 != -1) skimEvent->back().setLepton(muons,index4);
     }

    }else if(hypoType_==reco::SkimEvent::WWMUMU){//MUMU
     size_t i = index1;
     size_t j = index2;
     skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
     skimEvent->back().setEventInfo(iEvent);
                // Leptons
     skimEvent->back().setLepton(muons,i);
     skimEvent->back().setLepton(muons,j);

     for(size_t k=0;k<electrons->size();++k) {
      float delta1 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(i).p4());
      float delta2 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(j).p4());
      if(delta1 > 0.1 && delta2 > 0.1) 
       skimEvent->back().setExtraLepton(electrons,k);
     }
     for(size_t k=0;k<muons->size();++k) {
      if(k!=i && k!=j) 
       skimEvent->back().setExtraLepton(muons,k);
     }
     for(size_t k=0;k<softs->size();++k) {
      if( (softs->at(k).pt() != muons->at(i).pt() || softs->at(k).eta() != muons->at(i).eta())  && 
           (softs->at(k).pt() != muons->at(j).pt() || softs->at(k).eta() != muons->at(j).eta()) ) 
       skimEvent->back().setSoftMuon(softs,k);
     }
                // Everything else
     skimEvent->back().setTriggerBits(passBits);
     skimEvent->back().setJets(jetH);
     skimEvent->back().setJetRhoIso(rhoJetIso);
     skimEvent->back().setPFMet(pfMetH);
//                skimEvent->back().setTCMet(tcMetH);
//                skimEvent->back().setChargedMet(chargedMetH->get(0));
     skimEvent->back().setVertex(vtxH);
     if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
     if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
     if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
     else                  skimEvent->back().setTagJets(jetH);
     skimEvent->back().setChargedMetSmurf(doChMET(candsH,&muons->at(i),&muons->at(j)));
        //      skimEvent->back().setMvaMet(getMvaMet(&muons->at(i), &muons->at(j), lPV, *pfMetH));
     if(genParticles.isValid()) {
      skimEvent->back().setGenParticles(genParticles);
     }
     if (!(mcGenWeightTag_==edm::InputTag(""))) {
      skimEvent->back().setGenWeight(mcGenWeight);
     }
     if (!(mcGenEventInfoTag_==edm::InputTag(""))) {
      skimEvent->back().setGenInfo(GenInfoHandle);
     }


//      std::cout << " fin qui tutto bene! " << std::endl;
     if (type3 == 0 && index3 != -1) skimEvent->back().setLepton(electrons,index3);
     if (type3 == 1 && index3 != -1) skimEvent->back().setLepton(muons,index3);

     if (type4 == 0 && index4 != -1) skimEvent->back().setLepton(electrons,index4);
     if (type4 == 1 && index4 != -1) skimEvent->back().setLepton(muons,index4);

    }else{
//         throw cms::Exception("BadInput") 
//             << "ERROR: event type " << hypoType_ << " is not known" ;
//      std::cout << "ERROR: event type " << hypoType_ << " is not known" << std::endl ;
    }


    for (size_t jevent=0; jevent<skimEvent->size(); jevent++) {

      reco::SkimEvent* event = &skimEvent->at(jevent);

      addDYMVA(event);
    }


    iEvent.put(skimEvent);
}


reco::MET SkimEventProducer2L2N::doChMET(edm::Handle<reco::CandidateView> candsH,
        const reco::Candidate* cand1,const reco::Candidate* cand2){
    using namespace std;
    reco::Candidate::LorentzVector totalP4;
    for(reco::CandidateView::const_iterator it= candsH->begin(), ed =candsH->end(); it != ed; ++it){
        if( it->charge() == 0 ) continue;
        if(fabs(ROOT::Math::VectorUtil::DeltaR(it->p4(),cand1->p4())) <=0.1) continue;
        if(fabs(ROOT::Math::VectorUtil::DeltaR(it->p4(),cand2->p4())) <=0.1) continue;
        totalP4 += it->p4();
    }
    totalP4 +=cand1->p4();
    totalP4 +=cand2->p4();
    reco::Candidate::LorentzVector invertedP4(-totalP4);
    reco::MET met(invertedP4,reco::Candidate::Point(0,0,0));
    return met;
}


void SkimEventProducer2L2N::addDYMVA(reco::SkimEvent* event)
{
  float dymva0 = -999;
  float dymva1 = -999;

  if (event->nLep() >= 2) {

    size_t index           = 0;
    float  minPt           = 0;
    float  eta             = 4.7;
    int    applyCorrection = 1;
    int    applyID         = 4;

    float jet1pt  = event->leadingJetPt (index, minPt, eta, applyCorrection, applyID);
    float jet1phi = event->leadingJetPhi(index, minPt, eta, applyCorrection, applyID);

    double dPhiDiLepJet1 = fabs(event->dPhillLeadingJet(eta, applyCorrection, applyID));
    double dPhiJet1MET = fabs(deltaPhi(jet1phi, event->pfMetPhi()));
    double dPhillPfMET = fabs(event->dPhillPfMet());

    if (jet1pt < 15) {
      jet1pt        =   15;
      dPhiDiLepJet1 = -0.1;
      dPhiJet1MET   = -0.1;
    }

    float  px_rec = event->pfMet()*cos(event->pfMetPhi()) + event->pXll();       
    float  py_rec = event->pfMet()*sin(event->pfMetPhi()) + event->pYll();
    double recoil = sqrt(px_rec*px_rec + py_rec*py_rec);

    dymva0 = getDYMVA_v0->getValue(event->nCentralJets(30.0, eta, applyCorrection,applyID),
				   event->pfMet(),
				   event->chargedMetSmurf(),
				   jet1pt,
				   event->pfMetSignificance(),
				   dPhiDiLepJet1,
				   dPhiJet1MET,
				   event->mTHiggs(event->PFMET));

    dymva1 = getDYMVA_v1->getValue(event->nCentralJets(30.0, eta, applyCorrection, applyID),
				   event->projPfMet(),
				   event->projChargedMetSmurf(),
				   event->nGoodVertices(),
				   event->pTll(),
				   jet1pt,
				   event->pfMetMEtSig(),
				   dPhiDiLepJet1,
				   dPhillPfMET,
				   dPhiJet1MET,
				   recoil,
				   event->mTHiggs(event->PFMET));
  
  }

  event->addUserFloat("dymva0", dymva0);
  event->addUserFloat("dymva1", dymva1);
}


//------------------------------------------------------------------------------
// makeJets
//------------------------------------------------------------------------------
void SkimEventProducer2L2N::makeJets(std::vector<MetUtilities::JetInfo>    &iJetInfo,
				 const edm::Handle<pat::JetCollection> &jH,
				 reco::VertexCollection                &iVertices)
{
  iJetInfo.clear();

  pat::JetRefVector jrv;

  jrv.clear();

  for (size_t ijet=0; ijet<jH->size(); ++ijet) {

    jrv.push_back(pat::JetRef(jH,ijet));

    double eta                 = jrv[ijet]->correctedJet("Uncorrected", "none").eta();
    double energy              = jrv[ijet]->correctedJet("Uncorrected", "none").energy();
    double neutralHadronEnergy = jrv[ijet]->correctedJet("Uncorrected", "none").neutralHadronEnergy();
    double neutralEmEnergy     = jrv[ijet]->correctedJet("Uncorrected", "none").neutralEmEnergy();
    double chargedHadronEnergy = jrv[ijet]->correctedJet("Uncorrected", "none").chargedHadronEnergy();
    double chargedEmEnergy     = jrv[ijet]->correctedJet("Uncorrected", "none").chargedEmEnergy();
    //    int    nConstituents       = jrv[ijet]->correctedJet("Uncorrected", "none").nConstituents();
    int    chargedMultiplicity = jrv[ijet]->correctedJet("Uncorrected", "none").chargedMultiplicity();

    if (energy == 0)                                             continue;
    if (neutralHadronEnergy / energy > 0.99)                     continue; 
    if (neutralEmEnergy     / energy > 0.99)                     continue;
    //    if (nConstituents < 2)                                       continue;
    if (chargedHadronEnergy / energy <= 0    && fabs(eta) < 2.4) continue;
    if (chargedEmEnergy     / energy >  0.99 && fabs(eta) < 2.4) continue;
    if (chargedMultiplicity < 1              && fabs(eta) < 2.4) continue;

    double lNeuFrac = (energy > 0) ? (neutralEmEnergy + neutralHadronEnergy) / energy : 9999;

    MetUtilities::JetInfo pJetObject; 

    pJetObject.p4       = jrv[ijet]->p4(); 
    pJetObject.mva      = jrv[ijet]->userFloat("jetMva");
    pJetObject.neutFrac = lNeuFrac;

    iJetInfo.push_back(pJetObject);
  }
}


//------------------------------------------------------------------------------
// makeCandidates
//------------------------------------------------------------------------------
void SkimEventProducer2L2N::makeCandidates(std::vector<std::pair<LorentzVector,double> > &iPFInfo,
				       edm::Handle<reco::CandidateView>               cH,
				       reco::Vertex                                  *iPV)
{ 
  iPFInfo.clear();

  for (reco::CandidateView::const_iterator it=cH->begin(), ed=cH->end(); it!=ed; ++it) {

    double pDZ = -999;

    if (iPV != 0) {

      double bsx = iPV->x();
      double bsy = iPV->y();
      double bsz = iPV->z();

      double vx = it->vx();
      double vy = it->vy();
      double vz = it->vz();

      if (vx != 0 || vy != 0 || vz != 0) {

	double px = it->p4().px();
	double py = it->p4().py();
	double pz = it->p4().pz();
	double pt = it->p4().pt();
	
	pDZ = fabs((vz - bsz) - ((vx - bsx)*px + (vy - bsy)*py)/pt * pz/pt);

	if (pDZ == 0) pDZ = -999;
      }
    }

    std::pair<LorentzVector,double> pPFObject(it->p4(), pDZ);

    iPFInfo.push_back(pPFObject);
  }
}


//------------------------------------------------------------------------------
// makeVertices
//------------------------------------------------------------------------------
void SkimEventProducer2L2N::makeVertices(std::vector<Vector>    &iPVInfo,
				     reco::VertexCollection &iVertices)
{ 
  iPVInfo.clear();

  for (int i0 = 0; i0<(int)iVertices.size(); i0++) {

    const reco::Vertex *pVertex = &(iVertices.at(i0));

    Vector pVec;

    pVec.SetCoordinates(pVertex->x(), pVertex->y(), pVertex->z());

    iPVInfo.push_back(pVec);
  }
}


//------------------------------------------------------------------------------
// getMvaMet
//------------------------------------------------------------------------------
reco::PFMET SkimEventProducer2L2N::getMvaMet(const reco::Candidate *cand1,
					 const reco::Candidate *cand2,
					 reco::Vertex          *iPV,
					 reco::PFMETCollection  thePfMet)
{
  LorentzVector lVis1 = cand1->p4();
  LorentzVector lVis2 = cand2->p4();

  std::vector<LorentzVector > theLeptons;

  theLeptons.push_back(lVis1);
  theLeptons.push_back(lVis2);

  std::pair<LorentzVector,TMatrixD> lMVAMetInfo = fMVAMet->GetMet(theLeptons,
								  lJetInfo,
								  lPFInfo,
								  lVtxInfo,
								  false);

  reco::PFMET lDummy;

  reco::PFMET lMVAMet(lDummy.getSpecific(),
		      thePfMet.at(0).sumEt(),
		      lMVAMetInfo.first,
		      iPV->position());
  
  return lMVAMet;
}


SkimEventProducer2L2N::~SkimEventProducer2L2N() {
  delete getDYMVA_v0;
  delete getDYMVA_v1;
} 
void SkimEventProducer2L2N::beginJob() { }
void SkimEventProducer2L2N::endJob() { }
DEFINE_FWK_MODULE(SkimEventProducer2L2N);
