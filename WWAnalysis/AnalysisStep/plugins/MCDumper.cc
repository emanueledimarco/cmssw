#include "WWAnalysis/AnalysisStep/interface/MCDumper.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"



MCDumper::MCDumper(const edm::ParameterSet& cfg) {

    GenJetTag_ =  cfg.getParameter<edm::InputTag>("GenJetTag");
    GenMetTag_ =  cfg.getParameter<edm::InputTag>("GenMetTag");

    fs_ = edm::Service<TFileService>().operator->();
    tree_ = fs_->make<TTree>(cfg.getParameter<std::string>("@module_label").c_str(),cfg.getParameter<std::string>("@module_label").c_str());
    tree_->Branch("MHiggs"     ,&MHiggs_     ,"MHiggs/F" );
    tree_->Branch("PtHiggs"    ,&PtHiggs_    ,"PtHiggs/F" );
    tree_->Branch("EtaHiggs"   ,&EtaHiggs_   ,"EtaHiggs/F" );
    tree_->Branch("PhiHiggs"   ,&PhiHiggs_   ,"PhiHiggs/F" );
    tree_->Branch("PdgIdHiggs" ,&PdgIdHiggs_ ,"PdgIdHiggs/F" );

    tree_->Branch("V1pt"    ,&V1pt_    ,"V1pt/F" );
    tree_->Branch("V1eta"   ,&V1eta_   ,"V1eta/F" );
    tree_->Branch("V1phi"   ,&V1phi_   ,"V1phi/F" );
    tree_->Branch("V1pdgId" ,&V1pdgId_ ,"V1pdgId/F" );
    tree_->Branch("V2pt"    ,&V2pt_    ,"V2pt/F" );
    tree_->Branch("V2eta"   ,&V2eta_   ,"V2eta/F" );
    tree_->Branch("V2phi"   ,&V2phi_   ,"V2phi/F" );
    tree_->Branch("V2pdgId" ,&V2pdgId_ ,"V2pdgId/F" );

    tree_->Branch("f1V1pt"    ,&f1V1pt_    ,"f1V1pt/F" );
    tree_->Branch("f1V1eta"   ,&f1V1eta_   ,"f1V1eta/F" );
    tree_->Branch("f1V1phi"   ,&f1V1phi_   ,"f1V1phi/F" );
    tree_->Branch("f1V1pdgId" ,&f1V1pdgId_ ,"f1V1pdgId/F" );
    tree_->Branch("f2V1pt"    ,&f2V1pt_    ,"f2V1pt/F" );
    tree_->Branch("f2V1eta"   ,&f2V1eta_   ,"f2V1eta/F" );
    tree_->Branch("f2V1phi"   ,&f2V1phi_   ,"f2V1phi/F" );
    tree_->Branch("f2V1pdgId" ,&f2V1pdgId_ ,"f2V1pdgId/F" );

    tree_->Branch("f2V2pt"    ,&f2V2pt_    ,"f2V2pt/F" );
    tree_->Branch("f2V2eta"   ,&f2V2eta_   ,"f2V2eta/F" );
    tree_->Branch("f2V2phi"   ,&f2V2phi_   ,"f2V2phi/F" );
    tree_->Branch("f2V2pdgId" ,&f2V2pdgId_ ,"f2V2pdgId/F" );
    tree_->Branch("f1V2pt"    ,&f1V2pt_    ,"f1V2pt/F" );
    tree_->Branch("f1V2eta"   ,&f1V2eta_   ,"f1V2eta/F" );
    tree_->Branch("f1V2phi"   ,&f1V2phi_   ,"f1V2phi/F" );
    tree_->Branch("f1V2pdgId" ,&f1V2pdgId_ ,"f1V2pdgId/F" );

    tree_->Branch("jet1pt"    ,&jet1pt_    ,"jet1pt/F" );
    tree_->Branch("jet1eta"   ,&jet1eta_   ,"jet1eta/F" );
    tree_->Branch("jet1phi"   ,&jet1phi_   ,"jet1phi/F" );

    tree_->Branch("jet2pt"    ,&jet2pt_    ,"jet2pt/F" );
    tree_->Branch("jet2eta"   ,&jet2eta_   ,"jet2eta/F" );
    tree_->Branch("jet2phi"   ,&jet2phi_   ,"jet2phi/F" );

    tree_->Branch("jet3pt"    ,&jet3pt_    ,"jet3pt/F" );
    tree_->Branch("jet3eta"   ,&jet3eta_   ,"jet3eta/F" );
    tree_->Branch("jet3phi"   ,&jet3phi_   ,"jet3phi/F" );

    tree_->Branch("jet4pt"    ,&jet4pt_    ,"jet4pt/F" );
    tree_->Branch("jet4eta"   ,&jet4eta_   ,"jet4eta/F" );
    tree_->Branch("jet4phi"   ,&jet4phi_   ,"jet4phi/F" );

    tree_->Branch("jet5pt"    ,&jet5pt_    ,"jet5pt/F" );
    tree_->Branch("jet5eta"   ,&jet5eta_   ,"jet5eta/F" );
    tree_->Branch("jet5phi"   ,&jet5phi_   ,"jet5phi/F" );

    tree_->Branch("metpt"    ,&metpt_    ,"metpt/F" );
    tree_->Branch("metphi"   ,&metphi_   ,"metphi/F" );

}


void MCDumper::FindDaughterParticles(const reco::Candidate** pCurrent, std::vector<const reco::Candidate*>* pFinal) const {

  // Variable for decayed current particle
 const reco::Candidate* pCurrentNew = 0;

 for(unsigned int i = 0; i < (*pCurrent) -> numberOfDaughters(); ++i) {
  if( (pFinal != 0) && ((*pCurrent) -> daughter(i) -> pdgId() != (*pCurrent) -> pdgId()) )
   pFinal -> push_back((*pCurrent) -> daughter(i));

  if((*pCurrent) -> daughter(i) -> pdgId() == (*pCurrent) -> pdgId())
   pCurrentNew = (*pCurrent) -> daughter(i);
 }

  // Change the current particle into the decayed one
 if(pCurrentNew)
  (*pCurrent) = pCurrentNew;
}



void MCDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::GenJet> > genJetHandle;
    iEvent.getByLabel(GenJetTag_, genJetHandle);
    edm::View<reco::GenJet> genJets = *genJetHandle;

    edm::Handle<edm::View<reco::GenMET> > genMetHandle;
    iEvent.getByLabel(GenMetTag_, genMetHandle);
    edm::View<reco::GenMET> genMet = *genMetHandle;

    jet1pt_ = -100;
    jet1eta_= -100;
    jet1phi_= -100;

    jet2pt_ = -100;
    jet2eta_= -100;
    jet2phi_= -100;

    jet3pt_ = -100;
    jet3eta_= -100;
    jet3phi_= -100;

    jet4pt_ = -100;
    jet4eta_= -100;
    jet4phi_= -100;

    jet5pt_ = -100;
    jet5eta_= -100;
    jet5phi_= -100;

    // loop on jets
    for(unsigned int genJetIt = 0; genJetIt < genJets.size(); ++genJetIt)
    {
     reco::GenJet genJet = genJets.at(genJetIt);

     if (genJetIt == 0) {
      jet1pt_ = genJet.pt();
      jet1eta_= genJet.eta();
      jet1phi_= genJet.phi();
     }

     if (genJetIt == 1) {
      jet2pt_ = genJet.pt();
      jet2eta_= genJet.eta();
      jet2phi_= genJet.phi();
     }

     if (genJetIt == 2) {
      jet3pt_ = genJet.pt();
      jet3eta_= genJet.eta();
      jet3phi_= genJet.phi();
     }

     if (genJetIt == 3) {
      jet4pt_ = genJet.pt();
      jet4eta_= genJet.eta();
      jet4phi_= genJet.phi();
     }

     if (genJetIt == 4) {
      jet5pt_ = genJet.pt();
      jet5eta_= genJet.eta();
      jet5phi_= genJet.phi();
     }

    }



    //---- met ----
    metpt_ = -100;
    metphi_= -100;

    metpt_  = (genMet.at(0)).pt();
    metphi_ = (genMet.at(0)).phi();


    //---- Higgs ----

    edm::Handle<reco::GenParticleCollection> genParticlesGlobal;
    iEvent.getByLabel("genParticles",genParticlesGlobal);


    reco::GenParticleRefVector genParticles;
    genParticles.clear();
    for(size_t i=0; i<genParticlesGlobal->size(); ++i) {
     genParticles.push_back( reco::GenParticleRef(genParticlesGlobal,i) );
    }


    float mass = -100;
    float pt = -100;
    float eta = -100;
    float phi = -100;
    float HpdgId = -100;

    float V1pt = -100;
    float V1eta = -100;
    float V1phi = -100;
    float V1pdgId = -100;
    float V2pt = -100;
    float V2eta = -100;
    float V2phi = -100;
    float V2pdgId = -100;

    float f1V1pt = -100;
    float f1V1eta = -100;
    float f1V1phi = -100;
    float f1V1pdgId = -100;
    float f2V1pt = -100;
    float f2V1eta = -100;
    float f2V1phi = -100;
    float f2V1pdgId = -100;

    float f1V2pt = -100;
    float f1V2eta = -100;
    float f1V2phi = -100;
    float f1V2pdgId = -100;
    float f2V2pt = -100;
    float f2V2eta = -100;
    float f2V2phi = -100;
    float f2V2pdgId = -100;

    const reco::Candidate* mcH = 0;

  // loop over gen particles
    for(size_t gp=0; gp<genParticles.size();++gp){

     int pdgId  = genParticles[gp] -> pdgId();
     int status = genParticles[gp] -> status();

    //     Higgs  {25}       graviton {5000039}
     if( ((pdgId == 25) || (pdgId == 5000039)) && (status == 3) ) {
      mcH = &(*(genParticles[gp]));
      mass = mcH->mass();
      pt = mcH->pt();
      eta = mcH->eta();
      phi = mcH->phi();
      HpdgId = pdgId;

     std::vector<const reco::Candidate*> VfromHBuffer;
     FindDaughterParticles(&mcH,&VfromHBuffer);

     // H > VV
     if(VfromHBuffer.size() == 2) {
       const reco::Candidate* mcV1;
       const reco::Candidate* mcV2;

       mcV1 = VfromHBuffer.at(0);
       mcV2 = VfromHBuffer.at(1);

       V1pt    = mcV1->pt();
       V1eta   = mcV1->eta();
       V1phi   = mcV1->phi();
       V1pdgId = mcV1->pdgId();

       V2pt    = mcV2->pt();
       V2eta   = mcV2->eta();
       V2phi   = mcV2->phi();
       V2pdgId = mcV2->pdgId();

       const reco::Candidate* mcF1_fromV1;
       const reco::Candidate* mcF2_fromV1;
       const reco::Candidate* mcF1_fromV2;
       const reco::Candidate* mcF2_fromV2;

       bool isHWWok = true;
       std::vector<const reco::Candidate*> fFromV1Buffer;
       FindDaughterParticles(&mcV1,&fFromV1Buffer);
       if(fFromV1Buffer.size() == 2) {
        mcF1_fromV1 = fFromV1Buffer.at(0);
        mcF2_fromV1 = fFromV1Buffer.at(1);
        // If leptons, see if there is a photon emission
        if(abs(mcF1_fromV1 -> pdgId()) >= 11) {
         FindDaughterParticles(&mcF1_fromV1);
        }
        if(abs(mcF2_fromV1 -> pdgId()) >= 11) {
         FindDaughterParticles(&mcF2_fromV1);
        }
       }
       else {
        isHWWok  = false;
       }

       std::vector<const reco::Candidate*> fFromV2Buffer;
       FindDaughterParticles(&mcV2,&fFromV2Buffer);
       if(fFromV2Buffer.size() == 2) {
        mcF1_fromV2 = fFromV2Buffer.at(0);
        mcF2_fromV2 = fFromV2Buffer.at(1);
        // If leptons, see if there is a photon emission
        if(abs(mcF1_fromV2 -> pdgId()) >= 11) {
         FindDaughterParticles(&mcF1_fromV2);
        }
        if(abs(mcF2_fromV2 -> pdgId()) >= 11) {
         FindDaughterParticles(&mcF2_fromV2);
        }
       }
       else {
        isHWWok  = false;
       }

       if (isHWWok) {
        f1V1pt = mcF1_fromV1->pt();
        f1V1eta = mcF1_fromV1->eta();
        f1V1phi = mcF1_fromV1->phi();
        f1V1pdgId = mcF1_fromV1->pdgId();

        f2V1pt = mcF2_fromV1->pt();
        f2V1eta = mcF2_fromV1->eta();
        f2V1phi = mcF2_fromV1->phi();
        f2V1pdgId = mcF2_fromV1->pdgId();

        f1V2pt = mcF1_fromV2->pt();
        f1V2eta = mcF1_fromV2->eta();
        f1V2phi = mcF1_fromV2->phi();
        f1V2pdgId = mcF1_fromV2->pdgId();

        f2V2pt = mcF2_fromV2->pt();
        f2V2eta = mcF2_fromV2->eta();
        f2V2phi = mcF2_fromV2->phi();
        f2V2pdgId = mcF2_fromV2->pdgId();

       }
//      if ( (abs(mcF1_fromV1 -> pdgId()) == 13 || abs(mcF1_fromV1 -> pdgId()) == 14) && (abs(mcF1_fromV2 -> pdgId()) == 13 || abs(mcF1_fromV2 -> pdgId()) == 14) ) {
      }
     }
    } // loop over gen particles

    MHiggs_= mass;
    PtHiggs_ = pt;
    EtaHiggs_= eta;
    PhiHiggs_= phi;
    PdgIdHiggs_= HpdgId;

    V1pt_    = V1pt;
    V1eta_   = V1eta;
    V1phi_   = V1phi;
    V1pdgId_ = V1pdgId;

    V2pt_    = V2pt;
    V2eta_   = V2eta;
    V2phi_   = V2phi;
    V2pdgId_ = V2pdgId;

    f1V1pt_    = f1V1pt;
    f1V1eta_   = f1V1eta;
    f1V1phi_   = f1V1phi;
    f1V1pdgId_ = f1V1pdgId;
    f2V1pt_    = f2V1pt;
    f2V1eta_   = f2V1eta;
    f2V1phi_   = f2V1phi;
    f2V1pdgId_ = f2V1pdgId;

    f1V2pt_    = f1V2pt;
    f1V2eta_   = f1V2eta;
    f1V2phi_   = f1V2phi;
    f1V2pdgId_ = f1V2pdgId;
    f2V2pt_    = f2V2pt;
    f2V2eta_   = f2V2eta;
    f2V2phi_   = f2V2phi;
    f2V2pdgId_ = f2V2pdgId;

    tree_->Fill();

}

void MCDumper::beginJob() { }
MCDumper::~MCDumper() { } 
void MCDumper::endJob() { }
DEFINE_FWK_MODULE(MCDumper);


