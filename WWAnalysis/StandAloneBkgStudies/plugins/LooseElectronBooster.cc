// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Electron.h>


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include<vector>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


class LooseElectronBooster : public edm::EDProducer {
    public:
        explicit LooseElectronBooster(const edm::ParameterSet&);
        ~LooseElectronBooster();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

  edm::InputTag eleTag_;
  edm::InputTag elePairTag_;
  edm::InputTag jetTag_;
  edm::InputTag leadingElTag_;
  edm::InputTag metTag_;
  edm::InputTag triggerTag_;
  double drMin_;


  TriggerBitChecker* triggerEl8;
  TriggerBitChecker* triggerEl17;

  StringCutObjectSelector<reco::Candidate>* selectorEleTight;
  StringCutObjectSelector<reco::Candidate>* selectorEleLoose;

};

LooseElectronBooster::LooseElectronBooster(const edm::ParameterSet& iConfig) :
  eleTag_(iConfig.getParameter<edm::InputTag>("eleTag")),
  elePairTag_(iConfig.getParameter<edm::InputTag>("elePairTag")),
  jetTag_(iConfig.getParameter<edm::InputTag>("jetTag")),
  leadingElTag_(iConfig.getParameter<edm::InputTag>("leadingElTag")),
  metTag_(iConfig.getParameter<edm::InputTag>("metTag")),
  triggerTag_(iConfig.getParameter<edm::InputTag>("triggerTag")),
  drMin_(iConfig.getParameter<double>("drMin"))
{
  produces<pat::ElectronCollection>();  
  triggerEl8 = new TriggerBitChecker("HLT_Ele8_CaloIdL_CaloIsoVL_v*");
  triggerEl17 = new TriggerBitChecker("HLT_Ele17_CaloIdL_CaloIsoVL_v*");


  selectorEleTight = new StringCutObjectSelector<reco::Candidate>(iConfig.getParameter<std::string>("tightEleSelection"),true);
  selectorEleLoose = new StringCutObjectSelector<reco::Candidate>(iConfig.getParameter<std::string>("looseEleSelection"),true);

}


LooseElectronBooster::~LooseElectronBooster() {
  delete triggerEl8;
  delete triggerEl17;
}


// ------------ method called to produce the data  ------------
void LooseElectronBooster::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;


    //ESHandle<TransientTrackBuilder> theTTBuilder;
    //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    Handle<View<Candidate> > eles;
    iEvent.getByLabel(eleTag_,eles);

    Handle<View<Candidate> > elePairs;
    iEvent.getByLabel(elePairTag_,elePairs);

    Handle<View<Candidate> > jets;
    iEvent.getByLabel(jetTag_,jets);

    Handle<View<Candidate> > leadingEl;
    iEvent.getByLabel(leadingElTag_,leadingEl);

    Handle<PFMETCollection> pfmet;
    iEvent.getByLabel(metTag_,pfmet);
    
    Handle<TriggerResults> triggerResults;
    iEvent.getByLabel(triggerTag_,triggerResults);


    PFMET myMET= pfmet->front();


    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    double massPair(-999.);
    if(elePairs->size()) massPair = elePairs->begin()->mass();


    // ----- here is the real loop over the eles ----
    for(edm::View<reco::Candidate>::const_iterator el=eles->begin(); el!=eles->end(); ++el){    
      pat::Electron clone = *edm::RefToBase<reco::Candidate>(eles,el-eles->begin()).castTo<pat::ElectronRef>();

      clone.addUserFloat("pfmet",myMET.pt());
      float mt = sqrt(2*clone.pt()*myMET.pt()*(1-cos(clone.phi()-myMET.phi())));
      clone.addUserFloat("mt",mt);
      
      bool passLooseSel = (*selectorEleLoose)(*el);
      bool passTightSel = (*selectorEleTight)(*el);

      clone.addUserFloat("passTight",passTightSel);
      clone.addUserFloat("passLoose",passLooseSel);

      clone.addUserFloat("pairMass",massPair);

      double leadingElPt(-1);
      if(leadingEl->size())
	leadingElPt = leadingEl->front().pt();
      clone.addUserFloat("leadingLepPt",leadingElPt);
      

      clone.addUserInt("passEl8",triggerEl8->check(iEvent,*triggerResults));
      clone.addUserInt("passEl17",triggerEl17->check(iEvent,*triggerResults));

      if(jets->size()==0)  {
	//cout << "no jet for this event" << endl;
	pOut->push_back(clone);
      }
      else{
	double dr = ROOT::Math::VectorUtil::DeltaR(jets->begin()->momentum(),clone.momentum());
	if(dr>drMin_) pOut->push_back(clone);
      }
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void LooseElectronBooster::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void LooseElectronBooster::endJob() { } 


//define this as a plug-in
DEFINE_FWK_MODULE(LooseElectronBooster);



