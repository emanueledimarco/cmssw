// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "WWAnalysis/AnalysisStep/interface/MySingleDeposit.h"

#include<vector>

class PatMuonIsoAdder : public edm::EDProducer {
    public:
        explicit PatMuonIsoAdder(const edm::ParameterSet&);
        ~PatMuonIsoAdder();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag muonTag_;
        std::vector<MySingleDeposit> sources_;

};

PatMuonIsoAdder::PatMuonIsoAdder(const edm::ParameterSet& iConfig) :
    muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")) {

    std::vector<edm::ParameterSet> depPSets = iConfig.getParameter<std::vector<edm::ParameterSet> >("deposits");
    for (std::vector<edm::ParameterSet>::const_iterator it = depPSets.begin(), ed = depPSets.end(); it != ed; ++it) {
        sources_.push_back(MySingleDeposit(*it));
    }

    produces<pat::MuonCollection>();  
}


PatMuonIsoAdder::~PatMuonIsoAdder() {

    std::vector<MySingleDeposit>::iterator it, begin = sources_.begin(), end = sources_.end();
    for (it = begin; it != end; ++it) it->cleanup();

}


// ------------ method called to produce the data  ------------
void PatMuonIsoAdder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<reco::Candidate> > muons;
    iEvent.getByLabel(muonTag_,muons);

    std::vector<MySingleDeposit>::iterator it, begin = sources_.begin(), end = sources_.end();
    for (it = begin; it != end; ++it) it->open(iEvent, iSetup);

    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);

    // ----- here is the real loop over the muons ----
    for(edm::View<reco::Candidate>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){    
        pat::Muon clone = *edm::RefToBase<reco::Candidate>(muons,mu-muons->begin()).castTo<pat::MuonRef>();
        if(clone.type() == 8) continue;
        const reco::CandidateBaseRef musRef2(muons,mu-muons->begin());
        for (it = begin; it != end; ++it) clone.addUserFloat(it->getLabel(),it->compute(musRef2)); 

        pOut->push_back(clone);
    }
    iEvent.put(pOut);
}

void PatMuonIsoAdder::beginJob() { }
void PatMuonIsoAdder::endJob() { } 
DEFINE_FWK_MODULE(PatMuonIsoAdder);
