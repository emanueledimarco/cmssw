// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Electron.h>
#include "WWAnalysis/AnalysisStep/interface/MySingleDeposit.h"

#include<vector>

class PatElectronIsoAdder : public edm::EDProducer {
    public:
        explicit PatElectronIsoAdder(const edm::ParameterSet&);
        ~PatElectronIsoAdder();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag eleTag_;
        std::vector<MySingleDeposit> sources_;

};

PatElectronIsoAdder::PatElectronIsoAdder(const edm::ParameterSet& iConfig) :
    eleTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")) {

    std::vector<edm::ParameterSet> depPSets = iConfig.getParameter<std::vector<edm::ParameterSet> >("deposits");
    for (std::vector<edm::ParameterSet>::const_iterator it = depPSets.begin(), ed = depPSets.end(); it != ed; ++it) {
        sources_.push_back(MySingleDeposit(*it));
    }

    produces<pat::ElectronCollection>();  
}


PatElectronIsoAdder::~PatElectronIsoAdder() {

    std::vector<MySingleDeposit>::iterator it, begin = sources_.begin(), end = sources_.end();
    for (it = begin; it != end; ++it) it->cleanup();

}


// ------------ method called to produce the data  ------------
void PatElectronIsoAdder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<reco::Candidate> > eles;
    iEvent.getByLabel(eleTag_,eles);

    std::vector<MySingleDeposit>::iterator it, begin = sources_.begin(), end = sources_.end();
    for (it = begin; it != end; ++it) it->open(iEvent, iSetup);

    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    // ----- here is the real loop over the eles ----
    for(edm::View<reco::Candidate>::const_iterator el=eles->begin(); el!=eles->end(); ++el){    
        pat::Electron clone = *edm::RefToBase<reco::Candidate>(eles,el-eles->begin()).castTo<pat::ElectronRef>();
        const reco::CandidateBaseRef elesRef2(eles,el-eles->begin());
        for (it = begin; it != end; ++it) {
            clone.addUserFloat(it->getLabel(),it->compute(elesRef2)); 
        }

        pOut->push_back(clone);
    }
    iEvent.put(pOut);
}

void PatElectronIsoAdder::beginJob() { }
void PatElectronIsoAdder::endJob() { } 
DEFINE_FWK_MODULE(PatElectronIsoAdder);
