#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "WWAnalysis/AnalysisStep/interface/sorters.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

class BestZCandidateChooser : public edm::EDProducer {
    public:
        explicit BestZCandidateChooser(const edm::ParameterSet&);
        ~BestZCandidateChooser();

    private:
        edm::InputTag src;
        bool passesIsolation(const pat::Muon*, const pat::Muon*) const;
        bool passesIsolation(const pat::Electron*, const pat::Electron*) const;
        virtual void produce(edm::Event&, const edm::EventSetup&);
};




BestZCandidateChooser::BestZCandidateChooser(const edm::ParameterSet& iConfig) :
    src(iConfig.getParameter<edm::InputTag>("src"))
{
    produces<std::vector<reco::SkimEvent2L> >();
}

BestZCandidateChooser::~BestZCandidateChooser() {}


void BestZCandidateChooser::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<std::vector<reco::SkimEvent2L> > zcandH;
    iEvent.getByLabel(src,zcandH);
    std::vector<reco::SkimEvent2L> zcand = *zcandH;

    std::auto_ptr<std::vector<reco::SkimEvent2L> > out(new std::vector<reco::SkimEvent2L>());
        
    sorters::ZMassSorter zsorter;
    std::sort(zcand.begin(), zcand.end(), zsorter);
    if (zcand.size()>0) out->push_back(zcand[0]);

    iEvent.put(out);
}

DEFINE_FWK_MODULE(BestZCandidateChooser);
