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

class BestHiggsCandidateChooser : public edm::EDProducer {
    public:
        explicit BestHiggsCandidateChooser(const edm::ParameterSet&);
        ~BestHiggsCandidateChooser();

    private:
        edm::InputTag src;
        bool passesIsolation(const pat::Muon*, const pat::Muon*) const;
        bool passesIsolation(const pat::Electron*, const pat::Electron*) const;
        virtual void produce(edm::Event&, const edm::EventSetup&);
};




BestHiggsCandidateChooser::BestHiggsCandidateChooser(const edm::ParameterSet& iConfig) :
    src(iConfig.getParameter<edm::InputTag>("src"))
{
    produces<std::vector<reco::SkimEvent4L> >();
}

BestHiggsCandidateChooser::~BestHiggsCandidateChooser() {}


void BestHiggsCandidateChooser::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<std::vector<reco::SkimEvent4L> > higgsH;
    iEvent.getByLabel(src,higgsH);
    std::vector<reco::SkimEvent4L> higgs = *higgsH;

    std::auto_ptr<std::vector<reco::SkimEvent4L> > out(new std::vector<reco::SkimEvent4L>());
        
    sorters::HiggsZLeptonPtSorter hsorter;
    std::sort(higgs.begin(), higgs.end(), hsorter);
    if (higgs.size()>0) out->push_back(higgs[0]);

    iEvent.put(out);
}

DEFINE_FWK_MODULE(BestHiggsCandidateChooser);
