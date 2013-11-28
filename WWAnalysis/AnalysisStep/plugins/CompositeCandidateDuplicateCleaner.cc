#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include <cmath>
#include <algorithm>

namespace { const double MZ = 91.188; }

class CompositeCandidateDuplicateCleaner : public edm::EDProducer {
    public:
        explicit CompositeCandidateDuplicateCleaner(const edm::ParameterSet&);
        ~CompositeCandidateDuplicateCleaner();

    private:
        edm::InputTag src_;

        enum Algo { BestZ1Mass };
        Algo algo_;

        // interface
        virtual void produce(edm::Event&, const edm::EventSetup&);

        /// pair that uniquely identifies each candidate
        typedef std::pair<edm::ProductID, unsigned int> CandId;

        ///  get id from shallow clone or ptr clone
        CandId candId(const reco::Candidate &c) const {
            if (c.hasMasterClonePtr()) return CandId(c.masterClonePtr().id(), c.masterClonePtr().key());
            if (c.hasMasterClone()) return CandId(c.masterClone().id(), c.masterClone().key());
            return CandId();
        }

        bool sameLeptons(const reco::CompositeCandidate &c1, const reco::CompositeCandidate &c2) const ;
        bool firstIsBest(const reco::CompositeCandidate &first, const reco::CompositeCandidate &second) const ;

};




CompositeCandidateDuplicateCleaner::CompositeCandidateDuplicateCleaner(const edm::ParameterSet& iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
    std::string algo = iConfig.getParameter<std::string>("algo");

    if (algo == "BestZ1Mass") algo_ = BestZ1Mass;
    else throw cms::Exception("Configuration") << "Unknown algo '" << algo << "'. Supported values are: 'BestZ1Mass'.\n";

    produces<pat::CompositeCandidateCollection>();
}

CompositeCandidateDuplicateCleaner::~CompositeCandidateDuplicateCleaner() {}


void CompositeCandidateDuplicateCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace std;

    // Read Candidates from Event
    edm::Handle<edm::View<reco::CompositeCandidate> > src;
    iEvent.getByLabel(src_,src);

    // Prepare output
    auto_ptr<pat::CompositeCandidateCollection> out(new pat::CompositeCandidateCollection());

    // Prepare scoreboard
    std::vector<bool> good(src->size(), true);

    // Loop
    for (unsigned int i = 0, n = src->size(); i < n; ++i) {
        const reco::CompositeCandidate &c1 = (*src)[i];
        // skip those already flagged as bad
        if (!good[i]) continue;
        // search for dups among followers
        for (unsigned int j = i+1; j < n; ++j) {
            const reco::CompositeCandidate &c2 = (*src)[j];
            if (sameLeptons(c1,c2)) {
                // competing
                if (firstIsBest(c1,c2)) {
                    good[j] = false; // we win
                } else {
                    good[i] = false; // we loose
                    break;
                }
            }
        }
        // abort if we've been beaten by one of our followers
        if (!good[i]) continue;

        // OK, we're good to go

        // Make a copy, possibly converting to pat::CompositeCandidate on the fly
        if (typeid(c1) == typeid(pat::CompositeCandidate)) {
            out->push_back(dynamic_cast<const pat::CompositeCandidate&>(c1));
        } else {
            out->push_back(pat::CompositeCandidate(c1));
        }
    }

    iEvent.put(out);
}

bool CompositeCandidateDuplicateCleaner::sameLeptons(const reco::CompositeCandidate &c1, const reco::CompositeCandidate &c2) const {
    if (c1.numberOfDaughters() == 2 && c2.numberOfDaughters() == 2) {
        const reco::Candidate *c1_z1 = c1.daughter(0), *c1_z2 = c1.daughter(1);
        const reco::Candidate *c2_z1 = c2.daughter(0), *c2_z2 = c2.daughter(1);
        if (c1_z1->numberOfDaughters() == 2 && c1_z2->numberOfDaughters() == 2 &&
            c2_z1->numberOfDaughters() == 2 && c2_z2->numberOfDaughters() == 2) {
            CandId leps1[4], leps2[4];
            for (int i = 0; i < 4; ++i) {
                leps1[i] = candId(*c1.daughter(i/2)->daughter(i%2));
                leps2[i] = candId(*c2.daughter(i/2)->daughter(i%2));
            }
            std::sort(&leps1[0], &leps1[4]);
            std::sort(&leps2[0], &leps2[4]);
            for (int i = 0; i < 4; ++i) {
                if (leps1[i].first != leps2[i].first || leps1[i].second != leps2[i].second) {
                    return false;
                }
            }
            return true;
        } else {
            throw cms::Exception("NotSupported") << "sameLeptons only works for ZZ events for now\n";
        }
    } else throw cms::Exception("NotSupported") << "sameLeptons requires the candidates to have two daughters each\n";
}

bool CompositeCandidateDuplicateCleaner::firstIsBest(const reco::CompositeCandidate &c1, const reco::CompositeCandidate &c2) const {
    if (c1.numberOfDaughters() == 2 && c2.numberOfDaughters() == 2) {
        const reco::Candidate *c1_z1 = c1.daughter(0), *c1_z2 = c1.daughter(1);
        const reco::Candidate *c2_z1 = c2.daughter(0), *c2_z2 = c2.daughter(1);
        if (c1_z1->numberOfDaughters() == 2 && c1_z2->numberOfDaughters() == 2 &&
            c2_z1->numberOfDaughters() == 2 && c2_z2->numberOfDaughters() == 2) {
            switch (algo_) {
                case BestZ1Mass:
                    if (fabs(c1_z1->mass()-MZ) > fabs(c1_z2->mass()-MZ)) std::swap(c1_z1, c1_z2);
                    if (fabs(c2_z1->mass()-MZ) > fabs(c2_z2->mass()-MZ)) std::swap(c2_z1, c2_z2);
                    return fabs(c1_z1->mass()-MZ) <= fabs(c2_z1->mass()-MZ);
            }
            return true;
        } else {
            throw cms::Exception("NotSupported") << "firstIsBest only works for ZZ events for now\n";
        }
    } else throw cms::Exception("NotSupported") << "firstIsBest re= fabs(c2_z1->mass()-MZ);uires the candidates to have two daughters each\n";
}


DEFINE_FWK_MODULE(CompositeCandidateDuplicateCleaner);
