#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h" // needed, as IsoDeposit.h does not include it
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"

#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


class CompositeCandidateIsoAdder : public edm::EDProducer {
    public:
        explicit CompositeCandidateIsoAdder(const edm::ParameterSet&);
        ~CompositeCandidateIsoAdder();

    private:
        // interface
        virtual void produce(edm::Event&, const edm::EventSetup&);

        // definitions
        typedef std::pair<reco::isodeposit::Direction, int> CandDirAndId;

        /// A class that can veto around a list of candidates, with pdg-id dependent candidate vetos
        /// We make this thing non-copyable, since it contains a vector of pointers
        struct CompCandVeto : public reco::isodeposit::AbsVeto, boost::noncopyable {
            public:
                CompCandVeto(const edm::ParameterSet &) ;
                ~CompCandVeto() ;
                // initialization per event
                void setCands(const std::vector<CandDirAndId> * cands) { candsAndPdgId_ = cands; }
                // base class interface
                virtual void centerOn(double eta, double phi) {}
                virtual bool veto(double eta, double phi, float value) const ;
            private: 
                // data
                mutable std::map<int, reco::isodeposit::AbsVetos>  vetosForPdgId_;
                const std::vector<CandDirAndId> *                  candsAndPdgId_;
                // utilities
                void makeVetos(const std::vector<std::string> &strVetos, reco::isodeposit::AbsVetos &out);
        };

        /// A class that can compute the isolation for a single particle type
        class SingleObjectIsolation {
            public:
                /// Dummy default constructor for use with std::map
                SingleObjectIsolation() {}

                /// create from a cfg 
                SingleObjectIsolation(const edm::ParameterSet &cfg) ;

                /// read data from the event
                void setup(const edm::Event &iEvent) { iEvent.getByLabel(src_, map_); }

                /// set the candidates to veto around
                void setCands(const std::vector<CandDirAndId> * cands) { veto_->setCands(cands); }

                /// compute isolation for a single candidate 
                double isolate(const reco::Candidate &cand) const {
                    reco::isodeposit::AbsVetos vetos(1, veto_.get()); // only one veto, but very powerful one....
                    const reco::IsoDeposit & dep = (cand.hasMasterClonePtr() ? (*map_)[cand.masterClonePtr()] : (*map_)[cand.masterClone()]);
                    return dep.depositWithin(radius_, vetos, true);
                }
            private:
                edm::InputTag src_;
                edm::Handle<reco::IsoDepositMap> map_;
                double radius_;
                boost::shared_ptr<CompCandVeto> veto_; // shared ptr, so there are no issues with copies
                                    
        };

        /// A structure that contains a set of isolations for each particle type
        struct IsolationConfig : public std::map<int, SingleObjectIsolation> { 
            /// create from a cfg file
            IsolationConfig(const edm::ParameterSet &) ;
            /// label of the output
            std::string label;
            /// read data from the event
            void setup(const edm::Event &iEvent) {
                foreach (value_type &p, *this) { 
                    p.second.setup(iEvent); 
                }
            }
            /// set the candidates to veto around
            void setCands(const std::vector<CandDirAndId> * cands) { 
                foreach (value_type &p, *this) { 
                    p.second.setCands(cands); 
                }
            }
            
        };

        //=========== DATA ========
        /// Composite Candidates
        edm::InputTag src_;
        std::vector<IsolationConfig> configs_;

        //=========== UTILITIES ========
        void recursiveAddDaughters(const reco::Candidate *c, std::vector<CandDirAndId> &out) ;
        void recursiveDoIso(const reco::Candidate *c, const IsolationConfig &cfg, const char *postfix, pat::CompositeCandidate &out) ;
};


/// ============= CompositeCandidateIsoAdder::CompCandVeto ==================================
CompositeCandidateIsoAdder::CompCandVeto::CompCandVeto(const edm::ParameterSet &iConfig) :
        candsAndPdgId_(0) {
    makeVetos(iConfig.getParameter<std::vector<std::string> >("muonVetos"),     vetosForPdgId_[13]);
    makeVetos(iConfig.getParameter<std::vector<std::string> >("electronVetos"), vetosForPdgId_[11]);
    makeVetos(iConfig.getParameter<std::vector<std::string> >("photonVetos"),   vetosForPdgId_[22]);
}

CompositeCandidateIsoAdder::CompCandVeto::~CompCandVeto() {
    typedef std::map<int, reco::isodeposit::AbsVetos>::value_type idVetoPair;
    foreach(idVetoPair &vid, vetosForPdgId_) {
        foreach(reco::isodeposit::AbsVeto *v2, vid.second) delete v2;
    }
}

void CompositeCandidateIsoAdder::CompCandVeto::makeVetos(const std::vector<std::string> &vetos, reco::isodeposit::AbsVetos &out) {
    reco::isodeposit::EventDependentAbsVeto *evdep = 0; 
    foreach(const std::string &v, vetos) {
        out.push_back(IsoDepositVetoFactory::make(v.c_str(), evdep));
        if (evdep != 0) throw cms::Exception("Configuration") << "Used event-dependent veto " << v << "\n";
    }
}


bool CompositeCandidateIsoAdder::CompCandVeto::veto(double eta, double phi, float value) const {
    foreach(const CandDirAndId &cand, *candsAndPdgId_) {
        foreach(reco::isodeposit::AbsVeto *v, vetosForPdgId_[cand.second]) {
            v->centerOn(cand.first.eta(), cand.first.phi());
            if (v->veto(eta,phi,value)) return true;
        }
    }
    return false;
}

// ===================== CompositeCandidateIsoAdder::SingleObjectIsolation ============================
CompositeCandidateIsoAdder::SingleObjectIsolation::SingleObjectIsolation(const edm::ParameterSet& pset) :
    src_(pset.getParameter<edm::InputTag>("src")),
    radius_(pset.getParameter<double>("deltaR")),
    veto_(new CompCandVeto(pset)) {
}

// ===================== CompositeCandidateIsoAdder::SingleObjectIsolation ============================
CompositeCandidateIsoAdder::IsolationConfig::IsolationConfig(const edm::ParameterSet& cfg) :
    label(cfg.getParameter<std::string>("label"))
{
    const int   np = 3;
    const char *name[np] = { "muon", "electron", "photon" };
    const int   pid[np]  = {   13  ,     11    ,    22    };
    for (int i = 0; i < np; ++i) {
        std::string psetname = name[i] + std::string("Isolation");
        if (cfg.existsAs<edm::ParameterSet>(psetname)) {
            (*this)[pid[i]] = SingleObjectIsolation(cfg.getParameter<edm::ParameterSet>(psetname));
        }
    }
}



CompositeCandidateIsoAdder::CompositeCandidateIsoAdder(const edm::ParameterSet& iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
    produces<pat::CompositeCandidateCollection>();
    std::vector<edm::ParameterSet> psets = iConfig.getParameter<std::vector<edm::ParameterSet> >("deposits");
    for (unsigned int i = 0; i < psets.size(); ++i) {
        configs_.push_back(IsolationConfig(psets[i]));
    }
}

CompositeCandidateIsoAdder::~CompositeCandidateIsoAdder() {}


void CompositeCandidateIsoAdder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace std;

    // Read Candidates from Event
    edm::Handle<edm::View<reco::CompositeCandidate> > srcs;
    iEvent.getByLabel(src_,srcs);

    // Read IsoDeposit maps from Event
    foreach(IsolationConfig &ic, configs_) {
        ic.setup(iEvent);
    }

    // Prepare output
    auto_ptr<pat::CompositeCandidateCollection> out(new pat::CompositeCandidateCollection());
    foreach(const reco::CompositeCandidate &src, *srcs)  {
        // Make a copy, possibly converting to pat::CompositeCandidate on the fly
        if (typeid(src) == typeid(pat::CompositeCandidate)) {
            out->push_back(dynamic_cast<const pat::CompositeCandidate&>(src));
        } else {
            out->push_back(pat::CompositeCandidate(src));
        }

        // Get the copy, so we can add user data to it
        pat::CompositeCandidate &cc = out->back();

        // Figure out daughters for vetos
        vector<CandDirAndId> daughters;
        recursiveAddDaughters(&cc, daughters);

        // For each isolation that we want to compute
        foreach(IsolationConfig &ic, configs_) {

            // notify vetos of these daughters
            ic.setCands(&daughters);

            // now go compute isolation for the daughters, and add it as user data
            recursiveDoIso(&cc, ic, "", cc); 
        }
    }

    iEvent.put(out);
}

void CompositeCandidateIsoAdder::recursiveAddDaughters(const reco::Candidate *c, std::vector<CandDirAndId> &out) {
    for (unsigned int i = 0; i < c->numberOfDaughters(); ++i) {
        recursiveAddDaughters(c->daughter(i), out);
    }
    out.push_back(CandDirAndId(reco::isodeposit::Direction(c->eta(), c->phi()), abs(c->pdgId())));
}

void CompositeCandidateIsoAdder::recursiveDoIso(const reco::Candidate *c, const IsolationConfig &cfg, const char *postfix, pat::CompositeCandidate &out) {
    int id = abs(c->pdgId());
    IsolationConfig::const_iterator it = cfg.find(id);
    if (it != cfg.end()) {
        out.addUserFloat(cfg.label+postfix, it->second.isolate(*c));
    } 
    for (unsigned int i = 0; i < c->numberOfDaughters(); ++i) {
        char newpostfix[255]; sprintf(newpostfix,"%s_%d", postfix, i);
        recursiveDoIso(c->daughter(i), cfg, newpostfix, out);;
    }
}

DEFINE_FWK_MODULE(CompositeCandidateIsoAdder);
