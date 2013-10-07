#include <vector>
#include <algorithm>
#include <cmath>

#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// bunch of anonymous tools
namespace {
    typedef std::pair<edm::ProductID,size_t> edmID;
    template<typename T, typename C, typename F>  edmID id(const edm::Ref<T,C,F> &r) { return edmID(r.id(), r.key()); }
    template<typename T> edmID id(const edm::Ptr<T> &r) { return edmID(r.id(), r.key()); }
    template<typename T> edmID id(const edm::RefToBase<T> &r) { return edmID(r.id(), r.key()); }
    edmID id(const reco::Candidate *c) { return c->hasMasterClonePtr() ? id(c->masterClonePtr()) : id(c->masterClone()); }
    
    struct SortIndexByPt {
        public:
            SortIndexByPt(const edm::View<pat::PFParticle> &pho, bool byPt) : pho_(pho), byPt_(byPt) {}
            bool operator()(int i, int j) const { 
                if (byPt_ || (pho_[i].pt() > 4 || pho_[j].pt() > 4)) {
                    return pho_[i].pt() > pho_[j].pt(); 
                } else {
                    return deltaR(pho_[i], *pho_[i].overlaps("goodLepNoIso")[0]) < deltaR(pho_[j], *pho_[j].overlaps("goodLepNoIso")[0]);
                }
            }
        private:
            const edm::View<pat::PFParticle> &pho_;
            bool  byPt_;
    };
}

class FsrCollectorBase : public edm::EDProducer {
    public:
        FsrCollectorBase(const edm::ParameterSet &iConfig) ;
    protected:
        edm::InputTag photons_;
        StringCutObjectSelector<pat::PFParticle> photonSelection_;
        std::string photonMatch_;
        std::string isolationLabel_;
        double      isolationCut_;
        bool doFsrRecovery_;
        bool sortPhotonsByPt_;

        void associatePhotonsToLeptons(const edm::View<pat::PFParticle> &photons, std::map<edmID,std::vector<int> > &leptonToPhotonMap) ;
        double getRelIsoFSR(const reco::Candidate &lep, const pat::PFParticle &pho) const ;
        double getRelIsoNoFSR(const reco::Candidate &lep) const ;
        double getRelIsoFSR(const pat::Muon &lep, const pat::PFParticle &pho) const ;
        double getRelIsoNoFSR(const pat::Muon &lep) const ;
        double getRelIsoFSR(const pat::Electron &lep, const pat::PFParticle &pho) const ;
        double getRelIsoNoFSR(const pat::Electron &lep) const ;
};

FsrCollectorBase::FsrCollectorBase(const edm::ParameterSet &iConfig) :
    photons_(iConfig.getParameter<edm::InputTag>("photons")),
    photonSelection_(iConfig.getParameter<std::string>("photonSelection")),
    photonMatch_(iConfig.getParameter<std::string>("photonMatch")),
    isolationLabel_(iConfig.getParameter<std::string>("isolationLabel")),
    isolationCut_(iConfig.getParameter<double>("isolationCut")),
    doFsrRecovery_(iConfig.getParameter<bool>("doFsrRecovery")),
    sortPhotonsByPt_(iConfig.getParameter<bool>("sortPhotonsByPt"))
{
}

void FsrCollectorBase::associatePhotonsToLeptons(const edm::View<pat::PFParticle> &photons, 
                                                 std::map<edmID,std::vector<int> > &leptonToPhotonMap) {
    int idx = 0;
    for (edm::View<pat::PFParticle>::const_iterator it = photons.begin(), ed = photons.end(); it != ed; ++it, ++idx) {
        if (!photonSelection_(*it)) continue;
        const reco::CandidatePtrVector &leps = it->overlaps(photonMatch_);
        if (leps.empty()) continue;
        leptonToPhotonMap[id(leps[0])].push_back(idx);
    }       
}

class SkimEvent2LFsrCollector : public FsrCollectorBase {
    public:
        SkimEvent2LFsrCollector(const edm::ParameterSet &iConfig) ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag zll_;
        std::string   label_;
};

class SingleLeptonFsrCollector : public FsrCollectorBase {
    public:
        SingleLeptonFsrCollector(const edm::ParameterSet &iConfig) ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag leptons_;
        std::string   label_;
};

SkimEvent2LFsrCollector::SkimEvent2LFsrCollector(const edm::ParameterSet &iConfig) :
    FsrCollectorBase(iConfig),
    zll_(iConfig.getParameter<edm::InputTag>("zll")),
    label_(iConfig.existsAs<std::string>("label") ? iConfig.getParameter<std::string>("label") : std::string())
{
    produces<std::vector<reco::SkimEvent2L> >();
}

SingleLeptonFsrCollector::SingleLeptonFsrCollector(const edm::ParameterSet &iConfig) :
    FsrCollectorBase(iConfig),
    leptons_(iConfig.getParameter<edm::InputTag>("leptons")),
    label_(iConfig.getParameter<std::string>("label"))
{
    produces<edm::OwnVector<reco::Candidate> >();
}


void
SkimEvent2LFsrCollector::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<std::vector<reco::SkimEvent2L> > zll; 
    iEvent.getByLabel(zll_, zll);
    edm::Handle<edm::View<pat::PFParticle> > photons; 
    if (doFsrRecovery_) iEvent.getByLabel(photons_, photons);

    std::auto_ptr<std::vector<reco::SkimEvent2L> > out(new std::vector<reco::SkimEvent2L>());

    // Associate to each lepton a list of photons which have that lepton as closest one
    std::map<edmID,std::vector<int> > leptonToPhotonMap;
    if (doFsrRecovery_) associatePhotonsToLeptons(*photons, leptonToPhotonMap);
 
    // Now loop on Z candidates
    foreach(const reco::SkimEvent2L &z0, *zll) {
        // Retrieve FSR photons, make a list
        std::vector<int> myphotons;
        if (doFsrRecovery_)  {
            foreach(const int &i, leptonToPhotonMap[id(z0.daughter(0))]) myphotons.push_back(i);
            foreach(const int &i, leptonToPhotonMap[id(z0.daughter(1))]) myphotons.push_back(i);
            // now sort them by pt
            std::sort(myphotons.begin(), myphotons.end(), SortIndexByPt(*photons, sortPhotonsByPt_));
        }
        // now loop and pick the first one that works, if any
        const pat::PFParticle *match = 0; int imatch = -1;
        foreach(const int &i, myphotons) {
            const pat::PFParticle *pho = &photons->at(i);
            // Kinematic cuts
            reco::Particle::LorentzVector p4ll = z0.p4(), p4llg = p4ll + pho->p4();
            double mll = p4ll.M(), mllg = p4llg.M();
            if (mllg > 100 || mllg <  4) continue;
            if (fabs(mllg - 91.188) >= fabs(mll - 91.188)) continue;
            match = pho; imatch = i; break;
        }

        // Now test isolation
        double l1iso = getRelIsoNoFSR(z0.l(0)), l2iso = getRelIsoNoFSR(z0.l(1));

        if (match) {
            // must figure out if this was a photon created on the fly from a muon or not
            // since those were not in the isolation sum in the first place
            double minDR = std::min<double>(deltaR(z0.l(0), *match), deltaR(z0.l(1), *match));
            // if minDR is nearly zero, we can ignore this photon, since it was not in the isolation sums to begin with
            if (minDR > 3e-4) { 
                l1iso = getRelIsoFSR(z0.l(0), *match);
                l2iso = getRelIsoFSR(z0.l(1), *match);
            }
        }

        if (l1iso < isolationCut_ && l2iso < isolationCut_) {

            out->push_back(z0);
            if (match) {
                reco::ShallowClonePtrCandidate clone(photons->ptrAt(imatch), 0, match->p4(), match->vertex());
                clone.setPdgId(22); clone.setStatus(match->status());
                out->back().addDaughter(clone);
                out->back().setP4(out->back().p4() + clone.p4());
            }
            if (!label_.empty()) {
                out->back().addUserFloat(label_ + "[0]", l1iso);
                out->back().addUserFloat(label_ + "[1]", l2iso);
            }
        }
    }

    iEvent.put(out);
}

void
SingleLeptonFsrCollector::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<reco::Candidate> > leptons; 
    iEvent.getByLabel(leptons_, leptons);
    edm::Handle<edm::View<pat::PFParticle> > photons; 
    if (doFsrRecovery_) iEvent.getByLabel(photons_, photons);

    std::auto_ptr<edm::OwnVector<reco::Candidate> > out(new edm::OwnVector<reco::Candidate>());

    // Associate to each lepton a list of photons which have that lepton as closest one
    std::map<edmID,std::vector<int> > leptonToPhotonMap; 
    if (doFsrRecovery_) associatePhotonsToLeptons(*photons, leptonToPhotonMap);
 
    // Now loop on lepton candidates
    for (unsigned int i = 0, n = leptons->size(); i < n; ++i) {
        reco::CandidateBaseRef lep = leptons->refAt(i);
        out->push_back(*lep);
        const pat::PFParticle *match = 0;
        if (doFsrRecovery_)  {
            std::vector<int> myphotons;
            foreach(const int &j, leptonToPhotonMap[id(lep)]) myphotons.push_back(j);
            if (!myphotons.empty()) {
                std::sort(myphotons.begin(), myphotons.end(), SortIndexByPt(*photons, sortPhotonsByPt_));
                match = & (*photons)[myphotons.front()];
            }
        }
        double iso = (match ? getRelIsoFSR(*lep, *match) : getRelIsoNoFSR(*lep));
        if (abs(lep->pdgId()) == 13) {
            pat::Muon *mu = dynamic_cast<pat::Muon *>(&out->back());
            if (mu == 0) throw cms::Exception("CorruptData") << "Muon is not a pat::Muon but a " << typeid(*lep).name() << "\n";
            mu->addUserFloat(label_, iso);
        } else if (abs(lep->pdgId()) == 11) {
            pat::Electron *el = dynamic_cast<pat::Electron *>(&out->back());
            if (el == 0) throw cms::Exception("CorruptData") << "Electron is not a pat::Electron but a " << typeid(*lep).name() << "\n";
            el->addUserFloat(label_, iso);
        }
    }

    iEvent.put(out);
}

double FsrCollectorBase::getRelIsoFSR(const reco::Candidate &lep, const pat::PFParticle &pho) const {
    return abs(lep.pdgId()) == 13 ? 
                getRelIsoFSR(dynamic_cast<const pat::Muon     &>(lep), pho) :
                getRelIsoFSR(dynamic_cast<const pat::Electron &>(lep), pho);
}
double FsrCollectorBase::getRelIsoNoFSR(const reco::Candidate &lep) const {
    return abs(lep.pdgId()) == 13 ? 
                getRelIsoNoFSR(dynamic_cast<const pat::Muon     &>(lep)) :
                getRelIsoNoFSR(dynamic_cast<const pat::Electron &>(lep));
}
double FsrCollectorBase::getRelIsoFSR(const pat::Muon &lep, const pat::PFParticle &pho) const {
    double dr = deltaR(lep,pho);
    if (dr < 0.4 && dr >= 0.01 && pho.pt() > 0.5) {
        float chiso = lep.userFloat("pfChHadIso04");
        float nhiso = lep.userFloat("pfNHadIso04");
        float phiso = lep.userFloat("pfPhotonIso04");
        float chpu  = lep.userFloat("pfChHadPU04");
        if (phiso >= pho.pt()) phiso -= pho.pt();
        double iso = chiso + std::max<double>(0., nhiso + phiso - 0.5*chpu);
        
//         std::cout << "SkimEvent2LFSRCollector::getRelIsoFSR(muons) : "
//                   << lep.pt() << " " << lep.eta() << " " << lep.phi() << " | "
//                   << chiso << " " << nhiso << " " << phiso << " " << chpu << " : " << iso << " " << iso/lep.pt() << "\n";
        
        return iso/lep.pt();
    } else {
        return getRelIsoNoFSR(lep); 
    }

}
double FsrCollectorBase::getRelIsoNoFSR(const pat::Muon &lep) const {
    float chiso = lep.userFloat("pfChHadIso04");
    float nhiso = lep.userFloat("pfNHadIso04");
    float phiso = lep.userFloat("pfPhotonIso04");
    float chpu  = lep.userFloat("pfChHadPU04");
    double iso = chiso + std::max<double>(0., nhiso + phiso - 0.5*chpu);

//     std::cout << "SkimEvent2LFSRCollector::getRelIsoNoFSR(muons) : "
//               << lep.pt() << " " << lep.eta() << " " << lep.phi() << " | "
//               << chiso << " " << nhiso << " " << phiso << " " << chpu << " : " << iso << " " << iso/lep.pt() << "\n";
    
    return iso/lep.pt();
}

double FsrCollectorBase::getRelIsoFSR(const pat::Electron &lep, const pat::PFParticle &pho) const {
    double dr = deltaR(lep,pho);
    if (dr < 0.4 && (fabs(lep.superCluster()->eta()) < 1.479 || dr > 0.08)) {
        float chiso = lep.userFloat(isolationLabel_+"ChHad");
        float nhiso = lep.userFloat(isolationLabel_+"NHad");
        float phiso = lep.userFloat(isolationLabel_+"Pho");
        float eatot = lep.userFloat(isolationLabel_+"EAtot");
        float rho   = lep.userFloat(isolationLabel_+"Rho");
        //if (phiso < pho.pt()) throw cms::Exception("LogicError") << "Cannot subtract photon of pt " << pho.pt() << " from photon iso " << phiso << ", for electron of pt " << lep.pt() << ", eta " << lep.eta() << ", phi " << lep.phi() << ", deltaR = " << dr << std::endl;
        if (phiso >= pho.pt()) phiso -= pho.pt();
        double iso = chiso + std::max<double>(0., nhiso + phiso - rho*eatot);

        return iso/lep.pt();
    } else {
        return getRelIsoNoFSR(lep); 
    }
}
double FsrCollectorBase::getRelIsoNoFSR(const pat::Electron &lep) const {
    return lep.userFloat(isolationLabel_)/lep.pt();
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent2LFsrCollector);
DEFINE_FWK_MODULE(SingleLeptonFsrCollector);
