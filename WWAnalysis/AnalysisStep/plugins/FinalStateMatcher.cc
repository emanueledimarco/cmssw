#include <vector>
#include <algorithm>
#include <cmath>

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"


class RecoFinalStateMatcher : public edm::EDFilter {
    public:
        RecoFinalStateMatcher(const edm::ParameterSet &iConfig) ;
        virtual bool filter(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        struct GenLeptonWithFsr {
            reco::GenParticleRef initial;
            reco::GenParticleRef final;
            reco::GenParticleRefVector photons;
        };
        struct MatchStruct {
            MatchStruct(double aDist, int genIdx, int genStat, int recIdx):
                dist(aDist), gen(genIdx), rec(recIdx), genStatus(genStat) {}
                bool operator<(const MatchStruct &other) const { return dist < other.dist; }
                bool overlap(const MatchStruct &other) const { return gen == other.gen || rec == other.rec; }   
                bool isDropped() const { return dist < 0; }
                void drop() { dist = -1; }
            public:
                double dist;
                int    gen, rec; 
                int    genStatus; // mc status (to know if we matched to status == 3 or status == 1)
        };

        edm::InputTag genParticles_;
        double genPhotonPtMin_;

        edm::InputTag electrons_;
        StringCutObjectSelector<reco::Candidate> electronSelection_;
        double electronPtScale_, electronDetaScale_, electronDphiScale_;

        edm::InputTag muons_;
        StringCutObjectSelector<reco::Candidate>     muonSelection_;
        double muonPtScale_, muonDetaScale_, muonDphiScale_;

        edm::InputTag photons_;
        StringCutObjectSelector<reco::Candidate>     photonSelection_;
        double photonPtScale_, photonDetaScale_, photonDphiScale_;

        bool debug_;

        bool fillGenLepton(GenLeptonWithFsr &lepton);
        void match(const reco::GenParticle   &gen,  int i,
                   const reco::CandidateView &reco, const StringCutObjectSelector<reco::Candidate> &isgood, 
                   double ptScale, double detaScale, double dphiScale,
                   std::vector<MatchStruct> & out) ;
        void resolve(std::vector<MatchStruct>& out) ;
        void fillZBoson(const std::vector<GenLeptonWithFsr> &leptons, int il1, std::map<int,reco::CandidatePtr> &genToReco, reco::CompositeCandidate &cc) ;
};

RecoFinalStateMatcher::RecoFinalStateMatcher(const edm::ParameterSet &iConfig) :
    genParticles_(iConfig.getParameter<edm::InputTag>("genParticles")),
    genPhotonPtMin_(iConfig.getParameter<double>("genPhotonPtMin")),
    electrons_(iConfig.getParameter<edm::InputTag>("electrons")),
    electronSelection_(iConfig.getParameter<std::string>("electronSelection"), true),
    electronPtScale_(iConfig.getParameter<double>("electronDPtRelScale")),
    electronDetaScale_(iConfig.getParameter<double>("electronDEtaScale")),
    electronDphiScale_(iConfig.getParameter<double>("electronDPhiScale")),
    muons_(iConfig.getParameter<edm::InputTag>("muons")),
    muonSelection_(iConfig.getParameter<std::string>("muonSelection"), true),
    muonPtScale_(iConfig.getParameter<double>("muonDPtRelScale")),
    muonDetaScale_(iConfig.getParameter<double>("muonDEtaScale")),
    muonDphiScale_(iConfig.getParameter<double>("muonDPhiScale")),
    photons_(iConfig.getParameter<edm::InputTag>("photons")),
    photonSelection_(iConfig.getParameter<std::string>("photonSelection"), true),
    photonPtScale_(iConfig.getParameter<double>("photonDPtRelScale")),
    photonDetaScale_(iConfig.getParameter<double>("photonDEtaScale")),
    photonDphiScale_(iConfig.getParameter<double>("photonDPhiScale")),
    debug_(iConfig.getUntrackedParameter<bool>("debug", false))
{
    produces<reco::CompositeCandidateCollection>();
    produces<edm::Association<reco::GenParticleCollection> >();
}

bool
RecoFinalStateMatcher::filter(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<reco::GenParticleCollection> genParticles; 
    iEvent.getByLabel(genParticles_, genParticles);

    edm::Handle<reco::CandidateView> electrons; 
    edm::Handle<reco::CandidateView> muons; 
    edm::Handle<reco::CandidateView> photons; 
    iEvent.getByLabel(electrons_, electrons);
    iEvent.getByLabel(muons_, muons);
    iEvent.getByLabel(photons_, photons);
    if (debug_) std::cout << "Got electrons " << electrons_.encode()  << " from product id " << electrons.id() << std::endl;
    if (debug_) std::cout << "Got muons " << muons_.encode()  << " from product id " << muons.id() << std::endl;
    if (debug_) std::cout << "Got photons " << photons_.encode()  << " from product id " << photons.id() << std::endl;

    std::vector<int> electronMatches(electrons->size(), -1);
    std::vector<int> muonMatches(muons->size(), -1);
    std::vector<int> photonMatches(photons->size(), -1);

    std::vector<reco::GenParticleRef> zs;
    for (int i = 0, n = genParticles->size(); i < n; ++i) {
        const reco::GenParticle & gen = (*genParticles)[i];
        if (gen.pdgId() == 23 && gen.status() == 3 && gen.numberOfDaughters() > 0 && 
            (abs(gen.daughter(0)->pdgId()) == 11 || abs(gen.daughter(0)->pdgId()) == 13)) {
            zs.push_back(reco::GenParticleRef(genParticles, i));
        }
    }
    if (debug_) std::cout << "Found " << zs.size() << " Z bosons decaying to electrons or muons" << std::endl;
    if (zs.size() != 2) return false;
    if (zs[0]->mass() > zs[1]->mass()) std::swap(zs[1], zs[0]);

    std::vector<GenLeptonWithFsr> genLeptons(4);
    for (int i = 0; i < 4; ++i) {
        genLeptons[i].initial = zs[i/2]->daughterRef(i%2);
        if (!fillGenLepton(genLeptons[i])) return false;
    }
   
    // now we produce list of matches
    std::vector<MatchStruct> electronCandMatches, muonCandMatches, photonCandMatches;
    for (int i = 0; i < 4; ++i) {
        if (abs(genLeptons[i].initial->pdgId()) == 13) {
            if (debug_) std::cout << "Matching gen lepton " << i << " to muons: " << std::endl;
            match(*genLeptons[i].initial, genLeptons[i].final.key(), *muons, muonSelection_, muonPtScale_, muonDetaScale_, muonDphiScale_, muonCandMatches);
            match(*genLeptons[i].final,   genLeptons[i].final.key(), *muons, muonSelection_, muonPtScale_, muonDetaScale_, muonDphiScale_, muonCandMatches);
        } else {
            if (debug_) std::cout << "Matching gen lepton " << i << " to electrons: " << std::endl;
            match(*genLeptons[i].initial, genLeptons[i].final.key(), *electrons, electronSelection_, electronPtScale_, electronDetaScale_, electronDphiScale_, electronCandMatches);
            match(*genLeptons[i].final,   genLeptons[i].final.key(), *electrons, electronSelection_, electronPtScale_, electronDetaScale_, electronDphiScale_, electronCandMatches);
        }
    }
    if (debug_) std::cout << std::endl;

    if (debug_) std::cout << "Arbitrating electron matches: " << std::endl;
    resolve(electronCandMatches);
    if (debug_) std::cout << "Arbitrating muon matches: " << std::endl;
    resolve(muonCandMatches);

    if (debug_) std::cout << "Removing FSR photons for objects matched to status == 3 leptons" << std::endl;
    for (int i = 0, n = electronCandMatches.size(); i < n; ++i) {
        const MatchStruct & match = electronCandMatches[i]; 
        if (match.isDropped() || match.genStatus != 3) continue;
        for (int j = 0; j < 4; ++j) {
            if (match.gen == int(genLeptons[j].final.key())) {
                genLeptons[j].photons.clear();
                break;
            }
        }
    }
    for (int i = 0, n = muonCandMatches.size(); i < n; ++i) {
        const MatchStruct & match = muonCandMatches[i]; 
        if (match.isDropped() || match.genStatus != 3) continue;
        for (int j = 0; j < 4; ++j) {
            if (match.gen == int(genLeptons[j].final.key())) {
                genLeptons[j].photons.clear();
                break;
            }
        }
    }

    for (int i = 0; i < 4; ++i) {
        for (int k = 0, n = genLeptons[i].photons.size(); k < n; ++k) {
            if (debug_) std::cout << "Matching gen photon " << k << " from FSR of lepton " << i << " to photons: " << std::endl;
            reco::GenParticleRef phot = genLeptons[i].photons[k];
            match(*phot, phot.key(), *photons, photonSelection_, photonPtScale_, photonDetaScale_, photonDphiScale_, photonCandMatches);
        }
    }

    if (debug_) std::cout << "Arbitrating photon matches: " << std::endl;
    resolve(photonCandMatches);
 
    if (debug_) std::cout << std::endl;
    std::map<int,reco::CandidatePtr> genToReco; 
    // now we fill the reco-gen matches
    for (int i = 0, n = electronCandMatches.size(); i < n; ++i) {
        const MatchStruct & match = electronCandMatches[i]; 
        if (match.isDropped()) continue;
        if (debug_) printf("Final match for electron: gen %03d -> rec %03d\n", match.gen, match.rec);
        electronMatches[match.rec] = match.gen; genToReco[match.gen] = electrons->ptrAt(match.rec);
    }
    for (int i = 0, n = photonCandMatches.size(); i < n; ++i) {
        const MatchStruct & match = photonCandMatches[i]; 
        if (match.isDropped()) continue;
        if (debug_) printf("Final match for photon:   gen %03d -> rec %03d\n", match.gen,match.rec);
        photonMatches[match.rec] = match.gen; genToReco[match.gen] = photons->ptrAt(match.rec);
    }
    for (int i = 0, n = muonCandMatches.size(); i < n; ++i) {
        const MatchStruct & match = muonCandMatches[i]; 
        if (match.isDropped()) continue;
        if (debug_) printf("Final match for muon:     gen %03d -> rec %03d\n", match.gen,match.rec);
        muonMatches[match.rec] = match.gen; genToReco[match.gen] = muons->ptrAt(match.rec);
    }

    typedef edm::Association<reco::GenParticleCollection> GPA;
    std::auto_ptr<GPA> outmap(new GPA(genParticles));
    GPA::Filler filler(*outmap);
    if (debug_) std::cout << "Will fill for electrons: " << electrons.id() << std::endl; 
    bool fillMu = true, fillPho = true;
    if (muons.id() == electrons.id()) {
        for (int i = 0, n = muonMatches.size(); i < n; ++i) {
            if (muonMatches[i] != -1) electronMatches[i] = muonMatches[i];
        }
        fillMu = false;
    } else if (debug_) { std::cout << "Will fill for muons: " << muons.id() << std::endl; }
    if (photons.id() == electrons.id()) {
        for (int i = 0, n = photonMatches.size(); i < n; ++i) {
            if (photonMatches[i] != -1) electronMatches[i] = photonMatches[i];
        }
        fillPho = false;
    } else if (debug_) { std::cout << "Will fill for photons: " << photons.id() << std::endl; }
    filler.insert(electrons, electronMatches.begin(), electronMatches.end());
    if (fillMu) filler.insert(muons,     muonMatches.begin(), muonMatches.end());
    if (fillPho) filler.insert(photons,   photonMatches.begin(), photonMatches.end());
    filler.fill();
    iEvent.put(outmap);

    // now we try to fill the correct reco final state
    std::auto_ptr<std::vector<reco::CompositeCandidate> > out(new std::vector<reco::CompositeCandidate>(1));
    reco::CompositeCandidate &zz = out->front();
    reco::CompositeCandidate *z1 = new reco::CompositeCandidate(zs[0]->charge(), zs[0]->p4(), zs[0]->vertex(), zs[0]->pdgId());
    reco::CompositeCandidate *z2 = new reco::CompositeCandidate(zs[1]->charge(), zs[1]->p4(), zs[1]->vertex(), zs[1]->pdgId());
    fillZBoson(genLeptons, 0, genToReco, *z1);
    fillZBoson(genLeptons, 2, genToReco, *z2);
    zz.setP4(z1->p4() + z2->p4());
    zz.setStatus(std::min<int>(z1->status(), z2->status()));
    std::auto_ptr<reco::Candidate> pz1(z1);
    std::auto_ptr<reco::Candidate> pz2(z2);
    zz.addDaughter(pz1);
    zz.addDaughter(pz2);
    iEvent.put(out);
    
    return true;
}

bool
RecoFinalStateMatcher::fillGenLepton(GenLeptonWithFsr &lepton) {
    lepton.final = lepton.initial;
    while(lepton.final->status() != 1) {
        reco::GenParticleRef next;
        for (int j = 0, n = lepton.final->numberOfDaughters(); j < n; ++j) {
            reco::GenParticleRef dau = lepton.final->daughterRef(j);
            if (dau->pdgId() == lepton.initial->pdgId()) {
                next = dau;
            } else if (dau->pdgId() == 22 && dau->status() == 1) {
                if( dau->pt() > genPhotonPtMin_) lepton.photons.push_back(dau);
            } else {
                std::cerr << "Strange daughter " << dau->pdgId() << ", " << dau->status() << 
                    " of " << lepton.final->pdgId() << ", " << lepton.final->status() << std::endl;
                return false;
            }
        }
        if (next.isNull()) {
            std::cerr << "Missing daughter of lepton " << lepton.final->pdgId() << ", " << lepton.final->status() << std::endl;
            return false;
        }
        lepton.final = next;
    }
    return true; 
}

void
RecoFinalStateMatcher::match(const reco::GenParticle &gen, int i,
                   const reco::CandidateView &reco, const StringCutObjectSelector<reco::Candidate> &isgood,
                   double ptScale, double detaScale, double dphiScale,
                   std::vector<MatchStruct> & out) 
{
    if (debug_) printf("    Searching for matches for candidate %03d pdgId %+2d, status %1d with pt %7.3f, eta %+6.4f, phi %+6.4f\n", i, gen.pdgId(), gen.status(), gen.pt(), gen.eta(), gen.phi());  
    bool nmatch = 0;
    for (int j = 0, n = reco.size(); j < n; ++j) {
        const reco::Candidate &r = reco[j];
        if (!isgood(r)) continue;
        double deta = fabs(gen.eta() - r.eta()), dphi = fabs(::deltaPhi(gen.phi(), r.phi())), dpt = fabs(r.pt() - gen.pt())/gen.pt();
        double detan = deta/detaScale, dphin = dphi/dphiScale, dptn = dpt/ptScale;
        double dist = sqrt(detan*detan + dphin*dphin + dptn*dptn);
        if (detan < 1 && dphin < 1 && dptn < 1) {
            if (debug_) printf("        plausible reco object %03d pdgId %+2d, pt %7.3f, eta %+6.4f, phi %+6.4f: deta %7.5f (norm %7.5f),  dphi %7.5f (norm %7.5f), dPtRel %8.4f (norm %.4f), dist %.5f\n",
                                        j, r.pdgId(), r.pt(), r.eta(), r.phi(), deta, detan, dphi, dphin, dpt, dptn, dist);
            out.push_back(MatchStruct(dist, i, gen.status(), j));
            nmatch++;
        }
    }
    if (debug_ && nmatch == 0) {
        for (int j = 0, n = reco.size(); j < n; ++j) {
            const reco::Candidate &r = reco[j];
            if (!isgood(r)) continue;
            double deta = fabs(gen.eta() - r.eta()), dphi = fabs(::deltaPhi(gen.phi(), r.phi())), dpt = fabs(r.pt() - gen.pt())/gen.pt();
            double detan = deta/detaScale, dphin = dphi/dphiScale, dptn = dpt/ptScale;
            double dist = sqrt(detan*detan + dphin*dphin + dptn*dptn);
            printf("        too far   reco object %03d pdgId %+2d, pt %7.3f, eta %+6.4f, phi %+6.4f: deta %7.5f (norm %9.4f),  dphi %7.5f (norm %9.4f), dPtRel %8.4f (norm %9.4f), dist %9.4f\n",
                    j, r.pdgId(), r.pt(), r.eta(), r.phi(), deta, detan, dphi, dphin, dpt, dptn, dist);
        }
    }
}

void
RecoFinalStateMatcher::resolve(std::vector<MatchStruct> & out) 
{
    std::sort(out.begin(), out.end());
    for (int i = 0, n = out.size(); i < n; ++i) {
        if (out[i].isDropped()) continue;
        if (debug_) printf("    accepting match (%03d status %1d, %03d) with norm dist %.5f\n", out[i].gen, out[i].genStatus, out[i].rec, out[i].dist);
        for (int j = i+1; j < n; ++j) {
            if (out[j].isDropped()) continue;
            if (out[i].overlap(out[j]))  {
                if (debug_) printf("        disqualifying match (%03d status %1d, %03d) with norm dist %.5f\n", out[j].gen, out[j].genStatus, out[j].rec, out[j].dist);
                out[j].drop();
            }
        }
    }
}

void
RecoFinalStateMatcher::fillZBoson(const std::vector<GenLeptonWithFsr> &leptons, int il1, std::map<int,reco::CandidatePtr> &genToReco, reco::CompositeCandidate &cc) 
{
    int status = 2; // fully reconstructed
    for (int i = il1; i <= il1+1; ++i) {
        reco::GenParticleRef   gen = leptons[i].final;
        reco::CandidatePtr rec = genToReco[gen.key()];
        reco::ShallowClonePtrCandidate scc(rec, gen->charge(), gen->p4(), gen->vertex());
        scc.setPdgId(gen->pdgId());
        scc.setStatus(rec.isNonnull() ? 1 : 0);
        if (rec.isNull()) status = 0;
        cc.addDaughter(scc);
    }
    for (int i = il1; i <= il1+1; ++i) {
        reco::GenParticleRefVector phos = leptons[i].photons;
        for (int j = 0, n = phos.size(); j < n; ++j) {
            reco::GenParticleRef   gen = phos[j];
            reco::CandidatePtr rec = genToReco[gen.key()];
            reco::ShallowClonePtrCandidate scc(rec, gen->charge(), gen->p4(), gen->vertex());
            scc.setPdgId(gen->pdgId());
            scc.setStatus(rec.isNonnull() ? 1 : 0);
            if (rec.isNull() && status > 0) status = 1;
            cc.addDaughter(scc);
        }
    } 
    cc.setStatus(status);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoFinalStateMatcher);
