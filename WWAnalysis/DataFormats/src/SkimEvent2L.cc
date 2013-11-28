#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include <algorithm>
#include <cmath>

namespace utils2L {
    typedef StringObjectFunction<reco::Candidate,true>    Function;
    std::map<std::string,Function> functionCache2L_;
}

reco::SkimEvent2L::~SkimEvent2L() {}

reco::SkimEvent2L::SkimEvent2L(const reco::CompositeCandidate &src) :
    pat::CompositeCandidate(reco::CompositeCandidate(src.charge(), src.p4(), src.vertex(), src.pdgId(), src.status()))
{
    const reco::Candidate *srcl1 = src.daughter(0);
    const reco::Candidate *srcl2 = src.daughter(1);
    if (srcl1->pt() > srcl2->pt()) std::swap(srcl1, srcl2);
    addDaughter(*srcl1);
    addDaughter(*srcl2);
    for (unsigned int i = 2, n = src.numberOfDaughters(); i < n; ++i) {
        const reco::Candidate &pho = *src.daughter(i);
        if (pho.pdgId() != 22) throw cms::Exception("CorruptData") << "Daughter " << i << " of Z is not a photon!\n";
        addDaughter(pho);
    }
}

reco::SkimEvent2L::SkimEvent2L(const pat::CompositeCandidate &src) :
    pat::CompositeCandidate(src)
{
    clearDaughters();
    const reco::Candidate *srcl1 = src.daughter(0);
    const reco::Candidate *srcl2 = src.daughter(1);
    if (srcl1->pt() > srcl2->pt()) std::swap(srcl1, srcl2);
    addDaughter(*srcl1);
    addDaughter(*srcl2);
    for (unsigned int i = 2, n = src.numberOfDaughters(); i < n; ++i) {
        const reco::Candidate &pho = *src.daughter(i);
        if (pho.pdgId() != 22) throw cms::Exception("CorruptData") << "Daughter " << i << " of Z is not a photon!\n";
        addDaughter(pho);
    }
}

void reco::SkimEvent2L::setVertices(const edm::Handle<reco::VertexCollection> &vtxColl) 
{
    vertices_ = *vtxColl;
}

void reco::SkimEvent2L::setPFMet(const edm::Handle<reco::PFMETCollection> &pfMetColl)
{
    pfMet_ = pfMetColl->front();
}

void reco::SkimEvent2L::setPileupInfo(const edm::Handle<std::vector<PileupSummaryInfo> >& pu)
{
    for (size_t i = 0; i < pu->size(); i++) {
      if ((*pu)[i].getBunchCrossing() == 0) {
	numsimvertices_ = (*pu)[i].getPU_NumInteractions();
	numTrueInteractions_ = (*pu)[i].getTrueNumInteractions();
      }
    }
}

float reco::SkimEvent2L::lval(unsigned int il, const std::string &muExpr, const std::string &eleExpr) const 
{
    std::map<std::string,utils2L::Function>::const_iterator itMu = utils2L::functionCache2L_.find(muExpr);
    if (itMu == utils2L::functionCache2L_.end()) {
       itMu = utils2L::functionCache2L_.insert(std::make_pair(muExpr, utils2L::Function(muExpr))).first;
    }
    std::map<std::string,utils2L::Function>::const_iterator itEl = utils2L::functionCache2L_.find(eleExpr);
    if (itEl == utils2L::functionCache2L_.end()) {
       itEl = utils2L::functionCache2L_.insert(std::make_pair(eleExpr, utils2L::Function(eleExpr))).first;
    }
    const reco::Candidate &c = l(il);
    if      (typeid(c) == typeid(pat::Muon))     return itMu->second(c);
    else if (typeid(c) == typeid(pat::Electron)) return itEl->second(c);
    else throw cms::Exception("WrongType") << "Lepton " << il << " is of type " << typeid(c).name() << "\n";
}

float reco::SkimEvent2L::lval(unsigned int il, const std::string &expr) const
{
     std::map<std::string,utils2L::Function>::const_iterator itExpr = utils2L::functionCache2L_.find(expr);
     if (itExpr == utils2L::functionCache2L_.end()) {
        itExpr = utils2L::functionCache2L_.insert(std::make_pair(expr, utils2L::Function(expr))).first;
     }
     return itExpr->second(l(il));
}


reco::Particle::LorentzVector reco::SkimEvent2L::genllp4() const 
{
    reco::GenParticleRef gl1 = genl(0), gl2 = genl(1);
    if (gl1.isNonnull() && gl2.isNonnull()) {
        return gl1->p4() + gl2->p4();
    } else {
        return reco::Particle::LorentzVector();
    }
}

reco::GenParticleRef reco::SkimEvent2L::genz() const 
{
    reco::GenParticleRef gl1 = genl(0), gl2 = genl(1);
    // search for mothers = not leptons
    while (gl1.isNonnull() && gl1->numberOfMothers() != 0 && 11 <= abs(gl1->pdgId()) && abs(gl1->pdgId()) <= 15) {
        gl1 = gl1->motherRef(0);
    }
    while (gl2.isNonnull() && gl2->numberOfMothers() != 0 && 11 <= abs(gl2->pdgId()) && abs(gl2->pdgId()) <= 15) {
        gl2 = gl2->motherRef(0);
    }
    if (gl1.isNonnull() && gl2.isNonnull() && gl1 == gl2) {
        return gl1;
    } else {
        return reco::GenParticleRef();
    }
}


void reco::SkimEvent2L::setGenMatches(const edm::Association<reco::GenParticleCollection> &genMatch) 
{
    int n = numberOfDaughters();
    matches_.clear();
    for (int i = 0; i < n; ++i) {
        const reco::Candidate &proxy = *daughter(i);
        reco::GenParticleRef match;
        if (proxy.hasMasterClonePtr())   match = genMatch[proxy.masterClonePtr()];
        else if (proxy.hasMasterClone()) match = genMatch[proxy.masterClone()];
        matches_.push_back(match); 
    }
}
