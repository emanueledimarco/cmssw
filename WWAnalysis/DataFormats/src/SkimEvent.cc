#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "Math/VectorUtil.h"

#include <TVector3.h>

#include <iostream>
#include <algorithm>


const int reco::SkimEvent::channel() const { 
    switch(hypo_) {
        case WWMUMU: return 0;
        case WWELEL: return 1;
        case WWELMU: return 2;
        case WWMUEL: return 3;
    }
	return -1;
}

std::vector<std::string> reco::SkimEvent::jecFiles_;

const bool reco::SkimEvent::peaking() const {

    if (getMotherID(0).isNonnull() && getMotherID(1).isNonnull() && getMotherID(0) == getMotherID(1)) return true;

    return false;
}

const reco::GenParticle *reco::SkimEvent::genParticle(size_t i) const {
    return  isMuon(i) ? getMuon(i)->genLepton() : getElectron(i)->genLepton();
}

const reco::GenParticleRef reco::SkimEvent::getMotherID(size_t i) const {

    const reco::GenParticle *match = isMuon(i) ? getMuon(i)->genLepton() : getElectron(i)->genLepton();

    if( !match )                    return reco::GenParticleRef();
    if( !match->numberOfMothers() ) return reco::GenParticleRef();

    int pdgId = match->pdgId();

    reco::GenParticleRef mother = match->motherRefVector()[0];
    if( !mother.isNonnull() )        return reco::GenParticleRef();
    if( !mother->numberOfMothers() ) return reco::GenParticleRef();
    while(mother->pdgId() == pdgId) {
        mother = mother->motherRefVector()[0];
        if( !mother.isNonnull() )        return reco::GenParticleRef();
        if( !mother->numberOfMothers() ) return reco::GenParticleRef();
    }
    return mother;
}

const bool reco::SkimEvent::isHardMuID(size_t i) const {

    if( i >= leps_.size() ) return false;

    if( isElectron(i) )  return true;
    if( isMuon(i) )  {
        pat::Muon const * const mu = getMuon(i);
        return (mu->isGlobalMuon() && mu->isTrackerMuon());
    }
    return false;

}

const bool reco::SkimEvent::passesSmurfMuonID() const {

    switch(hypo_) {
        case WWELMU: return  isHardMuID(indexByPt(1)); break;
        case WWMUEL: return  isHardMuID(indexByPt(0)); break;
        case WWMUMU: return (isHardMuID(indexByPt(0)) && isHardMuID(indexByPt(1))); break;
        case WWELEL: return true; break;
        default    : return true; break;
    }
    return true;

}

struct indexValueStruct {
    indexValueStruct(const float &v, size_t i) : value(v), index(i) {}
    float value;
    size_t index;
};

bool lowToHigh(const indexValueStruct &a, const indexValueStruct &b) {
    return a.value < b.value;
}

bool highToLow(const indexValueStruct &a, const indexValueStruct &b) {
    return a.value > b.value;
}

// To be kept in synch with the enumerator definitions in SkimEvent.h file
const std::string reco::SkimEvent::hypoTypeNames[] = { "undefined", "WELNU", "WMUNU", "WWELEL", "WWELMU", "WWMUEL", "WWMUMU"};


std::string reco::SkimEvent::hypoTypeName(reco::SkimEvent::hypoType a){
  if(int(a) < int(hypoTypeSize) && int(a)>0) return hypoTypeNames[int(a)];
  return "undefined";
}



reco::SkimEvent::hypoType reco::SkimEvent::hypoTypeByName(const std::string &name){
  hypoType size = hypoTypeSize;
  int index = std::find(hypoTypeNames, hypoTypeNames+size, name)-hypoTypeNames;
  if(index == size) return undefined; // better this or throw() ?

  // cast
  return hypoType(index);
}



reco::SkimEvent::SkimEvent() : 
        hypo_(-1), sumPts_(0)/*, jec_(0), vtxPoint_(0,0,0) */{ }

reco::SkimEvent::SkimEvent(const reco::SkimEvent::hypoType &h) :
        hypo_(h), sumPts_(0)/*, jec_(0), vtxPoint_(0,0,0) */{ }



// set GenParticles
void reco::SkimEvent::setGenParticles(const edm::Handle<reco::GenParticleCollection> & h) {
   genParticles_.clear();
   for(size_t i=0; i<h->size(); ++i) {
    genParticles_.push_back( reco::GenParticleRef(h,i) );
   }
}


//EDM RefToBase implementation
void reco::SkimEvent::setLepton(const edm::Handle<edm::View<reco::RecoCandidate> > &h,size_t i){
  //std::cout << "setting lepton with collection ID: " << h->ptrAt(i).id() << std::endl;
  leps_.push_back( h->ptrAt(i) );
}

void reco::SkimEvent::setExtraLepton(const edm::Handle<edm::View<reco::RecoCandidate> > &h,size_t i){
  extraLeps_.push_back( h->ptrAt(i) );
//extraLeps_.push_back(refToCand(h,i));
}

void reco::SkimEvent::setSoftMuon(const edm::Handle<edm::View<reco::RecoCandidate> > &h,size_t i){
  softMuons_.push_back(h->ptrAt(i) );
//softMuons_.push_back(refToCand(h,i));
}

// Old implementation
// void reco::SkimEvent::setLepton(const pat::Electron& ele){
//   leps_.push_back(ele);
// }
// 
// void reco::SkimEvent::setLepton(const pat::Muon& mu){
//   leps_.push_back(mu);
// }
// 
// void reco::SkimEvent::setExtraLepton(const pat::Electron& ele){
//   extraLeps_.push_back(ele);
// }
// 
// void reco::SkimEvent::setExtraLepton(const pat::Muon& mu){
//   extraLeps_.push_back(mu);
// }
// 
// void reco::SkimEvent::setSoftMuon(const pat::Muon& mu){
//   softMuons_.push_back(mu);
// }


void reco::SkimEvent::setJets(const edm::Handle<pat::JetCollection> & jH) {

    jets_.clear();

    for(size_t i=0;i<jH->size();++i) 
      jets_.push_back(pat::JetRef(jH,i));

    //sortJetsByPt();

}

void reco::SkimEvent::setTagJets(const edm::Handle<pat::JetCollection> & jH) {

    tagJets_.clear();

    for(size_t i=0;i<jH->size();++i) 
      tagJets_.push_back(pat::JetRef(jH,i));

    //sortTagJetsByPt();

}

void reco::SkimEvent::setTCMet(const edm::Handle<reco::METCollection> & mH) {
    tcMet_ = reco::METRef(mH,0);
}

void reco::SkimEvent::setPFMet(const edm::Handle<reco::PFMETCollection> & mH) {
    pfMet_ = reco::PFMETRef(mH,0);
}

void reco::SkimEvent::setChargedMet(const reco::PFMET & chMET) {
    chargedMet_ = chMET;
}

void reco::SkimEvent::setVtxSumPts(const edm::Handle<edm::ValueMap<float> > &s) {

    for(size_t i=0;i<vtxs_.size();++i) sumPts_.push_back( (*s)[vtxs_[i]] );

}

void reco::SkimEvent::setVtxSumPt2s(const edm::Handle<edm::ValueMap<float> > &s) {

    for(size_t i=0;i<vtxs_.size();++i) sumPt2s_.push_back( (*s)[vtxs_[i]] );

}

void reco::SkimEvent::setVertex(const edm::Handle<reco::VertexCollection> & vtxH) {

    for(size_t i=0;i<vtxH->size();++i) vtxs_.push_back(reco::VertexRef(vtxH,i));
}


//Lepton variables

const int reco::SkimEvent::nLep(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = leps_.size(); 
    else for(size_t i=0;i<leps_.size();++i) if(leps_[i]->pt() > minPt) count++;
    return count;
} 

const int reco::SkimEvent::nExtraLep(float minPt) const { 
    int count = 0;
     if(minPt < 0) count = extraLeps_.size(); 
     else for(size_t i=0;i<extraLeps_.size();++i) if(extraLeps_[i]->pt() > minPt) count++;
//    for(size_t i=0;i<extraLeps_.size();++i) if(extraLeps_[i]->pt() > minPt && passesIP(extraLeps_[i]) ) count++;
    return count;
}


// Doesn't veto soft muons around jets by default
// i.e. vetoJets == -1
// if vetoJets != -1 then vetoJets is taken as the jet pt threshold
const int reco::SkimEvent::nSoftMu(float minPt, float vetoJets, float dRCut) const { 
    int count = 0;
    if(minPt < 0 && vetoJets < 0) count = softMuons_.size(); 
    else if(vetoJets < 0) {
        for(size_t i=0;i<softMuons_.size();++i) if(softMuons_[i]->pt() > minPt) count++;
    } else {
        for(size_t i=0;i<softMuons_.size();++i) {
            bool toCount = true;
            for(size_t j=0; j<jets_.size();++j){
                if( jetPt(j,true) <= vetoJets) continue;
                double dR = fabs(ROOT::Math::VectorUtil::DeltaR(softMuons_[i]->p4(),jets_[j]->p4()) );
                if(dR < dRCut) toCount = false;
            }
            if(toCount && softMuons_[i]->pt() > minPt) count++;
        }
    }
    return count;
} 


const int reco::SkimEvent::pdgId(size_t i) const { 
    if(i < leps_.size()) return leps_[i]->pdgId();
    else return -9999;
}


const bool reco::SkimEvent::isMuon(size_t  i) const {
  return isMuon(leps_[i]);
}

const bool reco::SkimEvent::isElectron(size_t  i) const {
  return isElectron(leps_[i]);
}

const bool reco::SkimEvent::isMuon(const refToCand &c) const {
  return (abs(c->pdgId()) == 13);
}

const bool reco::SkimEvent::isElectron(const refToCand &c) const {
  return (abs(c->pdgId()) == 11);
}

const pat::Muon *  reco::SkimEvent::getMuon(size_t  i) const {
  return getMuon(leps_[i]);
}

const pat::Electron  *  reco::SkimEvent::getElectron(size_t  i) const {
  return getElectron(leps_[i]);
}

const pat::Muon * reco::SkimEvent::getMuon(const refToCand &c) const {
  return static_cast<const pat::Muon*>(c.get());
    
}

const pat::Electron * reco::SkimEvent::getElectron(const refToCand &c) const {
  return static_cast<const pat::Electron*>(c.get());
}

const float reco::SkimEvent::pt(size_t i) const {
  if(i >= leps_.size()) return -9999.0;
  return leps_[i]->pt();
}


const int reco::SkimEvent::passCustom(size_t i, const std::string &muStr, const std::string &elStr) const {
  if(i >= leps_.size())    return 0;
  else if( isElectron(i) ) return StringCutObjectSelector<pat::Electron>(elStr,true)( *getElectron(i) );
  else if( isMuon(i)     ) return StringCutObjectSelector<pat::Muon    >(muStr,true)( *getMuon(i)     );
  else                     return 0;
}

const float reco::SkimEvent::leptBdt(size_t i) const {
  if(i >= leps_.size()) return -9999.0;
  if( isElectron(i) ) return getElectron(i)->userFloat("bdt");
  else                return 999999.;
}

const float reco::SkimEvent::leptLH(size_t i) const {
  if(i >= leps_.size()) return -9999.0;
  if( isElectron(i) ) return getElectron(i)->electronID("egammaIDLikelihood");
  else                return 999999.;
}

const float reco::SkimEvent::eta(size_t i) const {
  if(i < leps_.size()) return leps_[i]->eta();
  else return -9999.0;
}

const float reco::SkimEvent::etaSC(size_t i) const {
  if(i >- leps_.size()) return -9999.0;
  if( isElectron(i) ) return getElectron(i)->superCluster()->eta();
  else                return leps_[i]->eta();
}

const float reco::SkimEvent::nBrem(size_t i) const {
  if(i >- leps_.size()) return -9999.0;
  if( isElectron(i) ) return getElectron(i)->numberOfBrems();
  else                return 0;
}

const float reco::SkimEvent::phi(size_t i) const {
  if(i < leps_.size()) return leps_[i]->phi();
  else return -9999.0;
}


const int reco::SkimEvent::q(size_t i) const{
  if(i < leps_.size()) return leps_[i]->charge();
  else return -9999.0;
}


//Jet variables
void reco::SkimEvent::setupJEC(const std::string &l2File, const std::string &l3File, const std::string &residualFile) {

    jecFiles_.clear();
    jecFiles_.push_back(l2File);
    jecFiles_.push_back(l3File);
    if(residualFile.compare("")) jecFiles_.push_back(residualFile);

}


// void reco::SkimEvent::setupJEC(const JetCorrector *c) {
//     jec_ = 

const int reco::SkimEvent::nJets(float minPt, int applyCorrection,int applyID) const {
  return nCentralJets(minPt,99.9, applyCorrection,applyID);
}

const bool reco::SkimEvent::isThisJetALepton(pat::JetRef jet, float drCut) const {
    bool thisJetIsLepton(false);
    for(size_t j=0; j<leps_.size();++j){
        double dR = fabs(ROOT::Math::VectorUtil::DeltaR(jet->p4(),leps_[j]->p4()) );
        if(dR < drCut){ 
            thisJetIsLepton = true;
            break;
        }
    }

    return thisJetIsLepton;
}

const bool reco::SkimEvent::passJetID(pat::JetRef jet, int applyID) const{
  // no ID  
  if(applyID == 0) return true;
  
  // old ID
  else if(applyID == 1) {  
    unsigned int multiplicity = jet->neutralMultiplicity () + jet->chargedMultiplicity ();
    if(jet->neutralEmEnergyFraction() >=0.99 || 
            jet->neutralHadronEnergyFraction() >=0.99 ||
            multiplicity==0 ) return false;
    if(fabs(jet->eta())<2.4){
        if(jet->chargedEmEnergyFraction() >=0.99 ||
            jet->chargedHadronEnergyFraction() == 0 ||
            jet->chargedMultiplicity()==0 ) return false;
    }
  }
  
  // MVA ID loose
  else if(applyID == 4) {
    if(jet->userInt("jetId") >= 4) return true;
    else return false;
  }
  
  // MVA ID medium
  else if(applyID == 5) {
    if(jet->userInt("jetId") >= 6) return true;
    else return false;
  }  
  
  // MVA ID tight
  else if(applyID == 6) {
    if(jet->userInt("jetId") >= 7) return true;
    else return false;
  }  
  
  return false;
}

const float reco::SkimEvent::dPhiJetllInDegrees(size_t leadingIndex,float minPt,float eta,int applyCorrection,int applyID) const {
    return dPhiJetll(leadingIndex,minPt,eta,applyCorrection,applyID)/M_PI * 180.;
}

const float reco::SkimEvent::dPhiJetll(size_t leadingIndex,float minPt,float eta,int applyCorrection,int applyID) const {

    size_t count = 0, newIndex = 0;
    for(size_t i=0;i<jets_.size();++i) {
        if(!(passJetID(jets_[i],applyID)) ) continue;
        if( std::fabs(jets_[i]->eta()) >= eta) continue;
        if( jetPt(i,applyCorrection) <= minPt) continue;

        if(isThisJetALepton(jets_[i]))  continue;
        // When count becomes higher than leadingIndex, we've found the leadingIndexth jet that passes all the cuts
        if(++count > leadingIndex) {
            newIndex = i;
            break;
        }
    }
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4() + leps_[1]->p4(), jets_[newIndex]->p4()) );   

}

const pat::Jet* reco::SkimEvent::leadingJet(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(++count > index) return jets_[i].get();
    }
    return new pat::Jet();
}

const float reco::SkimEvent::leadingJetBtag(size_t index, std::string discriminator, float minPt,float eta,int applyCorrection,int applyID, float dzCut) const {

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(jets_[i]->hasUserFloat("dz") && fabs(jets_[i]->userFloat("dz")) > dzCut) continue;
      if(++count > index) return jets_[i]->bDiscriminator(discriminator);
    }
    return -9999.9;

}

const float reco::SkimEvent::leadingJetPhi(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(++count > index) return jets_[i]->phi();
    }
    return -9999.9;

}

const float reco::SkimEvent::leadingJetEta(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(++count > index) return jets_[i]->eta();
    }
    return -9999.9;

}

const int reco::SkimEvent::leadingJetId(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(++count > index) return jets_[i]->userInt("jetId");
    }
    return -9999.9;
}

const float reco::SkimEvent::leadingJetMva(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(++count > index) return jets_[i]->userFloat("jetMva");
    }
    return -9999.9;
}
const float reco::SkimEvent::leadingJetPt(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(++count > index) return jetPt(i,applyCorrection);
    }
    return -9999.9;
}


const float reco::SkimEvent::leadingVBFJetPhi(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {
    if(jets_.size() < 2) return -9999.;
    float jetEta1 = leadingJetEta(0,minPt,eta,applyCorrection,applyID);
    float jetEta2 = leadingJetEta(1,minPt,eta,applyCorrection,applyID);
    float jetEtaMax = std::max(jetEta1,jetEta2);
    float jetEtaMin = std::min(jetEta1,jetEta2);

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if( jets_[i]->eta() >= jetEtaMax ) continue;
      if( jets_[i]->eta() <= jetEtaMin ) continue;
      if(++count > index) return jets_[i]->phi();
    }
    return -9999.9;
}

const float reco::SkimEvent::leadingVBFJetEta(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {
    if(jets_.size() < 2) return -9999.;
    float jetEta1 = leadingJetEta(0,minPt,eta,applyCorrection,applyID);
    float jetEta2 = leadingJetEta(1,minPt,eta,applyCorrection,applyID);
    float jetEtaMax = std::max(jetEta1,jetEta2);
    float jetEtaMin = std::min(jetEta1,jetEta2);

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if( jets_[i]->eta() >= jetEtaMax ) continue;
      if( jets_[i]->eta() <= jetEtaMin ) continue;
      if(++count > index) return jets_[i]->eta();
    }
    return -9999.9;

}

const int reco::SkimEvent::leadingVBFJetId(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {
    if(jets_.size() < 2) return -9999.;
    float jetEta1 = leadingJetEta(0,minPt,eta,applyCorrection,applyID);
    float jetEta2 = leadingJetEta(1,minPt,eta,applyCorrection,applyID);
    float jetEtaMax = std::max(jetEta1,jetEta2);
    float jetEtaMin = std::min(jetEta1,jetEta2);

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if( jets_[i]->eta() >= jetEtaMax ) continue;
      if( jets_[i]->eta() <= jetEtaMin ) continue;
      if(++count > index) return jets_[i]->userInt("jetId");
    }
    return -9999.9;
}

const float reco::SkimEvent::leadingVBFJetMva(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {
    if(jets_.size() < 2) return -9999.;
    float jetEta1 = leadingJetEta(0,minPt,eta,applyCorrection,applyID);
    float jetEta2 = leadingJetEta(1,minPt,eta,applyCorrection,applyID);
    float jetEtaMax = std::max(jetEta1,jetEta2);
    float jetEtaMin = std::min(jetEta1,jetEta2);

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if( jets_[i]->eta() >= jetEtaMax ) continue;
      if( jets_[i]->eta() <= jetEtaMin ) continue;
      if(++count > index) return jets_[i]->userInt("jetMva");
    }
    return -9999.9;
}

const float reco::SkimEvent::leadingVBFJetPt(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {
    if(jets_.size() < 2) return -9999.;
    float jetEta1 = leadingJetEta(0,minPt,eta,applyCorrection,applyID);
    float jetEta2 = leadingJetEta(1,minPt,eta,applyCorrection,applyID);
    float jetEtaMax = std::max(jetEta1,jetEta2);
    float jetEtaMin = std::min(jetEta1,jetEta2);

    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if( jets_[i]->eta() >= jetEtaMax ) continue;
      if( jets_[i]->eta() <= jetEtaMin ) continue;
      if(++count > index) return jetPt(i,applyCorrection);
    }
    return -9999.9;
}

const int reco::SkimEvent::nCentralJets(float minPt,float eta,int applyCorrection,int applyID) const {

    int count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      count++;
    }
    return count;
}

const bool reco::SkimEvent::passesDPhillJet(float ptMin, float eta,int applyCorrection,int applyID) const {
    float dphi = 0, ptMax = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      float pt = jetPt(i,applyCorrection);
      if (pt > ptMax) {
        ptMax = pt;
        dphi  = fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4()+leps_[1]->p4(), jets_[i]->p4()) );
      }
    }
    return (ptMax <= ptMin || dphi / M_PI * 180. < 165.0);
}

const float reco::SkimEvent::dPhillLeadingJet(float eta,int applyCorrection,int applyID) const {
    float dphi = 0, ptMax = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      float pt = jetPt(i,applyCorrection);
      if (pt > ptMax) {
        ptMax = pt;
        dphi  = fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4()+leps_[1]->p4(), jets_[i]->p4()) );
      }
    }
    return dphi;
}

const int reco::SkimEvent::leadingJetIndex(size_t index, float minPt,float eta,int applyCorrection,int applyID) const {
  
    size_t count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if(++count > index) return i;
    }
    return -1;
}

const float reco::SkimEvent::dPhilljetjet(float eta,int applyCorrection,int applyID) const {
    float dphi = -1;
    int jet1 = leadingJetIndex(0,0,eta,applyCorrection,applyID);
    int jet2 = leadingJetIndex(1,0,eta,applyCorrection,applyID);

    if (jet1 != -1 && jet2 != -1) dphi  = fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4()+leps_[1]->p4(), jets_[jet1]->p4()+jets_[jet2]->p4()) );

    return dphi;
}



const float reco::SkimEvent::jetPt(size_t i, int applyCorrection) const {
    return jetPt(jets_[i].get(),applyCorrection);
}

const float reco::SkimEvent::jetPt(const pat::Jet *j, int applyCorrection) const {
//   if(applyCorrection) return jets_[i]->correctedJet("L3Absolute","none").pt();
  if(applyCorrection) return j->pt();
  else                return j->correctedJet("Uncorrected","none").pt();
}

const float reco::SkimEvent::tagJetPt(size_t i, int applyCorrection) const {
//   if(applyCorrection) return tagJets_[i]->correctedJet("L3Absolute","none").pt();
  if(applyCorrection) return tagJets_[i]->pt();
  else                return tagJets_[i]->correctedJet("Uncorrected","none").pt();
}



//Event variables
const float reco::SkimEvent::pfMet() const {

    if(pfMet_.isNonnull()) return pfMet_->pt();
    else return -9999.0;
}

const float reco::SkimEvent::pfMetPhi() const {

    if(pfMet_.isNonnull()) return pfMet_->phi();
    else return -9999.0;
}

const float reco::SkimEvent::tcMet() const {

    if(tcMet_.isNonnull()) return tcMet_->pt();
    else return -9999.0;
}

const float reco::SkimEvent::tcMetPhi() const {

    if(tcMet_.isNonnull()) return tcMet_->phi();
    else return -9999.0;
}

const float reco::SkimEvent::chargedMet() const {
    return chargedMet_.pt();
}

/*
const float reco::SkimEvent::minMet() const {
  return ((chargedMet() < pfMet()) ? chargedMet() : pfMet()) ;
}

const math::XYZTLorentzVector reco::SkimEvent::minMetP4() const {
    return ((chargedMet() < pfMet()) ? chargedMet_.p4() : pfMet_->p4()) ;
}
*/

const float reco::SkimEvent::tcMetX() const {

    if(tcMet_.isNonnull()) return tcMet_->px();
    else return -9999.0;
}

const float reco::SkimEvent::tcMetY() const {

    if(tcMet_.isNonnull()) return tcMet_->py();
    else return -9999.0;
}

const float reco::SkimEvent::mll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0]->p4() + leps_[1]->p4()).mass();
}

const float reco::SkimEvent::pTll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0]->p4() + leps_[1]->p4()).pt();
}

const float reco::SkimEvent::dPhill() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4(),leps_[1]->p4()) );   
}

const float reco::SkimEvent::dRll() const {
  if(leps_.size()!=2) return -9999.0;
  return ROOT::Math::VectorUtil::DeltaR(leps_[0]->p4(),leps_[1]->p4());
}

const float reco::SkimEvent::dEtall() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(leps_[0]->eta() - leps_[1]->eta());
}

const float reco::SkimEvent::etall() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0]->p4() + leps_[1]->p4()).eta();
}

const float reco::SkimEvent::yll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0]->p4() + leps_[1]->p4()).Rapidity();
}

const float reco::SkimEvent::dPhillMet(metType metToUse) const {
  switch (metToUse) {
    case TCMET:  return dPhillTcMet();
    case PFMET:  return dPhillPfMet();
    case CHMET:  return dPhillChargedMet();
      //case MINMET: return dPhillMinMet();
  }
	return 0;
}

const float reco::SkimEvent::dPhillPfMet() const {
  if(leps_.size()!=2 || pfMet_.isNull()) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4()+leps_[1]->p4(),pfMet_->p4()) );
}

const float reco::SkimEvent::dPhillTcMet() const {
  if(leps_.size()!=2 || tcMet_.isNull()) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4()+leps_[1]->p4(),tcMet_->p4()) );
}

const float reco::SkimEvent::dPhillChargedMet() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4()+leps_[1]->p4(),chargedMet_.p4()) );
}

/*
const float reco::SkimEvent::dPhillMinMet() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0]->p4()+leps_[1]->p4(),minMetP4()) );
}
*/

const float reco::SkimEvent::mTHiggs(metType metToUse) const {
    // AN 2011/155, v2, p19
    return sqrt( 2 * pTll() * met(metToUse) * ( 1 - cos(dPhillMet(metToUse)) ) );
    //version 2 from guillelmo's talk
    /*
    return sqrt( mll()*mll() + 
                 2*( sqrt(pTll()*pTll()+mTll()*mTll()) * tcMet()  -
                     (pXll()+tcMetX())*(pXll()+tcMetX()) - 
                     (pYll()+tcMetY())*(pYll()+tcMetY()) ) );
    */
}

const float reco::SkimEvent::pXll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0]->p4() + leps_[1]->p4()).px();
}

const float reco::SkimEvent::pYll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0]->p4() + leps_[1]->p4()).py();
}

const float reco::SkimEvent::mTll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0]->p4() + leps_[1]->p4()).mt();
}

const float reco::SkimEvent::mT(size_t i, metType metToUse) const {
  if(i>=leps_.size()) return -9999.0;
  return sqrt(2*pt(i)*met(metToUse)*(1 - cos(dPhilMet(i, metToUse))));
}

const float reco::SkimEvent::met(metType metToUse) const {
    switch (metToUse) {
        case TCMET:  return tcMet();
        case PFMET:  return pfMet();
        case CHMET:  return chargedMet();
	  //case MINMET: return minMet();
    }
	return 0;
}

const float reco::SkimEvent::projMet(metType metToUse) const {
    float dphi = dPhilMet(metToUse);
    if(dphi < M_PI/2.) return met(metToUse)*sin(dphi);
    else               return met(metToUse);       
}

const float reco::SkimEvent::projPfMet() const {
    float dphi = dPhilPfMet();
    if(dphi < M_PI/2.) return pfMet()*sin(dphi);
    else               return pfMet();       
}

const float reco::SkimEvent::projTcMet() const {
    float dphi = dPhilTcMet();
    if(dphi < M_PI/2.) return tcMet()*sin(dphi);
    else               return tcMet();       
}

const float reco::SkimEvent::projChargedMet() const {
    float dphi = dPhilChargedMet();
    if(dphi < M_PI/2.) return chargedMet()*sin(dphi);
    else               return chargedMet();       
}

const float reco::SkimEvent::projChargedMetSmurf() const {
    float dphi = dPhilChargedMetSmurf();
    if(dphi < M_PI/2.) return chargedMetSmurf()*sin(dphi);
    else               return chargedMetSmurf();       
}


/*
const float reco::SkimEvent::projMinMet() const {
    float dphi = dPhilMinMet();
    if(dphi < M_PI/2.) return minMet()*sin(dphi);
    else               return minMet();       
}
*/

const float reco::SkimEvent::dPhilMet(metType metToUse) const {
    switch (metToUse) {
        case TCMET:  return dPhilTcMet();
        case PFMET:  return dPhilPfMet();
        case CHMET:  return dPhilChargedMet();
	  //case MINMET: return dPhilMinMet();
    }
	return 0;
}   


const float reco::SkimEvent::dPhilTcMet() const {
    float smallestDphi = 9999.;
    for(size_t l=0; l<leps_.size();++l){
        float dphi = dPhilTcMet(l);
        if( dphi < smallestDphi) smallestDphi = dphi;
    }
    return smallestDphi;
}

const float reco::SkimEvent::dPhilPfMet() const {
    float smallestDphi = 9999.;
    for(size_t l=0; l<leps_.size();++l){
        float dphi = dPhilPfMet(l);
        if( dphi < smallestDphi) smallestDphi = dphi;
    }
    return smallestDphi;
}

const float reco::SkimEvent::dPhilChargedMet() const {
    float smallestDphi = 9999.;
    for(size_t l=0; l<leps_.size();++l){
        float dphi = dPhilChargedMet(l);
        if( dphi < smallestDphi) smallestDphi = dphi;
    }
    return smallestDphi;
}

const float reco::SkimEvent::dPhilChargedMetSmurf() const {
    float smallestDphi = 9999.;
    for(size_t l=0; l<leps_.size();++l){
        float dphi = dPhilChargedMetSmurf(l);
        if( dphi < smallestDphi) smallestDphi = dphi;
    }
    return smallestDphi;
}

/*
const float reco::SkimEvent::dPhilMinMet() const {
    float smallestDphi = 9999.;
    for(size_t l=0; l<leps_.size();++l){
        float dphi = dPhilMinMet(l);
        if( dphi < smallestDphi) smallestDphi = dphi;
    }
    return smallestDphi;
}
*/

const float reco::SkimEvent::dPhilMet(size_t i, metType metToUse) const {
    if( i >= leps_.size() ) return -9999.0;
    switch (metToUse) {
        case TCMET:  return dPhilTcMet(i);
        case PFMET:  return dPhilPfMet(i);
        case CHMET:  return dPhilChargedMet(i);
	  //case MINMET: return dPhilMinMet(i);
    }
	return -9999.0;
}   

const float reco::SkimEvent::dPhilTcMet(size_t i) const {
    if( i >= leps_.size() ) return -9999.0;
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(tcMet_->p4(),leps_[i]->p4()) );
}

const float reco::SkimEvent::dPhilPfMet(size_t i) const {
    if( i >= leps_.size() ) return -9999.0;
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(pfMet_->p4(),leps_[i]->p4()) );
}

const float reco::SkimEvent::dPhilChargedMet(size_t i) const {
    if( i >= leps_.size() ) return -9999.0;
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(chargedMet_.p4(),leps_[i]->p4()) );
}

const float reco::SkimEvent::dPhilChargedMetSmurf(size_t i) const {
    if( i >= leps_.size() ) return -9999.0;
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(chargedMetSmurf_.p4(),leps_[i]->p4()) );
}

/*
const float reco::SkimEvent::dPhilMinMet(size_t i) const {
    if( i >= leps_.size() ) return -9999.0;
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(minMetP4(),leps_[i]->p4()) );
}
*/


const float reco::SkimEvent::nTracks() const {

    return 0;
}

/*
const float reco::SkimEvent::dPhiJetPfMet() const {

    return 0;
}

const float reco::SkimEvent::dPhiJetTcMet() const {

    return 0;
}
*/

const float reco::SkimEvent::cosThetaStar(size_t i) const {

    return 0;
}

const float reco::SkimEvent::afb(size_t i) const {

    return 0;
}


const bool reco::SkimEvent::leptEtaCut(float maxAbsEtaMu,float maxAbsEtaEl) const{
  bool check0(true);
  bool check1(true);

  if(abs(leps_[0]->pdgId())==11 && fabs(leps_[0]->eta())>=maxAbsEtaEl) check0=false;
  if(abs(leps_[0]->pdgId())==13 && fabs(leps_[0]->eta())>=maxAbsEtaMu) check0=false;

  if(abs(leps_[1]->pdgId())==11 && fabs(leps_[1]->eta())>=maxAbsEtaEl) check1=false;
  if(abs(leps_[1]->pdgId())==13 && fabs(leps_[1]->eta())>=maxAbsEtaMu) check1=false;

  return (check0 && check1);
}

// ... in spite my egregious programming 
void reco::SkimEvent::setTriggerBits( const std::vector<bool> &bits) {
    
    passesSingleMuData_ = bits[0];
    passesSingleElData_ = bits[1];
    passesDoubleMuData_ = bits[2];
    passesDoubleElData_ = bits[3];
    passesMuEGData_     = bits[4];
    passesSingleMuMC_   = bits[5];
    passesSingleElMC_   = bits[6];
    passesDoubleMuMC_   = bits[7];
    passesDoubleElMC_   = bits[8];
    passesMuEGMC_       = bits[9];

}

const bool reco::SkimEvent::triggerBitsCut( SkimEvent::primaryDatasetType pdType) const{

    if (pdType == MC) return true;

    if( hypo() == WWMUMU ) {
        if      ( pdType == DoubleMuon ) return ( passesDoubleMuData_ );
        else if ( pdType == SingleMuon ) return ( !passesDoubleMuData_ && passesSingleMuData_ );
        //else if ( pdType == MC         ) return ( passesDoubleMuMC_ || passesSingleMuMC_ );
    } else if( hypo() == WWMUEL || hypo() == WWELMU ) {
        if      ( pdType == SingleMuon ) return ( passesSingleMuData_ );
        else if ( pdType == MuEG       ) return ( !passesSingleMuData_ && passesMuEGData_ );
        //else if ( pdType == MC         ) return ( passesSingleMuMC_ || passesMuEGMC_ );
    } else if( hypo() == WWELEL ) {
        if      ( pdType == DoubleElectron ) return ( passesDoubleElData_ );
        //else if ( pdType == MC             ) return ( passesDoubleElMC_ );
    }

    return false;

}

const bool reco::SkimEvent::guillelmoTrigger( SkimEvent::primaryDatasetType pdType ) const {

    //Guillelmo's Implementation:
    if(pdType == MC) return true;

    if(pdType == MuEG)                  {    return (  passesMuEGData_ );
    } else if(pdType == DoubleMuon)     {    return ( !passesMuEGData_ &&  passesDoubleMuData_ );
    } else if(pdType == SingleMuon)     {    return ( !passesMuEGData_ && !passesDoubleMuData_ &&  passesSingleMuData_ );
    } else if(pdType == DoubleElectron) {    return ( !passesMuEGData_ && !passesDoubleMuData_ && !passesSingleMuData_ &&  passesDoubleElData_ );
    } else if(pdType == SingleElectron) {    return ( !passesMuEGData_ && !passesDoubleMuData_ && !passesSingleMuData_ && !passesDoubleElData_ && passesSingleElData_ );
    }

    return false;
}

bool reco::SkimEvent::passTriggerSingleMu(size_t i, bool isData) const{ 
  bool result(false);

  if( !isMuon(i) ) return false;
  pat::Muon const * const mu = getMuon(i);

  const pat::TriggerObjectStandAlone * match = mu->triggerObjectMatchByCollection("hltL3MuonCandidates");
  if(isData){
    if(match) result= (match->hasPathName("HLT_Mu24_v*",false) || match->hasPathName("HLT_IsoMu17_v*",false) );}
  else{
    if(match) result=(match->hasPathName("HLT_Mu21_v*",false) && match->pt()>24.0);
  }

  return result;
}

bool reco::SkimEvent::passTriggerDoubleMu(size_t i, bool isData) const{
  using namespace std;
  bool result(false);
  
  if( !isMuon(i) ) return false;
  pat::Muon const * const mu = getMuon(i);

  const pat::TriggerObjectStandAlone * match = mu->triggerObjectMatchByCollection("hltL3MuonCandidates");
  if(isData){
    if(match) result=match->hasPathName("HLT_DoubleMu7_v*",false);  }
  else{
    if(match) result=(match->hasPathName("HLT_DoubleMu5_v*",false) && match->pt()>7.0);  
  }
  return result;
}

bool reco::SkimEvent::passTriggerDoubleEl(size_t i, bool isData) const{ 
  bool result(false);
  
  if( !isElectron(i) ) return false;
  pat::Electron const * const el = getElectron(i);
    
  if(isData){  
    if(el->triggerObjectMatchesByPath("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*").size() ||
       el->triggerObjectMatchesByPath("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*").size() ) 
      result=true;
  }else{
    if(el->triggerObjectMatchesByPath("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*").size() )
      result=true;
  }
  
  return result;

}
bool reco::SkimEvent::passTriggerElMu(size_t i, bool isData) const{ 
  /*
    vector< std::string > pathNames = match->pathNames(false);
    for(unsigned int i=0; i<pathNames.size();++i){
    cout << "match path name: " << pathNames[i] << endl;
    }
  */
  bool result(false);
  using namespace std;
  if( isMuon(i) ) {
    pat::Muon const * const mu = getMuon(i);
    const pat::TriggerObjectStandAlone * match = mu->triggerObjectMatchByCollection("hltL3MuonCandidates");
    
    if(match){
      if(isData){
	bool res1 = match->hasPathName("HLT_Mu8_Ele17_CaloIdL_v*",false);
	bool res2 = match->hasPathName("HLT_Mu17_Ele8_CaloIdL_v*",false);
	result=( res1 || res2);}
      else{
	bool res1 = (match->hasPathName("HLT_Mu5_Ele17_v*",false) && match->pt()>8.0);
	bool res2 = (match->hasPathName("HLT_Mu11_Ele8_v*",false) && match->pt()>17.0);
	result=( res1 || res2);
      }    
    }
    return result;
  }

  if( isElectron(i) ) {
    pat::Electron const * const el = getElectron(i);
    if(isData){
      if(el->triggerObjectMatchesByPath("HLT_Mu8_Ele17_CaloIdL*").size() ||
	 el->triggerObjectMatchesByPath("HLT_Mu17_Ele8_CaloIdL*").size() ) 
	result=true;}
    else{
      const pat::TriggerObjectStandAlone * match1=
	el->triggerObjectMatchByPath("HLT_Mu5_Ele17_v*",true);
      const pat::TriggerObjectStandAlone * match2=
	el->triggerObjectMatchByPath("HLT_Mu11_Ele8_v*",true);     
      result=( match1 || match2 );
    }        
  }

  return result;

}


// TO BE FIXED: This guy should take the list of trigger from some 
// sort of configuration file. 
const bool reco::SkimEvent::triggerMatchingCut(SkimEvent::primaryDatasetType pdType) const{
  if( pdType == MC ) return true;

  using namespace std;
  bool result(false);  

  if(hypo()==WWMUMU){
    if(pdType==DoubleMuon){ //configuration (1)
      result=(passTriggerDoubleMu(0) && passTriggerDoubleMu(1));}
    if(pdType==SingleMuon) {//configuration (2)
      result=(  (passTriggerSingleMu(0) || passTriggerSingleMu(1)) &&
		!(passTriggerDoubleMu(0) && passTriggerDoubleMu(1)) );
    }
    if(pdType==MC){
      result=( (passTriggerDoubleMu(0,false) && passTriggerDoubleMu(1,false)) ||
	       (passTriggerSingleMu(0,false) || passTriggerSingleMu(1,false)) );
    }
  }

  if(hypo()==WWMUEL || hypo()==WWELMU){
    if(pdType==SingleMuon) //configuration (3)
      result=(  passTriggerSingleMu(0) || passTriggerSingleMu(1) );
    if(pdType==MuEG){       //configuration (4)
      result=( (passTriggerElMu(0) && passTriggerElMu(1)) && 
	       !(passTriggerSingleMu(0) || passTriggerSingleMu(1)) );    
    }
    if(pdType==MC){
      result=( (passTriggerSingleMu(0,false) || passTriggerSingleMu(1,false)) ||
	       (passTriggerElMu(0,false) && passTriggerElMu(1,false))  );	       
    }
  }


  if(hypo()==WWELEL){
    if(pdType==DoubleElectron)//configuration (5)
      result= (passTriggerDoubleEl(0) && passTriggerDoubleEl(1));
    if(pdType==MC)
      result= (passTriggerDoubleEl(0,false) && passTriggerDoubleEl(1,false));
  }

  return result;
}


//Iso Functions
const float reco::SkimEvent::tkIso(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return getElectron(i)->dr03TkSumPt();
    } else if ( isMuon(i) ) {
        return getMuon(i)->isolationR03().sumPt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::ecalIso(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return getElectron(i)->dr03EcalRecHitSumEt();
    } else if ( isMuon(i) ) {
        return getMuon(i)->isolationR03().emEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::hcalIso(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return getElectron(i)->userFloat("hcalFull");
    } else if ( isMuon(i) ) {
        return getMuon(i)->isolationR03().hadEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::getRho(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;
    
    if( isElectron(i) ) {
        return getElectron(i)->userFloat("rhoEl");
    } else if ( isMuon(i) ) {
        return getMuon(i)->userFloat("rhoMu");
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::allIso(size_t i) const {
    
    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return getElectron(i)->userFloat("eleSmurfPF04");
    } else if( isMuon(i) ) {
        return getMuon(i)->userFloat("muSmurfPF");
    } else {
        std::cout << " Do I ever friggin get here?" << std::endl;
        return 9999.0;
    }
//     if( isElectron(i) && isEB(i) ) {
//         return tkIso(i) + std::max((float)0,ecalIso(i)-1) + hcalIso(i) - getRho(i) * 3.14159265 * 0.3 * 0.3;
//     } else if( (isElectron(i) && !isEB(i)) || isMuon(i) ) {
//         return tkIso(i) + ecalIso(i) + hcalIso(i) - getRho(i) * 3.14159265 * 0.3 * 0.3;
//     } else {
//         std::cout << " Do I ever friggin get here?" << std::endl;
//         return 9999.0;
//     }
}

const float reco::SkimEvent::tkVeto(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return 0;
    } else if ( isMuon(i) ) {
        return getMuon(i)->isolationR03().trackerVetoPt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::ecalVeto(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return 0;
    } else if ( isMuon(i) ) {
        return getMuon(i)->isolationR03().emVetoEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::hcalVeto(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return 0;
    } else if ( isMuon(i) ) {
        return getMuon(i)->isolationR03().hadVetoEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::allVeto(size_t i) const {
    
    if( i >= leps_.size() ) return -9999.0;

    if( isElectron(i) ) {
        return 0;
    } else if( isMuon(i) ) {
        return tkVeto(i) + ecalVeto(i) + hcalVeto(i);
    } else {
        return -9999.0;
    }
}

const size_t reco::SkimEvent::indexByPt(size_t i) const {

    if( i >= leps_.size() ) return 9999;
    std::vector<indexValueStruct> a;

    for(size_t j=0;j<leps_.size();++j) a.push_back(indexValueStruct(pt(j),j));
    std::sort(a.begin(),a.end(),highToLow);

    return a[i].index;
}

const size_t reco::SkimEvent::indexByIso(size_t i) const {
    
    if( i >= leps_.size() ) return 9999;
    std::vector<indexValueStruct> a;

    for(size_t j=0;j<leps_.size();++j) a.push_back(indexValueStruct(allIso(j),j));
    std::sort(a.begin(),a.end(),lowToHigh);

    return a[i].index;
}

const bool reco::SkimEvent::isEB(size_t i) const {

    if( isElectron(i) ) {
        return getElectron(i)->isEB();
    } else {
        return false;
    }
}

const bool reco::SkimEvent::isEE(size_t i) const {

    if( isElectron(i) ) {
        return getElectron(i)->isEE();
    } else {
        return false;
    }
}

const float reco::SkimEvent::tkPt(size_t i) const {

    if( isElectron(i) ) {
        return getElectron(i)->gsfTrack()->pt();
    } else if( isMuon(i) && !isSTA(i)) {
        return getMuon(i)->track()->pt();
    } else {
        return -9999.9;
    }

}

const bool reco::SkimEvent::passesVtxSel(size_t i) const {
    if(i >= vtxs_.size()) {
        return false;
    } else {
        return ( vtxs_[i]->isValid() && !vtxs_[i]->isFake()         && 
                 vtxs_[i]->ndof()                            >= 4.0 && 
            fabs(vtxs_[i]->position().Rho())                  < 2.0 &&
            fabs(vtxs_[i]->z())                               < 24.0  );
    }
    return false;
}

const int reco::SkimEvent::nGoodVertices() const {

    int count = 0;
    for(size_t i=0;i<vtxs_.size();++i) 
        if( passesVtxSel(i) ) count++;

    return count;
}

const bool reco::SkimEvent::hasGoodVertex() const {

    return (nGoodVertices() > 0);
}

const reco::Vertex reco::SkimEvent::highestPtVtx() const {
    if(vtxs_.size() == 0) return reco::Vertex();
    if(sumPts_.size() == 0) return *vtxs_[0];
    double sum = 0;
    size_t high = 0;
    for(size_t i=0;i<vtxs_.size();++i) {
        if( sumPts_[i] > sum && passesVtxSel(i) ) {
            high = i;
            sum = sumPts_[i];
        }
    }
    return *vtxs_[high];
}

const bool reco::SkimEvent::passesIP() const {
    return (passesIP(leps_[0]) && passesIP(leps_[1]));
}

const bool reco::SkimEvent::passesIP(const refToCand &c) const {

    if( isElectron(c) ) {

        if( fabs(getElectron(c)->userFloat("ip2")) >= 0.03) return false;

    } else if( isMuon(c) && !isSTA(c) ) {

        if( fabs(getMuon(c)->userFloat("tip2")) >= 0.01) return false;
        if( fabs(getMuon(c)->userFloat("dzPV")) >= 0.05) return false;

    }
    
    return true;

}

const double reco::SkimEvent::d0Reco(size_t i) const {
  double dxyPV = 9999;
  if( isElectron(i) ) {
    dxyPV = getElectron(i)->userFloat("dxyPV");
  } else if( isMuon(i) && !isSTA(i) ) {
    dxyPV = getMuon(i)->userFloat("dxyPV");
  } 
  return dxyPV;
}

const double reco::SkimEvent::dZReco(size_t i) const {
  double dzPV = 9999;
  if( isElectron(i) ) {
    dzPV = getElectron(i)->userFloat("dzPV");
  } else if( isMuon(i) && !isSTA(i) ) {
    dzPV = getMuon(i)->userFloat("dzPV");
  } 
  return dzPV;
}

const bool reco::SkimEvent::passesConversion(size_t i) const {
  if( isElectron(i) ) {
    pat::Electron const * const e = getElectron(i);
    if( fabs(e->userFloat("convValueMapProd:dist")) < 0.02 &&
	fabs(e->userFloat("convValueMapProd:dcot")) < 0.02 ) {
      return false;
    }

    if( e->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() > 0 ) {
      return false;
    }
    return true;
  } else if( isMuon(i) ) {
    return true;
  } else {
    return false;
  }

}

const bool reco::SkimEvent::isSTA(size_t i) const {
  return isSTA(leps_[i]);
}

const bool reco::SkimEvent::isSTA(const refToCand &c) const {
  if( isMuon(c) ) {
    return (getMuon(c)->type() == 8);
  } else {
    return false;
  }
}

const float reco::SkimEvent::highestHardBDisc(const float& maxPt, std::string discriminator, int applyID, float dzCut) const {

    float disc=-9999.9;

    for(size_t i=0;i<tagJets_.size();++i) {      
        if( tagJetPt(i,true) < maxPt ) continue;
        if(!(passJetID(tagJets_[i],applyID)) ) continue;
        if(isThisJetALepton(tagJets_[i])) continue;
        if(jets_[i]->hasUserFloat("dz") && fabs(jets_[i]->userFloat("dz")) > dzCut) continue;
        if( tagJets_[i]->bDiscriminator(discriminator) > disc ) disc = tagJets_[i]->bDiscriminator(discriminator);
    }

    return disc;

}

const float reco::SkimEvent::highestBDiscRange(const float& minPt, const float& maxPt, std::string discriminator, int applyID, float dzCut) const {

    float disc=-9999.9;

    for(size_t i=0;i<tagJets_.size();++i) {      
        if( tagJetPt(i,true) > maxPt ) continue;
        if( tagJetPt(i,false) <= minPt ) continue;
        if(!(passJetID(tagJets_[i],applyID)) ) continue;
        if(isThisJetALepton(tagJets_[i])) continue;
        if(jets_[i]->hasUserFloat("dz") && fabs(jets_[i]->userFloat("dz")) > dzCut) continue;
        if( tagJets_[i]->bDiscriminator(discriminator) > disc ) disc = tagJets_[i]->bDiscriminator(discriminator);
    }

    return disc;

}

const float reco::SkimEvent::highestSoftBDisc(const float& maxPt, std::string discriminator, int applyID, float dzCut) const {

    float disc=-9999.9;

    for(size_t i=0;i<tagJets_.size();++i) {      
        if( tagJetPt(i,true) > maxPt ) continue;
        if(!(passJetID(tagJets_[i],applyID)) ) continue;
        if(isThisJetALepton(tagJets_[i])) continue;
        if(jets_[i]->hasUserFloat("dz") && fabs(jets_[i]->userFloat("dz")) > dzCut) continue;
        if( tagJets_[i]->bDiscriminator(discriminator) > disc ) disc = tagJets_[i]->bDiscriminator(discriminator);
    }

    return disc;

}

const int reco::SkimEvent::bTaggedJetsBetween(const float& minPt, const float& maxPt, const float& cut, std::string discriminator, int applyID , float dzCut) const {

    int count=0;

    for(size_t i=0;i<tagJets_.size();++i) {     
        if( tagJetPt(i,true) > maxPt ) continue;
        if( tagJetPt(i,false) <= minPt ) continue;
        if(!(passJetID(tagJets_[i],applyID)) ) continue;
        if( tagJets_[i]->bDiscriminator(discriminator) <= cut ) continue;	
        if(jets_[i]->hasUserFloat("dz") && fabs(jets_[i]->userFloat("dz")) > dzCut) continue;
        if(isThisJetALepton(tagJets_[i])) continue;
        count++;
    }

    return count;
}

const int reco::SkimEvent::bTaggedJetsUnder(const float& maxPt, const float& cut, std::string discriminator, int applyID, float dzCut) const {

    int count=0;

    for(size_t i=0;i<tagJets_.size();++i) {     
        if( tagJetPt(i,true) > maxPt ) continue;
        if(!(passJetID(tagJets_[i],applyID)) ) continue;
        if( tagJets_[i]->bDiscriminator(discriminator) <= cut ) continue;	
        if(jets_[i]->hasUserFloat("dz") && fabs(jets_[i]->userFloat("dz")) > dzCut) continue;
        if(isThisJetALepton(tagJets_[i])) continue;
        count++;
    }

    return count;
}

const int reco::SkimEvent::bTaggedJetsOver(const float& maxPt, const float& cut, std::string discriminator, int applyID, float dzCut) const {

    int count=0;

    for(size_t i=0;i<tagJets_.size();++i) {
        if( tagJetPt(i,true) <= maxPt ) continue;
        if(!(passJetID(tagJets_[i],applyID)) ) continue;
        if( tagJets_[i]->bDiscriminator(discriminator) <= cut ) continue;
        if(jets_[i]->hasUserFloat("dz") && fabs(jets_[i]->userFloat("dz")) > dzCut) continue;
        if(isThisJetALepton(tagJets_[i])) continue;
        count++;
    }

    return count;
}

const bool reco::SkimEvent::isMuTriggered(size_t i) const {

    if( isMuon(i) ) {
        pat::Muon const * const mu = getMuon(i);
        return ( !mu->triggerObjectMatchesByPath("HLT_Mu9").empty() || !mu->triggerObjectMatchesByPath("HLT_Mu15_v1").empty() );
    } else {
        return false;
    }


}


const int reco::SkimEvent::mitType() const {

    if( abs(pdgIdByPt(0)) == 11 ) {
        if( abs(pdgIdByPt(1)) == 11 ) return 11;
        if( abs(pdgIdByPt(1)) == 13 ) return 12;
    } else if( abs(pdgIdByPt(0)) == 13 ) {
        if( abs(pdgIdByPt(1)) == 11 ) return 13;
        if( abs(pdgIdByPt(1)) == 13 ) return 10;
    } else {
        return -1;
    }
	return -1;
}


const float reco::SkimEvent::nearestJet(int i,float minPt, float eta, bool applyCorrection , int applyID ) const {

    if (i >= (int)leps_.size() || i < -1) return -9999.9;

    float dR = 9999;
    for(size_t j=0;j<jets_.size();++j) {
        if(!(passJetID(jets_[i],applyID)) ) continue;
        if( std::fabs(jets_[j]->eta()) >= eta) continue;
        if( jetPt(j,applyCorrection) <= minPt) continue;
        if(isThisJetALepton(jets_[j]))  continue;
    
        float tempdR;
        if(i != -1) {
            tempdR = fabs(ROOT::Math::VectorUtil::DeltaR(jets_[j]->p4(),leps_[i]->p4()) );
            if( tempdR < dR ) dR = tempdR;
        } else {
            for(size_t k=0; k<leps_.size();++k){
               tempdR = fabs(ROOT::Math::VectorUtil::DeltaR(jets_[j]->p4(),leps_[k]->p4()) );
               if( tempdR < dR ) dR = tempdR;
            }
        }
    }
    return dR;
}

const pat::JetRef reco::SkimEvent::matchedJet(size_t i, float minDr) const {
    pat::JetRef ret;
    if (i >= leps_.size()) return ret;

    float dR = minDr;
    for(size_t j=0;j<jets_.size();++j) {
        float tempdR = fabs(ROOT::Math::VectorUtil::DeltaR(jets_[j]->p4(),leps_[i]->p4()) );
        if( tempdR < dR ) { 
            dR = tempdR;
            ret = jets_[j];
        }
    }
    return ret;
}

const float reco::SkimEvent::matchedJetPt(size_t i, float minDr, bool applyCorrection) const {
    if (i >= leps_.size()) return -9999.9;

    float dR = minDr, pt = 0;
    for(size_t j=0;j<jets_.size();++j) {
        float tempdR = fabs(ROOT::Math::VectorUtil::DeltaR(jets_[j]->p4(),leps_[i]->p4()) );
        if( tempdR < dR ) { 
            dR = tempdR;
            pt = jetPt(j,applyCorrection);
        }
    }
    return pt;
}



// ============== Matt's methods ===============================
// const int reco::SkimEvent::nExtraLepMatt(float minPt) const { 
//     int count = 0;
//     for(size_t i=0;i<extraLeps_.size();++i) {
//         if( extraLeps_[i]->pt() <= minPt) continue;
//         if( fabs(extraLeps_[i]->pdgId()) == 11 ) {
//             pat::Electron e = static_cast<const pat::Electron&>(extraLeps_[i]);
//             if( fabs(e.eta()) >= 2.5 ) continue;
//             if( fabs(e.userFloat("dxyPV")) >= 0.020 ) continue;
//             if( fabs(e.userFloat("dzPV") ) >= 1.0 ) continue;
//             if( fabs(e.userFloat("convValueMapProd:dist")) < 0.2 &&
//                 fabs(e.userFloat("convValueMapProd:dcot")) < 0.2 ) continue;
//             if( e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() > 0 ) continue;
//             if(!(( e.isEB() && e.sigmaIetaIeta() < 0.01 &&
//                   fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.06 &&
//                   fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.004 &&
// //                   e.hadronicOverEm() < 0.04 && 
//                   (e.dr03TkSumPt() + std::max(e.dr03EcalRecHitSumEt()-1,(float)0) + e.dr03HcalTowerSumEt() - e.userFloat("rhoEl")*3.14159*0.3*0.3)/e.pt() < 0.1) ||
//                   ( !e.isEB() && e.sigmaIetaIeta() < 0.03  && 
//                   fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.03 &&
//                   fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.007 &&
//                   (e.dr03TkSumPt() + e.dr03EcalRecHitSumEt() + e.dr03HcalTowerSumEt() - e.userFloat("rhoEl")*3.14159265*0.3*0.3)/e.pt() < 0.1 //&&
//                   /*e.hadronicOverEm() < 0.025 */))) continue;
//         } else if ( fabs(extraLeps_[i]->pdgId()) == 13 ) {
//             pat::Muon m = static_cast<const pat::Muon&>(extraLeps_[i]);
//             if( fabs(m.userFloat("dxyPV")) >= 0.020 ) continue;
//             if( fabs(m.userFloat("dzPV") ) >= 1.0 ) continue;
//             if( !(m.isGlobalMuon() && m.isTrackerMuon() && 
//                   m.innerTrack()->found() > 10 &&
//                   m.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
//                   m.globalTrack()->normalizedChi2() < 10 &&
//                   m.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
//                   m.numberOfMatches() > 1 && fabs(m.track()->ptError() / m.pt()) < 0.10 )) continue;
//             if( (m.isolationR03().emEt  + m.isolationR03().hadEt + m.isolationR03().sumPt - m.userFloat("rhoMu") * 3.14159265 * 0.3 * 0.3) / m.pt() >= 0.15 ) continue;
//         }
//         count++;
//     }
//     return count;
// } 
// const int reco::SkimEvent::nSoftMuMatt(float minPt) const { 
//   int count = 0;
//   for(size_t i=0;i<softMuons_.size();++i) {
//     if(softMuons_[i].pt() <= minPt ) continue;
//     if( (softMuons_[i].isolationR03().emEt  +
// 	 softMuons_[i].isolationR03().hadEt +
// 	 softMuons_[i].isolationR03().sumPt ) / softMuons_[i].pt() < 0.1 &&
// 	softMuons_[i].pt() > 20.) continue;
//     count++;
//   }
//   return count;
// } 
// const bool reco::SkimEvent::passesIDV1(size_t i) const {
// 
//     if( fabs(leps_[i]->pdgId()) == 11 ) {
//         const pat::Electron & e = static_cast<const pat::Electron&>(leps_[i]);
//         return (( e.isEB() && e.sigmaIetaIeta() < 0.01 &&
//                 fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.06 &&
//                 fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.004 //&&
//                 /*e.hadronicOverEm() < 0.04*/) ||
//                 ( !e.isEB() && e.sigmaIetaIeta() < 0.03  && 
//                 fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.03 &&
//                 fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.007 //&&
//                 /*e.hadronicOverEm() < 0.025*/ ));
//     } else if( fabs(leps_[i]->pdgId()) == 13 ) {
//         const pat::Muon & m = static_cast<const pat::Muon&>(leps_[i]);
//         return (m.isGlobalMuon() && m.isTrackerMuon() && 
//                 m.innerTrack()->found() > 10 &&
//                 m.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
// //                 ( (float)m.innerTrack()->hitPattern().trackerLayersWithMeasurement() / (float)(
// //                     m.innerTrack()->hitPattern().trackerLayersWithoutMeasurement() + 
// //                     m.innerTrack()->hitPattern().trackerLayersWithMeasurement() + 
// //                     m.innerTrack()->trackerExpectedHitsInner().numberOfLostHits() +
// //                     m.innerTrack()->trackerExpectedHitsOuter().numberOfLostHits() ) >= 0.75 ) &&
// //                 ( (float)m.innerTrack()->hitPattern().numberOfValidHits() /  (float)
// //                     (m.innerTrack()->hitPattern().numberOfHits() - 
// //                      m.innerTrack()->hitPattern().numberOfInactiveHits()) > 0.92 ) &&
//                 m.globalTrack()->normalizedChi2() < 10 &&
//                 m.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
//                 m.numberOfMatches() > 1 && fabs(m.track()->ptError() / m.pt()) < 0.10 );
//     } else {
//         return false;
//     }
// 
// }



// ============== Boris' shit ===============================
// Sorry dude, commented out temporarily
// const bool reco::SkimEvent::passMuID0() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       answer = answer && (mu.genParticleRef().isNonnull() && fabs(mu.genParticleRef()->pdgId())==13);
//     }
//   }
//   return answer;
// }
// 
// 
// const bool reco::SkimEvent::passMuID1() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       answer = answer && (mu.isGlobalMuon() || (mu.isTrackerMuon() && mu.muonID("TMLastStationTight")));
//     }
//   }
//   return answer;
// }
// 
// const bool reco::SkimEvent::passMuID2() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       answer = answer && (mu.isTrackerMuon());
//     }
//   }
//   return answer;
// }
// 
// const bool reco::SkimEvent::passMuID3() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       answer = answer && (mu.innerTrack()->found()>10);
//     }
//   }
//   return answer;
// }
// 
// const bool reco::SkimEvent::passMuID4() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       answer = answer && (mu.innerTrack()->hitPattern().numberOfValidPixelHits()>0);
//     }
//   }
//   return answer;
// }
// 
// const bool reco::SkimEvent::passMuID5() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       if(mu.isGlobalMuon())
// 	answer = answer && (mu.globalTrack()->normalizedChi2()<10);
//     }
//   }
//   return answer;
// }
// 
// const bool reco::SkimEvent::passMuID6() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       if(mu.isGlobalMuon())
// 	answer = answer && (mu.globalTrack()->hitPattern().numberOfValidMuonHits()>0);
//     }
//   }
//   return answer;
// }
// 
// const bool reco::SkimEvent::passMuID7() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       if(!mu.isGlobalMuon())
// 	answer = answer && (mu.numberOfMatches()>1);
//     }
//   }
//   return answer;
// }
// 
// const bool reco::SkimEvent::passMuID8() const {
//   bool answer=true;
//   for(unsigned int i=0; i<=1; i++){
//     if( fabs(leps_[i]->pdgId()) == 13 ) {
//       const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
//       answer = answer && (mu.track()->ptError()/mu.pt()<0.10);
//     }
//   }
//   return answer;
// }

// New emanuele gamma mr star thingy
const float reco::SkimEvent::mRStar() const {
  float A = leps_[0]->p();
  float B = leps_[1]->p();
  float az = leps_[0]->pz();
  float bz = leps_[1]->pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(leps_[0]->px(),leps_[0]->py(),0.0);
  jbT.SetXYZ(leps_[1]->px(),leps_[1]->py(),0.0);

  float temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  return temp;
}

const float reco::SkimEvent::gamma() const {
  float A = leps_[0]->p();
  float B = leps_[1]->p();
  float az = leps_[0]->pz();
  float bz = leps_[1]->pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(leps_[0]->px(),leps_[0]->py(),0.0);
  jbT.SetXYZ(leps_[1]->px(),leps_[1]->py(),0.0);

  float ATBT = (jaT+jbT).Mag2();
  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));

  return 1./sqrt(1.-mybeta*mybeta);
}

const float reco::SkimEvent::gammaMRStar() const {

    return gamma() * mRStar();

}

const int reco::SkimEvent::nJetVBF(float minPt,float eta,int applyCorrection, int applyID) const {
    if(jets_.size() < 2) return -9999.;
    float jetEta1 = leadingJetEta(0,minPt,eta,applyCorrection,applyID);
    float jetEta2 = leadingJetEta(1,minPt,eta,applyCorrection,applyID);
    float jetEtaMax = std::max(jetEta1,jetEta2);
    float jetEtaMin = std::min(jetEta1,jetEta2);

    int count = 0;
    for(size_t i=0;i<jets_.size();++i) {
      if(!(passJetID(jets_[i],applyID)) ) continue;
      if( std::fabs(jets_[i]->eta()) >= eta) continue;
      if( jetPt(i,applyCorrection) <= minPt) continue;

      if(isThisJetALepton(jets_[i]))  continue;
      if( jets_[i]->eta() >= jetEtaMax ) continue;
      if( jets_[i]->eta() <= jetEtaMin ) continue;
      count++;
    }
    return count;

}

const float reco::SkimEvent::mjj(float pt ,float eta,int applyCorrection, int applyID) const {
    if(jets_.size() < 2) return -9999.;
    return (leadingJet(0,pt,eta,applyCorrection,applyID)->p4() + 
            leadingJet(1,pt,eta,applyCorrection,applyID)->p4()).mass();

}

const float reco::SkimEvent::dEtajj(float pt ,float eta,int applyCorrection, int applyID) const {
    if(jets_.size() < 2) return -9999.;
    return fabs(leadingJetEta(0,pt,eta,applyCorrection,applyID) - 
                leadingJetEta(1,pt,eta,applyCorrection,applyID));

}

const float reco::SkimEvent::zeppenfeld(size_t a,float pt ,float eta,int applyCorrection, int applyID) const {
    if(jets_.size() < 2) return -9999.;
    return -9999.;

}



void reco::SkimEvent::FindDaughterParticles(const reco::Candidate** pCurrent, std::vector<const reco::Candidate*>* pFinal) const {
  
  // Variable for decayed current particle
  const reco::Candidate* pCurrentNew = 0;
  
  for(unsigned int i = 0; i < (*pCurrent) -> numberOfDaughters(); ++i) {
    if( (pFinal != 0) && ((*pCurrent) -> daughter(i) -> pdgId() != (*pCurrent) -> pdgId()) )
      pFinal -> push_back((*pCurrent) -> daughter(i));
    
    if((*pCurrent) -> daughter(i) -> pdgId() == (*pCurrent) -> pdgId())
      pCurrentNew = (*pCurrent) -> daughter(i);
  }
  
  // Change the current particle into the decayed one
  if(pCurrentNew)
    (*pCurrent) = pCurrentNew;
}





const float reco::SkimEvent::getFinalStateMC() const {

//  std::cout << " getFinalStateMC " << std::endl;
  float finalState = -1;
  // 0 = mm
  // 1 = ee
  // 2 = tt

  const reco::Candidate* mcV = 0;
  const reco::Candidate* mcF1_fromV;
  const reco::Candidate* mcF2_fromV;
  // loop over gen particles
  for(size_t gp=0; gp<genParticles_.size();++gp){
    const reco::Candidate* pMother = 0;
    if(genParticles_[gp] -> mother()) {
      pMother = genParticles_[gp] -> mother();
    }

    int pdgId  = genParticles_[gp] -> pdgId();
    int status = genParticles_[gp] -> status();
    int charge = genParticles_[gp] -> charge();
    int motherPdgId = 0;
    if(genParticles_[gp] -> mother()) {
      motherPdgId = pMother -> pdgId();
    }

    // Z {23}
    if( (pdgId == 23) && (status == 3) ) {
      mcV = &(*(genParticles_[gp]));
    }
    
  } // loop over gen particles
    
  if (mcV != 0) {
  // find fermions from vector boson decay
   std::vector<const reco::Candidate*> fFromVBuffer;
   FindDaughterParticles(&mcV, &fFromVBuffer);
    
   if (fFromVBuffer.size() == 2) {
    mcF1_fromV = fFromVBuffer.at(0); 
    mcF2_fromV = fFromVBuffer.at(1);
  
   // If leptons, see if there is a photon emission
    if(abs(mcF1_fromV -> pdgId()) >= 11) {
     FindDaughterParticles(&mcF1_fromV);
    }
    if(abs(mcF2_fromV -> pdgId()) >= 11) {
     FindDaughterParticles(&mcF2_fromV);
    }

   // mm
    if ( abs(mcF1_fromV -> pdgId()) == 13 && abs(mcF2_fromV -> pdgId()) == 13 ) {
     finalState = 0;
    }
   
   // ee
    if ( abs(mcF1_fromV -> pdgId()) == 11 && abs(mcF2_fromV -> pdgId()) == 11 ) {
     finalState = 1;
    }
 
   // tt
    if ( abs(mcF1_fromV -> pdgId()) == 15 && abs(mcF2_fromV -> pdgId()) == 15 ) {
     finalState = 2;
    }
   }
  }

  return finalState;
}









