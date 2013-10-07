#ifndef AnalysisDataFormats_SkimEvent2L_h
#define AnalysisDataFormats_SkimEvent2L_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <string>
#include <vector>
#include <cmath>
#include <map>
#include "WWAnalysis/DataFormats/interface/hzziso.h"

namespace reco {

   class SkimEvent2L : public pat::CompositeCandidate {


        public:
            SkimEvent2L() {}
            virtual ~SkimEvent2L() ;
            explicit SkimEvent2L(const pat::CompositeCandidate &src) ;
            explicit SkimEvent2L(const reco::CompositeCandidate &src) ;

            const reco::Candidate & l(unsigned int il) const {
                const reco::Candidate &proxy = * daughter(il);
                return  proxy.hasMasterClonePtr() ? * proxy.masterClonePtr() :
                       (proxy.hasMasterClone()    ? * proxy.masterClone() :
                        proxy);
            }

            float luserFloat(unsigned int il, const std::string &label)  const {
                const reco::Candidate &c = l(il);
                if      (typeid(c) == typeid(pat::Muon))     return (dynamic_cast<const pat::Muon &>(c)).userFloat(label);
                else if (typeid(c) == typeid(pat::Electron)) return (dynamic_cast<const pat::Electron &>(c)).userFloat(label);
                else throw cms::Exception("WrongType") << "Lepton " << il << " is of type " << typeid(c).name() << "\n";
            }
            int luserInt(unsigned int il, const std::string &label)  const {
                const reco::Candidate &c = l(il);
                if      (typeid(c) == typeid(pat::Muon))     return (dynamic_cast<const pat::Muon &>(c)).userInt(label);
                else if (typeid(c) == typeid(pat::Electron)) return (dynamic_cast<const pat::Electron &>(c)).userInt(label);
                else throw cms::Exception("WrongType") << "Lepton " << il << " is of type " << typeid(c).name() << "\n";
            }



            float lval(unsigned int il, const std::string &muExpr, const std::string &eleExpr) const;
            float lval(unsigned int il, const std::string &expr) const;

            unsigned int numvertices() const { return vertices_.size(); }
            unsigned int numsimvertices() const { return numsimvertices_; }
            float numTrueInteractions() const {return numTrueInteractions_;}
            float met() const { return pfMet_.et(); }
            void setVertices(const edm::Handle<reco::VertexCollection> &);
            void setPFMet(const edm::Handle<reco::PFMETCollection> &);
            void setPileupInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &);

            /// Di-lepton mass, even if there's a FSR photon
            const double mll() const { return (daughter(0)->p4() + daughter(1)->p4()).M(); }

            float rho() const {
                return lval(0, "userFloat(\"rhoMu\")", "userFloat(\"rhoEl\")");
            }

            float getEcalIso(unsigned int l) const {
                if (l != 0 && l != 1) return false;
                float ecaliso = luserFloat(l, "ecalZZ4L");
                return HZZIso::getEcalIso(ecaliso, rho(), daughter(l)->pt(), fabs(daughter(l)->eta()), abs(daughter(l)->pdgId()) == 13);
            }

            float getHcalIso(unsigned int l) const {
                if (l != 0 && l != 1) return false;
                float hcaliso = luserFloat(l, "hcalZZ4L");
                return HZZIso::getHcalIso(hcaliso, rho(), daughter(l)->pt(), fabs(daughter(l)->eta()), abs(daughter(l)->pdgId()) == 13);
            }

            float combinedIso(unsigned int l) const {
                if (l != 0 && l != 1) return false;
                float tkiso = luserFloat(l, "tkZZ4L");
                return tkiso + getEcalIso(l) + getHcalIso(l);
            }
            
            float combinedPairRelativeIso() const {
                return (combinedIso(0) / daughter(0)->pt()) + (combinedIso(1) / daughter(1)->pt());
            }

            bool passTrig(unsigned int l) const {
                if (l != 0 && l != 1) return false;
                std::string mustr = "triggerObjectMatchesByPath(\"HLT_Mu17_Mu8_v*\").size()";
                std::string elstr = "triggerObjectMatchesByPath(\"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*\").size()";
                return lval(l, mustr, elstr) > 0;
            }
            
            bool passID(unsigned int l) const {
                if (l != 0 && l != 1) return false;
                bool isMu = abs(daughter(l)->pdgId()) == 13;
                if (!isMu) {
                    float elid = lval(l, "electronID(\"cicTight\")");
                    bool passcic = (elid == 1 || elid == 3 || elid == 5 || elid == 7 || elid == 9 || elid == 11 || elid == 13 || elid == 15);
                    bool passconv = lval(l, "gsfTrack.trackerExpectedHitsInner.numberOfHits")<=1;
                    return (passcic && passconv);
                }
                else {
                    return (lval(l, "isGlobalMuon")>0 && lval(l, "track.numberOfValidHits")>10);
                }
        
            }            

            bool passIP(unsigned int l) const {
                if (l != 0 && l != 1) return false;
                return (luserFloat(l, "ip"))/(luserFloat(l, "ipErr")) < 4;
            }

            reco::GenParticleRef genz() const ;
            reco::GenParticleRef genl(int il) const { return (unsigned(il) < matches_.size() ? matches_[il] : reco::GenParticleRef()); }
            reco::Particle::LorentzVector genllp4() const ;

            void setGenMatches(const edm::Association<reco::GenParticleCollection> &genMatch) ;
        private:
            std::vector<reco::Vertex>  vertices_;
            reco::PFMET   pfMet_;
            unsigned int numsimvertices_;
	    float numTrueInteractions_;
            std::vector<reco::GenParticleRef> matches_; 
    };

}

#endif

