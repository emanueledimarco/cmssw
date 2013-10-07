//
// $Id: L1Extra2PAT.cc,v 1.4 2011/03/31 09:59:33 gpetrucc Exp $
//

/**
  \class    pat::L1Extra2PAT L1Extra2PAT.h "MuonAnalysis/MuonAssociators/interface/L1Extra2PAT.h"
  \brief    Simple copy of l1extra to pat::TriggerObjectStandAlone
            
  \author   Giovanni Petrucciani
  \version  $Id: L1Extra2PAT.cc,v 1.4 2011/03/31 09:59:33 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace pat {

  class L1Extra2PAT : public edm::EDProducer {
    public:
      explicit L1Extra2PAT(const edm::ParameterSet & iConfig);
      virtual ~L1Extra2PAT() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::InputTag src_;
      int           bx_;
      edm::InputTag srcIso_, srcNonIso_, srcJC_, srcJF_, srcJT_;

      template<typename T>
      void getAndCopy(const edm::Event & iEvent, const edm::InputTag src, pat::TriggerObjectStandAloneCollection &out) ;

      // we use overloading here since template specialization at class scope is not permitted
      void copyMoreVars(const reco::LeafCandidate &l1, pat::TriggerObjectStandAlone & patl1) {  }
      void copyMoreVars(const l1extra::L1MuonParticle &l1, pat::TriggerObjectStandAlone & patl1) {
          patl1.setStatus(l1.gmtMuonCand().quality()); patl1.setPdgId(l1.gmtMuonCand().charge()*13);
      }
      void copyMoreVars(const l1extra::L1EmParticle &l1, pat::TriggerObjectStandAlone & patl1) {
          patl1.setPdgId(22); // photon, as I don't have the charge at L1
      }
  } ; // class
      
} // namespace

pat::L1Extra2PAT::L1Extra2PAT(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    bx_(iConfig.getParameter<uint32_t>("nbx")/2),
    srcIso_(src_.label(), "Isolated", src_.process()),
    srcNonIso_(src_.label(), "NonIsolated", src_.process()),
    srcJC_(src_.label(), "Central", src_.process()),
    srcJF_(src_.label(), "Forward", src_.process()),
    srcJT_(src_.label(), "Tau", src_.process())
{
   produces<pat::TriggerObjectStandAloneCollection>();
}

void 
pat::L1Extra2PAT::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
    std::auto_ptr<pat::TriggerObjectStandAloneCollection> out(new pat::TriggerObjectStandAloneCollection());

    getAndCopy<l1extra::L1MuonParticle>(iEvent, src_, *out);
    getAndCopy<l1extra::L1EmParticle>(iEvent, srcIso_, *out);
    getAndCopy<l1extra::L1EmParticle>(iEvent, srcNonIso_, *out);
    getAndCopy<l1extra::L1JetParticle>(iEvent, srcJF_, *out);
    getAndCopy<l1extra::L1JetParticle>(iEvent, srcJF_, *out);
    getAndCopy<l1extra::L1JetParticle>(iEvent, srcJF_, *out);

    iEvent.put(out);
}

template<typename T>
void
pat::L1Extra2PAT::getAndCopy(const edm::Event & iEvent, const edm::InputTag src, pat::TriggerObjectStandAloneCollection &out) 
{
    edm::Handle<std::vector<T> > handle;
    iEvent.getByLabel(src, handle);
    for (typename std::vector<T>::const_iterator it = handle->begin(), ed = handle->end(); it != ed; ++it) {
        if (abs(it->bx()) > bx_) continue;
        pat::TriggerObjectStandAlone patl1(*it);
        patl1.setCollection(src);
        patl1.addPathName("L1");
        copyMoreVars(*it, patl1);
        out.push_back(patl1);
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(L1Extra2PAT);
