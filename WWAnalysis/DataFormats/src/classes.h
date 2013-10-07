#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"
#include "DataFormats/Common/interface/RefVector.h"
#include <DataFormats/Common/interface/BaseHolder.h>
#include "CMGTools/External/interface/PileupJetIdentifier.h"

namespace
{
  struct dictionary {

    reco::SkimEvent dummy1;
    std::vector<reco::SkimEvent> dummy2;
    edm::Wrapper<reco::SkimEvent> dummy3;
    edm::Wrapper<std::vector<reco::SkimEvent> > dummy4;

    reco::SkimEvent4L dummy4L1;
    reco::SkimEvent4L::Pair dummy4L1p;
    std::vector<reco::SkimEvent4L> dummy4L2;
    edm::Wrapper<reco::SkimEvent4L> dummy4L3;
    edm::Wrapper<std::vector<reco::SkimEvent4L> > dummy4L4;

    reco::SkimEvent2L dummy2L1;
    std::vector<reco::SkimEvent2L> dummy2L2;
    edm::Wrapper<reco::SkimEvent2L> dummy2L3;
    edm::Wrapper<std::vector<reco::SkimEvent2L> > dummy2L4;


    edm::ValueMap<reco::PFMET> mwlDummy01;
    edm::Wrapper<edm::ValueMap<reco::PFMET> > mwlDummy02;
    edm::RefVector<std::vector<reco::Vertex> > mwlDummy03;
    edm::Wrapper<edm::RefVector<std::vector<reco::Vertex> > > mwlDummy04;
    //std::vector<edm::RefToBase<reco::RecoCandidate> > reftobaseRecoCand;
    edm::RefToBase<reco::RecoCandidate> refToBaseRecoCandidateDummy;
    edm::RefToBaseVector<reco::RecoCandidate> refToBaseVectorRecoCandidateDummy;
    edm::Ptr<reco::RecoCandidate> ptrRecoCandidateDummy;
    std::vector<edm::Ptr<reco::RecoCandidate> > vectorPtrRecoCandidateDummy;
    //edm::reftobase::BaseVectorHolder<reco::RecoCandidate> reftobaseBaseVectorHolderRecoCandidateDummy;

    //edm::reftobase::BaseHolder<reco::RecoCandidate> reftobaseBaseHolderRecoCandidateDummy;

    edm::Wrapper<edm::ValueMap<StoredPileupJetIdentifier> > DummyStoredPileupJetIdentifier ; 
    std::vector<StoredPileupJetIdentifier> vectorDummyStoredPileupJetIdentifier ;


//     SimpleJetCorrector dummy9;
//     FactorizedJetCorrector dummy10;
//     std::vector<SimpleJetCorrector*> dummy6;
//     std::vector<FactorizedJetCorrector::LevelTypes> dummy5;
//     std::vector<FactorizedJetCorrector::VarTypes> dummy8;
//     std::vector<std::vector<FactorizedJetCorrector::VarTypes> > dummy7;
// 
//     edm::Wrapper<SimpleJetCorrector > dummy91;
//     edm::Wrapper<FactorizedJetCorrector > dummy101;
//     edm::Wrapper<std::vector<SimpleJetCorrector*> > dummy61;
//     edm::Wrapper<std::vector<FactorizedJetCorrector::LevelTypes> > dummy51;
//     edm::Wrapper<std::vector<FactorizedJetCorrector::VarTypes> > dummy81;
//     edm::Wrapper<std::vector<std::vector<FactorizedJetCorrector::VarTypes> > > dummy71;
// 
//   <class name="SimpleJetCorrector"/>
//   <class name="FactorizedJetCorrector"/>
//   <class name="std::vector<SimpleJetCorrector*>"/>
//   <class name="std::vector<FactorizedJetCorrector::LevelTypes>"/>
//   <class name="std::vector<FactorizedJetCorrector::VarTypes>"/>
//   <class name="std::vector<std::vector<FactorizedJetCorrector::VarTypes> >"/>
//   <class name="edm::Wrapper<SimpleJetCorrector >"/>
//   <class name="edm::Wrapper<FactorizedJetCorrector >"/>
//   <class name="edm::Wrapper<std::vector<SimpleJetCorrector*> >"/>
//   <class name="edm::Wrapper<std::vector<FactorizedJetCorrector::LevelTypes> >"/>
//   <class name="edm::Wrapper<std::vector<FactorizedJetCorrector::VarTypes> >"/>
//   <class name="edm::Wrapper<std::vector<std::vector<FactorizedJetCorrector::VarTypes> > >"/>


  };
}

