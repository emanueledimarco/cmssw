#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "TMath.h"

class ElectronPFIsoSingleTypeMapProd : public edm::EDProducer {
public:
  explicit ElectronPFIsoSingleTypeMapProd(const edm::ParameterSet&);
  ~ElectronPFIsoSingleTypeMapProd();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag vtxLabel_;
  edm::InputTag eleLabel_;
  edm::InputTag pfLabel_;
  std::vector<int> pfTypes_;
  double deltaR_;
  double neutralThreshold_;
  bool directional_;
  bool radial_;

};



ElectronPFIsoSingleTypeMapProd::ElectronPFIsoSingleTypeMapProd(const edm::ParameterSet& iConfig) :
  vtxLabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxLabel")),
  eleLabel_(iConfig.getUntrackedParameter<edm::InputTag>("eleLabel")),
  pfLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfLabel")),
  pfTypes_(iConfig.getUntrackedParameter<std::vector<int> >("pfTypes")),
  deltaR_(iConfig.getUntrackedParameter<double>("deltaR")),
  neutralThreshold_(iConfig.getUntrackedParameter<double>("neutralThreshold")),
  directional_(iConfig.getUntrackedParameter<bool>("directional")),
  radial_(iConfig.getUntrackedParameter<bool>("radial")) {

    produces<edm::ValueMap<float> >().setBranchAlias("pfTypeElIso");
}

void ElectronPFIsoSingleTypeMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::GsfElectronCollection> eleH;
  iEvent.getByLabel(eleLabel_,eleH);

  //  edm::Handle<reco::PFCandidateCollection> pfH;
  edm::Handle< std::vector< edm::FwdPtr<reco::PFCandidate> > > pfH;
  iEvent.getByLabel(pfLabel_,pfH);

  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByLabel(vtxLabel_,vtxH);

  std::vector<float> isoV;
  std::auto_ptr<edm::ValueMap<float> > isoM(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler isoF(*isoM);

  for(size_t i=0; i<eleH->size();++i) {
    const reco::GsfElectron &ele = eleH->at(i);

//     Double_t zLepton = 0.0;
//     if(ele.track().isNonnull()) zLepton = ele.track()->dz(vtxH->at(0).position());

    Double_t ptSum =0.;  
    Double_t ptSum_radial =0.;  
    math::XYZVector isoAngleSum;
    std::vector<math::XYZVector> coneParticles;

    for (size_t j=0; j<pfH->size();j++) {   
      const edm::FwdPtr<reco::PFCandidate> &pf = pfH->at(j);

      // consider only the types requested
      bool skip=true;
      for(size_t t=0; t<pfTypes_.size(); t++) {
        if(pf->particleId() == pfTypes_[t]) {
          skip=false;
          break;
        }
      }
      if(skip) continue;

      // exclude electron // this would be good, but for synch with smurfs removed
      //       if(pf->trackRef().isNonnull() && ele.closestCtfTrackRef().isNonnull() &&
      //          pf->trackRef() == ele.closestCtfTrackRef()) continue;

      // exclude electron
      //       if(pf->gsfTrackRef().isNonnull() && ele.gsfTrack().isNonnull() &&
      //          pf->gsfTrackRef() == ele.gsfTrack()) continue;


      // pt cut applied to neutrals
      if(!pf->trackRef().isNonnull() && pf->pt() < neutralThreshold_) continue;

      // ignore the pf candidate if it is too far away in Z
      //       Double_t deltaZ = 0.0;
      //       if(pf->trackRef().isNonnull()) {
      //         deltaZ = fabs(pf->trackRef()->dz(vtxH->at(0).position()) - zLepton);
      //       }
      //       if (deltaZ >= 0.1) continue;


      Double_t dr = ROOT::Math::VectorUtil::DeltaR(ele.momentum(), pf->momentum());
      // add the pf pt if it is inside the extRadius and outside the intRadius
      if ( dr < deltaR_ && dr >= 0.0 ) {

        // remove the eventual duplicates in pf-reco if the pfcand does not make a jpsi with ele
//         if (pf->particleId() == reco::PFCandidate::e) {
//           math::XYZTLorentzVector eleP4 = ele.p4();
//           math::XYZTLorentzVector pfP4 = pf->p4();
//           float invmass = (eleP4 + pfP4).mass();
//           if(invmass<2.5 || invmass>3.6) continue;
//         }

        // vetoes optimization from:
        // https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=154207

        // supercluster veto on electrons failing missing hits cut
        if (ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits()>0 && pf->mva_nothing_gamma() > 0.99 && 
            ele.superCluster().isNonnull() && pf->superClusterRef().isNonnull() && 
            ele.superCluster() == pf->superClusterRef()) continue;

        // dR Veto for Gamma: no-one in EB, dR > 0.08 in EE
        if (pf->particleId() == reco::PFCandidate::gamma && fabs(ele.superCluster()->eta())>1.479 
            && ROOT::Math::VectorUtil::DeltaR(ele.momentum(), pf->momentum()) < 0.08) continue;
        
        // charged hadron: no-one in EB, dR > 0.015 in EE
        if(pf->trackRef().isNonnull() && fabs(ele.superCluster()->eta())>1.479
           && ROOT::Math::VectorUtil::DeltaR(ele.momentum(), pf->momentum()) < 0.015) continue; 

        // neutral hadron: no-one in EB, InnerCone (One Tower = dR < 0.07) Veto for non-gamma neutrals
        // if (!pf->trackRef().isNonnull() && pf->particleId() == reco::PFCandidate::h0 && fabs(ele.superCluster()->eta())>1.479 
        //      && ROOT::Math::VectorUtil::DeltaR(ele.momentum(), pf->momentum()) < 0.07 ) continue; 
        
        // scalar sum
        ptSum += pf->pt();            
        
        // directional sum
        math::XYZVector transverse( pf->eta() - ele.eta()
                                    , reco::deltaPhi(pf->phi(), ele.phi())
                                    , 0);
        transverse *= pf->pt() / transverse.rho();
        if (transverse.rho() > 0) {
          isoAngleSum += transverse;
          coneParticles.push_back(transverse);
        }

        
        // for radial isolation:
        ptSum_radial += pf->pt() * ( 1. - 10.*dr*deltaR_ );


      } //if deltaR

    } // for pf candidates

    if (directional_) {
      double directionalPT = 0;
      for (unsigned int iPtcl = 0; iPtcl < coneParticles.size(); ++iPtcl)
        directionalPT += pow(TMath::ACos( coneParticles[iPtcl].Dot(isoAngleSum) / coneParticles[iPtcl].rho() / isoAngleSum.rho() ),2) * coneParticles[iPtcl].rho();
      isoV.push_back(directionalPT);
    } else if( radial_ ) {
      isoV.push_back( ptSum_radial );
    } else isoV.push_back(ptSum);


  } // for electrons

  isoF.insert(eleH,isoV.begin(),isoV.end());

  isoF.fill();
  iEvent.put(isoM);

}

ElectronPFIsoSingleTypeMapProd::~ElectronPFIsoSingleTypeMapProd() { }
void ElectronPFIsoSingleTypeMapProd::beginJob() { }
void ElectronPFIsoSingleTypeMapProd::endJob() { }
DEFINE_FWK_MODULE(ElectronPFIsoSingleTypeMapProd);
