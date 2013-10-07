// -*- C++ -*-
//
// Package:    PatrickMETProducer
// Class:      PatrickMETProducer
// 
/**\class PatrickMETProducer PatrickMETProducer.cc WWAnalysis/PatrickMETProducer/src/PatrickMETProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu Apr 14 23:15:22 CEST 2011
// $Id: PatrickMETProducer.cc,v 1.2 2011/04/23 15:01:13 mangano Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

//
// class declaration
//

class PatrickMETProducer : public edm::EDProducer 
{
public:
  explicit PatrickMETProducer(const edm::ParameterSet&);
  ~PatrickMETProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
  bool isJetFromPU(const reco::VertexCollection& vertexs, reco::PFJetCollection::const_iterator jet) const;

  // ----------member data ---------------------------
  double ptJetCut;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PatrickMETProducer::PatrickMETProducer(const edm::ParameterSet& iConfig)
{
   produces<reco::PFCandidateCollection>();  
   ptJetCut = iConfig.getParameter<double>("ptJetCut");  
}


PatrickMETProducer::~PatrickMETProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PatrickMETProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  auto_ptr< reco::PFCandidateCollection >  pOutput( new reco::PFCandidateCollection ); 

  Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel( "particleFlow", pfCandidates);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("offlinePrimaryVertices",vertices);



  Handle<reco::PFJetCollection> pfJets;
  iEvent.getByLabel( "ak5PFJets", pfJets);

  set< const reco::PFCandidate* >  allPUconsts;
  for(reco::PFJetCollection::const_iterator jet=pfJets->begin(); jet!=pfJets->end(); ++jet){
    if(jet->pt()>ptJetCut && 
       ((fabs(jet->eta())>2.4 &&	       
	 jet->nConstituents()>1 &&
	 jet->neutralEmEnergyFraction()<0.99 &&
	 jet->neutralHadronEnergyFraction()<0.99) ||
	(fabs(jet->eta())<=2.4 &&	       
	 jet->nConstituents()>1 &&
	 jet->neutralEmEnergyFraction()<0.99 &&
	 jet->neutralHadronEnergyFraction()<0.99 &&
	 jet->chargedEmEnergyFraction()<0.99 &&
	 jet->chargedHadronEnergyFraction()<0.99 &&
	 jet->chargedMultiplicity()>1)
	)){
      //cout << "I'm going to add the components of this jet to the output collection" << endl;
      vector< reco::PFCandidatePtr >  jetConsts = jet->getPFConstituents();
      //cout << "initial allPUconsts size: " << allPUconsts.size() << endl;
      //cout << "adding: " << jetConsts.size() << endl;
      for(vector< reco::PFCandidatePtr >::const_iterator constituent=jetConsts.begin(); 
	  constituent!=jetConsts.end();++constituent)
	allPUconsts.insert(constituent->get());
      
      //cout << "new allPUconsts size: " << allPUconsts.size() << endl;
    }
  }

  set< const reco::PFCandidate* >  allconsts;
  for(reco::PFCandidateCollection::const_iterator cand = pfCandidates->begin(); cand!=pfCandidates->end(); ++cand)
    allconsts.insert( &(*cand));

  //cout << "allPUconsts size: " << allPUconsts.size() << endl;
  //cout << "allconsts size: " << allconsts.size() << endl;
  vector<const reco::PFCandidate*> result;
  set_difference(allconsts.begin(),allconsts.end(),allPUconsts.begin(),allPUconsts.end(),
		 std::inserter(result, result.end()));
  //cout << "result size: " << result.size() <<endl;

  for(vector<const reco::PFCandidate*>::const_iterator it=result.begin(); it!=result.end(); ++it)
    pOutput->push_back( (**it));

  //cout << "pOutput size: " << pOutput->size() << endl;
  iEvent.put(pOutput);

 
}


bool PatrickMETProducer::isJetFromPU(const reco::VertexCollection& vertexs, 
				     reco::PFJetCollection::const_iterator jet) const {
  using namespace std;
  if(jet->chargedMultiplicity()==0) return false; //cannot associate to any PV
  reco::TrackRefVector  jetTracks = jet->getTrackRefs();

  double jetZ(0);
  int goodTrackCounter(0);
  for(reco::TrackRefVector::const_iterator track=jetTracks.begin(); track!=jetTracks.end(); ++track){
    if( (*track)->hitPattern().hasValidHitInFirstPixelBarrel() ){
      goodTrackCounter++;
      jetZ += (*track)->vz();
    }
  }
  // the averageZ is a temporary/quick solution. We should do some sort of vertexing of the jet's track
  jetZ = jetZ/goodTrackCounter;

  double dist(9999.);
  int selectedIndex(-1);
  int index(0);
  for(reco::VertexCollection::const_iterator vtx = vertexs.begin(); vtx !=vertexs.end(); ++vtx){
    float tmpDist = fabs(jetZ-vtx->position().z());
    if(tmpDist<dist){
      dist=tmpDist;
      selectedIndex=index;
    }
    index++;
  }
  //cout << "selectedIndex: " << selectedIndex << endl;
  if(selectedIndex!=0 && dist<0.5) {
    //cout << "dist with PU vtx: " << dist << endl;
    return true;
  }else  {
    //cout << "dist with Sig vtx: " << dist << endl;
    return false;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
PatrickMETProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PatrickMETProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PatrickMETProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PatrickMETProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PatrickMETProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PatrickMETProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PatrickMETProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatrickMETProducer);
