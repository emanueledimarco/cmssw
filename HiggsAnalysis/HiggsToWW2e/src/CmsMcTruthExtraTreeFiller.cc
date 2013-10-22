//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsMcTruthExtraTreeFiller
//      Simple class for dumping MC truth info to an ntuple. 
//
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthExtraTreeFiller.h"


#include <string>

struct CmsMcTruthTreeFillerData {
  CmsTree *cmstree;
  
};

//----------------
// Constructors --
//----------------
CmsMcTruthExtraTreeFiller::CmsMcTruthExtraTreeFiller(CmsTree *cmstree):
  privateData_(new CmsMcTruthTreeFillerData)
{
  privateData_->cmstree=cmstree; 
  saveLHE_ = false;
}

//--------------
// Destructor --
//--------------
CmsMcTruthExtraTreeFiller::~CmsMcTruthExtraTreeFiller() {
  delete privateData_;
}

//-------------
// Methods   --
//-------------
void CmsMcTruthExtraTreeFiller::writeCollectionToTree(edm::InputTag mcTruthCollection, std::vector<std::string>* lheComments,
						       const edm::Event& iEvent, int range, bool firstEvent) {

  edm::Handle< reco::GenParticleCollection > genParticleHandle;
  try { iEvent.getByLabel(mcTruthCollection, genParticleHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("HWWTreeDumper") << "Can't get MC Truth Collection: " << mcTruthCollection; }
  const reco::GenParticleCollection *genParticleCollection = genParticleHandle.product();

  vector<float> pMC,thetaMC,etaMC,phiMC,energyMC;
  vector<int>   idMC,statusMC;

  vector<int> mothSt, grandmothSt;
  vector<int> mothId, grandmothId;

  int countParticles = 0;

  reco::GenParticleCollection::const_iterator genPart;
  for(genPart=genParticleCollection->begin(); genPart!=genParticleCollection->end(); genPart++) {
    
    const reco::Candidate & cand = *genPart;

    countParticles++;
    
    if( countParticles <= (range+1) ) continue; // we are only interested in the particles that were not stored in CmsMcTruthFiller class. That one has (range+1) particles 

    if( cand.status() == 1 && ( fabs(cand.pdgId()) == 11 || fabs(cand.pdgId()) == 13  ) ) { // store the information only for status 1 muons (|id| == 13) or electrons (|id| == 11)

      // Fill candidate information
      pMC.push_back(cand.p());
      thetaMC.push_back(cand.theta());
      etaMC.push_back(cand.eta());
      phiMC.push_back(cand.phi());
      energyMC.push_back(cand.energy());
      idMC.push_back(cand.pdgId());
      statusMC.push_back(cand.status());

      const reco::Candidate *mom = cand.mother(); // take the mother
      
      // Avoid self-mom (grandmom) particles
      // E.g. : if we have the case
      //  ____________________
      // |                    |
      // | id  : 13           |
      // | st  : 1            |
      // | mother st  : 2     |
      // | mother id  : 15    |
      // | gmother st  : 3    |
      // | gmother id  : 15   |
      // |____________________|
      //
      // it will become
      //  ____________________
      // |                    |
      // | id  : 13           |
      // | st  : 1            |
      // | mother st  : 2     |
      // | mother id  : 15    |
      // | gmother st  : 3    |
      // | gmother id  : -24  |
      // |____________________|
      //
      // where we took the mother of the tau particle

      // First check cand-mother
      int canId = cand.pdgId();
      int momId = mom->pdgId();

      if (canId == momId){
	
	int tempId = momId;
	while (tempId == canId){
	  mom = mom->mother();
	  tempId = mom->pdgId();
	}

      }// end check cand-mother

      // Fill mother information
      mothSt.push_back(mom->status());
      mothId.push_back(mom->pdgId());

      // Now check mother-grandmother
      mom = mom->mother(); // take the grandmother
      int gmomId = mom->pdgId();
      
      if(momId == gmomId){

	int tempId = gmomId;
	while(tempId == momId){
	  mom = mom->mother();
	  tempId = mom->pdgId();
	}

      } // end check mother-grandmother

      // Fill grandmother information
      grandmothSt.push_back(mom->status());
      grandmothId.push_back(mom->pdgId());

    }

  }// for genParticles


  // LHE Event
  // to be written in the ntuple
  if(false) {
    lheComments->clear();
    edm::Handle<LHEEventProduct> product;
    iEvent.getByLabel("source", product);
    std::vector<std::string>::const_iterator c_begin = product->comments_begin();
    std::vector<std::string>::const_iterator c_end = product->comments_end();
    for( std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
      lheComments->push_back(*cit);
    }
  }
  
  std::string theName = "tHMc"; // name changed wrt standard collection of gen particles -- was "Mc"
  std::string indName = "n"+theName;
  
  // common variables as for standard gen particle collection
  privateData_->cmstree->column(indName.c_str(),(int)pMC.size(),0,theName.c_str());
  privateData_->cmstree->column(("p"+theName).c_str(),      pMC,      indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("theta"+theName).c_str(),  thetaMC,  indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("eta"+theName).c_str(),    etaMC,    indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("phi"+theName).c_str(),    phiMC,    indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("energy"+theName).c_str(), energyMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("id"+theName).c_str(),     idMC,     indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("status"+theName).c_str(), statusMC, indName.c_str(), 0, theName.c_str());

  // new variables for tH analysis 
  privateData_->cmstree->column(("mothSt"+theName).c_str(),      mothSt,      indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("grandmothSt"+theName).c_str(), grandmothSt, indName.c_str(), 0, theName.c_str());

  privateData_->cmstree->column(("mothId"+theName).c_str(),      mothId,      indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("grandmothId"+theName).c_str(), grandmothId, indName.c_str(), 0, theName.c_str());

  if(firstEvent && saveLHE_) privateData_->cmstree->getTree()->Branch("commentLHE", &(*lheComments));

}
  

  
