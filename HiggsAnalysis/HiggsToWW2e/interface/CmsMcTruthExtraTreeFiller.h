// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HtoWWTreeDumper
// Description:
//      Class CmsMcTruthTreeFiller
//      Simple class for dumping RECO (or AOD) contents to an ntuple
//
//-----------------------------------------------------------------------

#ifndef CmsMcTruthExtraTreeFiller_h
#define CmsMcTruthExtraTreeFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include <string>
#include <TTree.h>


class CmsMcTruthTreeFillerData;

class CmsMcTruthExtraTreeFiller {

 public:

  // Constructor
  CmsMcTruthExtraTreeFiller(CmsTree *);

  // Destructor
  virtual ~CmsMcTruthExtraTreeFiller();

  // Write the content of the collection
  void writeCollectionToTree( edm::InputTag mcTruthCollection, std::vector<std::string>* lheComments, const edm::Event& iEvent, int range=100, bool firstEvent=false );

  // Modifiers
  void saveLHEComments(bool what) {saveLHE_ = what; }

 private:
 
  CmsMcTruthTreeFillerData *privateData_;
  bool saveLHE_;

};

#endif // CmsMcTruthExtraTreeFiller_h
