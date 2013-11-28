#ifndef HWWUtils_h
#define HWWUtils_h
// -*- C++ -*-
//
// Package:    VBFProcessFilter
// Class:      VBFProcessFilter
// 
/* 
 
 Description: filter events based on the Pythia ProcessID and the Pt_hat
 Implementation: inherits from generic EDFilter
 
 */
//
// $Id: HWWUtils.h,v 1.1 2008/02/14 16:07:21 ceballos Exp $
//
//
// system include files
#include <memory>

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "TLorentzVector.h"

typedef reco::CaloJetCollection::const_iterator HWWjetIt ;

void 
setMomentum (TLorentzVector & myvector, 
             const reco::Candidate & gen) ;

std::pair<HWWjetIt,HWWjetIt>	
findTagJets (HWWjetIt begin, HWWjetIt end,
             double jetPtMin, double jetEtaMax) ;

double deltaPhi (double phi1, double phi2) ;

#endif
