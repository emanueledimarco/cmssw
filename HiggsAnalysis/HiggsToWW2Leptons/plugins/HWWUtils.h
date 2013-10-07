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
// $Id: VBFUtils.h,v 1.3 2008/02/05 13:08:53 tancini Exp $
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
