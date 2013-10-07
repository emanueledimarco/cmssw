// -*- C++ -*-
//
// Package:    HtoWWElectrons
// Class:      hwwEleTrackerIsolation
// 
/*
   Description: <one line class summary>
   Electron isolation using tracker info

   Implementation:
 
*/
//
// Original Author:  Chiara Rovelli
//
//


#ifndef EleTrackerIsolationAlgo_h
#define EleTrackerIsolationAlgo_h

// my includes
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"

//CLHEP
#include <CLHEP/Vector/LorentzVector.h>

using namespace edm;
using namespace std;
using namespace reco;

class EleTrackerIsolationAlgo{
 public:
  
  //constructors
  EleTrackerIsolationAlgo();
  EleTrackerIsolationAlgo(const GsfElectron *gsfEle, const TrackCollection trackColl);

  //methods
  void setExtRadius (float extRadius);
  void setIntRadius (float intRadius);
  
  float getPtTracks () const;
  float minDeltaR (float minPt) const;
  float minDeltaR_withVeto (float minPt) const;
  bool isIsolated (float ptCut = 0.05) const;

  //destructor 
  ~EleTrackerIsolationAlgo();
  
 private:

  const GsfElectron*_myGsfEle;  	  
  const TrackCollection _tracks;
  
  float _extRadius;
  float _intRadius;
};

#endif
