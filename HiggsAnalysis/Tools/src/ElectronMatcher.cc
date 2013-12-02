#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "HiggsAnalysis/Tools/interface/ElectronMatcher.hh"

ElectronMatcher::~ElectronMatcher() { }

reco::GsfElectronRef ElectronMatcher::matchByGsfTrack(reco::GsfElectronRef theEle) {

  reco::GsfTrackRef matchTrkRef = theEle->get<reco::GsfTrackRef>();

  for(int i=0; i<(int)collection_.size(); ++i) {
    reco::GsfElectronRef ele = collection_.refAt(i).castTo<reco::GsfElectronRef>();
    reco::GsfTrackRef trkRef = ele->get<reco::GsfTrackRef>();
    if(trkRef==matchTrkRef) {
      return ele;
    }
  }
  // if not found, return the reference to the electron itself
  return theEle;
}
