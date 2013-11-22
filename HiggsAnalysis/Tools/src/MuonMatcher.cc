#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "HiggsAnalysis/Tools/interface/MuonMatcher.hh"

MuonMatcher::~MuonMatcher() { }

reco::MuonRef MuonMatcher::matchByTrack(reco::MuonRef theMuon) {

  reco::TrackRef matchTrkRef = theMuon->track();

  for(int i=0; i<(int)collection_.size(); ++i) {
    reco::MuonRef mu = collection_.refAt(i).castTo<reco::MuonRef>();
    reco::TrackRef trkRef = mu->track();
    if(trkRef==matchTrkRef) {
      return mu;
    }
  }
  // if not found, return the reference to the electron itself
  return theMuon;
}
