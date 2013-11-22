#ifndef MuonMatcher_h
#define MuonMatcher_h

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class MuonMatcher{

public:

  //! ctor 
  MuonMatcher(edm::View<reco::Muon> coll) { collection_ = coll; } ;
  //!dtor 
  virtual ~MuonMatcher();

  reco::MuonRef matchByTrack(reco::MuonRef theMuon);

private:
  edm::View<reco::Muon> collection_;
};

#endif

