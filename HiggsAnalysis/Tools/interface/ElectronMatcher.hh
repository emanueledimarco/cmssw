#ifndef ElectronMatcher_h
#define ElectronMatcher_h

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

class ElectronMatcher{

public:

  //! ctor 
  ElectronMatcher(edm::View<reco::GsfElectron> coll) { collection_ = coll; } ;
  //!dtor 
  virtual ~ElectronMatcher();

  reco::GsfElectronRef matchByGsfTrack(reco::GsfElectronRef theEle);

private:
  edm::View<reco::GsfElectron> collection_;
};

#endif

