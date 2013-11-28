#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoMuon/TrackerSeedGenerator/plugins/CollectionCombiner.h"

typedef CollectionCombiner< reco::SuperClusterCollection > SuperClusterCombiner;

DEFINE_FWK_MODULE(SuperClusterCombiner);

