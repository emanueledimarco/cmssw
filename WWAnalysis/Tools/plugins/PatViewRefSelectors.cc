#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef SingleObjectSelector<
          edm::View<pat::Muon>,
          StringCutObjectSelector<reco::Candidate>
       > PATMuonViewRefSelector;
typedef SingleObjectSelector<
          edm::View<pat::Electron>,
          StringCutObjectSelector<reco::Candidate>
       > PATElectronViewRefSelector;

DEFINE_FWK_MODULE(PATMuonViewRefSelector);
DEFINE_FWK_MODULE(PATElectronViewRefSelector);
