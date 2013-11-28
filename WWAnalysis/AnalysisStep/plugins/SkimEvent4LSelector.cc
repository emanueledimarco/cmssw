#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"

typedef SingleObjectSelector<
            std::vector<reco::SkimEvent4L>,
            StringCutObjectSelector<reco::SkimEvent4L,true>,
            std::vector<reco::SkimEvent4L>
        > SkimEvent4LSelector;

typedef SingleObjectSelector<
            std::vector<reco::SkimEvent2L>,
            StringCutObjectSelector<reco::SkimEvent2L,true>,
            std::vector<reco::SkimEvent2L>
        > SkimEvent2LSelector;

DEFINE_FWK_MODULE(SkimEvent4LSelector);
DEFINE_FWK_MODULE(SkimEvent2LSelector);
