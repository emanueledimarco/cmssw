#include "FWCore/Framework/interface/MakerMacros.h"
#include "WWAnalysis/AnalysisStep/interface/StringSortCollectionSelector.h"
#include "WWAnalysis/AnalysisStep/interface/SortSelector.h"
#include "WWAnalysis/AnalysisStep/interface/StringComparator.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"

typedef SortSelector<
          edm::View<reco::Candidate>,
  StringComparator<reco::Candidate, true>,
          reco::CandidateCollection
       > CandViewSorter;

typedef SortSelector<
          edm::View<reco::SkimEvent2L>,
  StringComparator<reco::Candidate, true>,
          std::vector<reco::SkimEvent2L>
       > SkimEvent2LSorter;

typedef SortSelector<
          edm::View<reco::SkimEvent4L>,
  StringComparator<reco::Candidate, true>,
          std::vector<reco::SkimEvent4L>
       > SkimEvent4LSorter;


DEFINE_FWK_MODULE(CandViewSorter);
DEFINE_FWK_MODULE(SkimEvent2LSorter);
DEFINE_FWK_MODULE(SkimEvent4LSorter);
