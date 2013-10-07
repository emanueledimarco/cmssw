#ifndef SORTERS_H
#define SORTERS_H

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"

namespace sorters {
    struct ZCompositeCandidateMassSorter {

        bool operator() (const reco::CompositeCandidate& i, const reco::CompositeCandidate& j) { return (fabs(i.mass()-91.188) < fabs(j.mass()-91.188));}

    };

    struct ZMassSorter {

        bool operator() (const reco::SkimEvent2L& i, const reco::SkimEvent2L& j) { return (fabs(i.mass()-91.188) < fabs(j.mass()-91.188));}

    };

    struct HiggsZLeptonPtSorter {

        bool operator() (const reco::SkimEvent4L& i, const reco::SkimEvent4L& j) { return i.lpt(1, 0) + i.lpt(1, 1) > j.lpt(1, 0) + j.lpt(1, 1);}

    };

}


#endif
