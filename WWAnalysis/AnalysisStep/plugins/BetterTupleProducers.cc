#include "WWAnalysis/AnalysisStep/interface/BetterTupleProducer.icc"

#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef BetterTupleProducer<pat::Electron> BetterElectronTupleProducer;
typedef BetterTupleProducer<pat::Muon> BetterMuonTupleProducer;
typedef BetterTupleProducer<reco::SkimEvent> SkimTupleProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BetterElectronTupleProducer);
DEFINE_FWK_MODULE(BetterMuonTupleProducer);
DEFINE_FWK_MODULE(SkimTupleProducer);
