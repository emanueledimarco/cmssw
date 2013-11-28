#include "FWCore/Framework/interface/MakerMacros.h"

#include "WWAnalysis/AnalysisStep/interface/CreateEventHists.h"
#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"
typedef edm::AnalyzerWrapper<CreateEventHists> CreateEventHistsEDMWrapped;
DEFINE_FWK_MODULE(CreateEventHistsEDMWrapped);
