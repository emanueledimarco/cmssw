import FWCore.ParameterSet.Config as cms

fourLeptonBlinder = cms.EDFilter("FourLeptonBlinder",
   electrons = cms.InputTag("boostedElectrons"), 
   muons = cms.InputTag("boostedMuons"), 
   electronSelection = cms.string("pt > 5"),
   muonSelection = cms.string("pt > 3"),
   blindMC   = cms.bool(False),
   blindData = cms.bool(True),
   #firstBlindedDataRun = cms.uint32(190000),
   firstBlindedDataRun = cms.uint32(0),
)
