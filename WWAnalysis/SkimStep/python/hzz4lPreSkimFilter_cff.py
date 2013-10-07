import FWCore.ParameterSet.Config as cms

hzz4lPreFilterMuons = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("(isGlobalMuon || numberOfMatches > 0) && track.isNonnull && pt > 3"),
)

hzz4lPreFilterElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string("pt > 0"),
)
hzz4lPreFilterLeps = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(
        cms.InputTag("hzz4lPreFilterMuons"), 
        cms.InputTag("hzz4lPreFilterElectrons"), 
   )
)


hzz4lPreFilterDiLep = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('hzz4lPreFilterLeps hzz4lPreFilterLeps'),
    cut = cms.string('deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.05'),
    checkCharge = cms.bool(False)
)

hzz4lPreFilterDiLep178Filter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hzz4lPreFilterDiLep"),
    cut = cms.string('min(daughter(0).pt,daughter(1).pt) >  8 && ' +
                     'max(daughter(0).pt,daughter(1).pt) > 17'),
    filter = cms.bool(True),
)

hzz4lPreFilterMll40Filter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("hzz4lPreFilterDiLep"),
    cut = cms.string('mass > 40 && daughter(0).pdgId = -daughter(1).pdgId'), ## and SF/OS
    filter = cms.bool(True),
)

hzz4lPreFilterTriLepFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("hzz4lPreFilterLeps"),
    minNumber = cms.uint32(3),
)

hzz4lPreFilterSeq  = cms.Sequence(
    (hzz4lPreFilterMuons + hzz4lPreFilterElectrons) *
    hzz4lPreFilterLeps *
    hzz4lPreFilterTriLepFilter
)

