import FWCore.ParameterSet.Config as cms
from WWAnalysis.Filters.genFilterWW2L2Nu_cfi import genFilterWW2L2Nu as temp


def addPreYieldFilter(process,isMC,addBorisFilter,addVV=False):
    process.preYieldFilter = cms.Sequence()

    process.totalKinematicsFilter = cms.EDFilter('TotalKinematicsFilter',
        src       = cms.InputTag("genParticles"),
        tolerance = cms.double(0.5),
        verbose   = cms.untracked.bool(False)                                   
    )

    if isMC: process.preYieldFilter += process.totalKinematicsFilter

    process.genFilterWW2L2Nu = temp.clone()
    if addVV: process.preYieldFilter += process.genFilterWW2L2Nu

    process.nonSTAMuons = cms.EDFilter("MuonRefSelector",
        cut = cms.string("type!=8"),
        src = cms.InputTag("muons"),
        filter = cms.bool(False)
    )
    
    process.cleanRecoTaus = cms.EDFilter("PFTauSelector",
        src = cms.InputTag("hpsPFTauProducer"),
        discriminators = cms.VPSet(
            cms.PSet( 
              discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
              selectionCut=cms.double(0.5)
            )
        )
    )
    
    process.allLeps = cms.EDProducer("CandViewMerger",
        src = cms.VInputTag(
            cms.InputTag("gsfElectrons"), 
            cms.InputTag("nonSTAMuons"), 
            cms.InputTag("cleanRecoTaus")
        )
    )
    
    process.noTauLeps = cms.EDProducer("CandViewMerger",
        src = cms.VInputTag(
            cms.InputTag("gsfElectrons"), 
            cms.InputTag("nonSTAMuons"), 
        )
    )
    
    process.allDiLep = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string('allLeps noTauLeps'),
        cut = cms.string(
            'deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.05 && ' + 
            'min(daughter(0).pt,daughter(1).pt) >  8 && ' +
            'max(daughter(0).pt,daughter(1).pt) > 18'
        ),
        checkCharge = cms.bool(False)
    )
    
    process.countDiLeps  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("allDiLep"), minNumber = cms.uint32(1))
    process.preYieldFilter += ( process.nonSTAMuons * (process.cleanRecoTaus * process.allLeps + process.noTauLeps) * process.allDiLep * process.countDiLeps )
    
    process.genLepFromW10 = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        filter = cms.bool(True),
        cut = cms.string("(abs(pdgId)==11 || abs(pdgId)==13) && abs(mother.mother.pdgId)==24 &&" +
                         "abs(eta)<2.5 && pt>10.")
    )
    
    process.genLep10CountFilter=cms.EDFilter("CandViewCountFilter", 
                      src =cms.InputTag("genLepFromW10"), 
                      minNumber = cms.uint32(2)
    )
    
    process.genLepFromW20 = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genLepFromW10"),
        filter = cms.bool(True),
        cut = cms.string("pt>20.")
    )
    
    process.genFilter = cms.Sequence(process.genLepFromW10*process.genLep10CountFilter*process.genLepFromW20)
    
    if addBorisFilter:
        process.preYieldFilter.replace(
            process.nonSTAMuons,
            process.nonSTAMuons*
            process.genFilter
        )

