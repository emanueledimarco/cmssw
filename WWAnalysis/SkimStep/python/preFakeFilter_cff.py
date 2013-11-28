import FWCore.ParameterSet.Config as cms

def addPreFakeFilter(process,isMC):
    ### Pre-filters for fake rate studies
    from HLTrigger.HLTfilters.triggerResultsFilter_cfi import triggerResultsFilter
    process.hltFilter4Fakes = triggerResultsFilter.clone(
        l1tResults = '',
        hltResults = cms.InputTag( "TriggerResults", "", "HLT"),
        throw = True,
        triggerConditions = [ 'HLT_Mu8_v*',
                              'HLT_Mu15_v*',
                              'HLT_Mu24_v*',
                              'HLT_Ele8_v*',
                              'HLT_Ele8_CaloIdL_CaloIsoVL_v*',
                              'HLT_Ele17_CaloIdL_CaloIsoVL_v*' ]
    )
    
    process.jetsPt15 = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("ak5PFJets"),
        cut = cms.string("pt > 15 && numberOfDaughters >= 3"),
        filter = cms.bool(True),
    )
    process.leptonPlusJet = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string('noTauLeps jetsPt15'),
        cut = cms.string(
            'deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.7 && ' +
            'daughter(0).pt >= 10' 
        ),
        checkCharge = cms.bool(False),
    )
    process.countLeptonPlusJet = process.countDiLeps.clone(src = "leptonPlusJet")
    
    process.diLeptons4Veto = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string('noTauLeps@+ noTauLeps@-'),
#        cut = cms.string('abs(daughter(0).pdgId) == abs(daughter(1).pdgId) && (mass < 12 || abs(mass-91.1876) < 15)')
        cut = cms.string('abs(daughter(0).pdgId) == abs(daughter(1).pdgId) && (mass < 12)')
    )
    process.diLeptons4VetoFilter = process.countDiLeps.clone(src = "diLeptons4Veto")
    
    process.metVeto20 = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("pfMet"),
        cut = cms.string("pt < 20"),
        filter = cms.bool(True),
    )
    process.recoW4Veto = cms.EDProducer("CandViewShallowCloneCombiner",
        decay = cms.string('noTauLeps pfMet'),
        cut = cms.string("sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi))) > 20"),
        checkCharge = cms.bool(False)
    )
    process.recoW4VetoFilter = process.countDiLeps.clone(src = "recoW4Veto")
    
    process.preFakeFilter = cms.Sequence()
    if not isMC:
        process.preFakeFilter += process.hltFilter4Fakes 

    process.preFakeFilter += (
        process.metVeto20 +
        process.nonSTAMuons * process.noTauLeps *
        ( process.jetsPt15 * process.leptonPlusJet * process.countLeptonPlusJet +
          process.diLeptons4Veto * ~process.diLeptons4VetoFilter +
          process.recoW4Veto     * ~process.recoW4VetoFilter )
    )

