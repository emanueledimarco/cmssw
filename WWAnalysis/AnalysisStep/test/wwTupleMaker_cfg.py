import FWCore.ParameterSet.Config as cms

process = cms.Process("EWKTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
#     fileNames = cms.untracked.vstring('file:data/REPLACEME.root'),
#     fileNames = cms.untracked.vstring('file:test.root'),
#     fileNames = cms.untracked.vstring('file:wwEvents.root'),
    fileNames = cms.untracked.vstring('file:wEvents.root'),
)
process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_conditionsInEdm_*_*",  "drop *_MEtoEDMConverter_*_*")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.wweleltup = cms.EDAnalyzer("SkimTupleProducer", 
    src = cms.untracked.InputTag("wwelel"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("nlep"), quantity = cms.untracked.string("nLep()")),
#         cms.PSet( tag = cms.untracked.string("neel"), quantity = cms.untracked.string("nExtraEl()")),
#         cms.PSet( tag = cms.untracked.string("nel"), quantity = cms.untracked.string("nEl()")),
#         cms.PSet( tag = cms.untracked.string("pdgid1"), quantity = cms.untracked.string("pdgId(0)")),
#         cms.PSet( tag = cms.untracked.string("pdgid2"), quantity = cms.untracked.string("pdgId(1)")),
#         cms.PSet( tag = cms.untracked.string("pt1"), quantity = cms.untracked.string("pt(0)")),
#         cms.PSet( tag = cms.untracked.string("pt2"), quantity = cms.untracked.string("pt(1)")),
#         cms.PSet( tag = cms.untracked.string("eta1"), quantity = cms.untracked.string("eta(0)")),
#         cms.PSet( tag = cms.untracked.string("eta2"), quantity = cms.untracked.string("eta(1)")),
#         cms.PSet( tag = cms.untracked.string("phi1"), quantity = cms.untracked.string("phi(0)")),
#         cms.PSet( tag = cms.untracked.string("phi2"), quantity = cms.untracked.string("phi(1)")),
#         cms.PSet( tag = cms.untracked.string("q1"), quantity = cms.untracked.string("q(0)")),
#         cms.PSet( tag = cms.untracked.string("q2"), quantity = cms.untracked.string("q(1)")),
        cms.PSet( tag = cms.untracked.string("nj20"), quantity = cms.untracked.string("nJets(20)")),
        cms.PSet( tag = cms.untracked.string("ncj20"), quantity = cms.untracked.string("nCentralJets(20)")),
#         cms.PSet( tag = cms.untracked.string("nj30"), quantity = cms.untracked.string("nJets(30)")),
#         cms.PSet( tag = cms.untracked.string("ncj30"), quantity = cms.untracked.string("nCentralJets(30)")),
        cms.PSet( tag = cms.untracked.string("pf"), quantity = cms.untracked.string("pfMet()")),
        cms.PSet( tag = cms.untracked.string("tc"), quantity = cms.untracked.string("tcMet()")),
#         cms.PSet( tag = cms.untracked.string("dphilpf1"), quantity = cms.untracked.string("dPhilPfMet(0)")),
#         cms.PSet( tag = cms.untracked.string("dphiltc1"), quantity = cms.untracked.string("dPhilTcMet(0)")),
#         cms.PSet( tag = cms.untracked.string("dphilpf2"), quantity = cms.untracked.string("dPhilPfMet(1)")),
#         cms.PSet( tag = cms.untracked.string("dphiltc2"), quantity = cms.untracked.string("dPhilTcMet(1)")),
        cms.PSet( tag = cms.untracked.string("ppf"), quantity = cms.untracked.string("projPfMet()")),
        cms.PSet( tag = cms.untracked.string("ptc"), quantity = cms.untracked.string("projTcMet()")),
#         cms.PSet( tag = cms.untracked.string("ppf1"), quantity = cms.untracked.string("projPfMet(0)")),
#         cms.PSet( tag = cms.untracked.string("ptc1"), quantity = cms.untracked.string("projTcMet(0)")),
#         cms.PSet( tag = cms.untracked.string("ppf2"), quantity = cms.untracked.string("projPfMet(1)")),
#         cms.PSet( tag = cms.untracked.string("ptc2"), quantity = cms.untracked.string("projTcMet(1)")),
        cms.PSet( tag = cms.untracked.string("mll"), quantity = cms.untracked.string("mll()")),
        cms.PSet( tag = cms.untracked.string("dphill"), quantity = cms.untracked.string("dPhill()")),
        cms.PSet( tag = cms.untracked.string("ptll"), quantity = cms.untracked.string("pTll()")),
        cms.PSet( tag = cms.untracked.string("drll"), quantity = cms.untracked.string("dRll()")),
#         cms.PSet( tag = cms.untracked.string("pfmt1"), quantity = cms.untracked.string("pfMT(0)")),
#         cms.PSet( tag = cms.untracked.string("pfmt2"), quantity = cms.untracked.string("pfMT(1)")),
#         cms.PSet( tag = cms.untracked.string("tcmt1"), quantity = cms.untracked.string("tcMT(0)")),
#         cms.PSet( tag = cms.untracked.string("tcmt2"), quantity = cms.untracked.string("tcMT(1)")),
#         cms.PSet( tag = cms.untracked.string("tcmll"), quantity = cms.untracked.string("tcMTll()")),
#         cms.PSet( tag = cms.untracked.string("pfmll"), quantity = cms.untracked.string("pfMTll()")),
        #iso vars

        cms.PSet( tag = cms.untracked.string("tkLow"),        quantity = cms.untracked.string("tkIsoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("tkHigh"),       quantity = cms.untracked.string("tkIsoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("ecalLow"),      quantity = cms.untracked.string("ecalIsoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("ecalHigh"),     quantity = cms.untracked.string("ecalIsoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("hcalLow"),      quantity = cms.untracked.string("hcalIsoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("hcalHigh"),     quantity = cms.untracked.string("hcalIsoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("isoLow"),       quantity = cms.untracked.string("allIsoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("isoHigh"),      quantity = cms.untracked.string("allIsoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("tkVetoLow"),    quantity = cms.untracked.string("tkVetoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("tkVetoHigh"),   quantity = cms.untracked.string("tkVetoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("ecalVetoLow"),  quantity = cms.untracked.string("ecalVetoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("ecalVetoHigh"), quantity = cms.untracked.string("ecalVetoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("hcalVetoLow"),  quantity = cms.untracked.string("hcalVetoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("hcalVetoHigh"), quantity = cms.untracked.string("hcalVetoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("isoVetoLow"),   quantity = cms.untracked.string("allVetoByIso(0)")),
        cms.PSet( tag = cms.untracked.string("isoVetoHigh"),  quantity = cms.untracked.string("allVetoByIso(1)")),
        cms.PSet( tag = cms.untracked.string("pdgidLow"),     quantity = cms.untracked.string("pdgIdByIso(0)")),
        cms.PSet( tag = cms.untracked.string("pdgidHigh"),    quantity = cms.untracked.string("pdgIdByIso(1)")),
        cms.PSet( tag = cms.untracked.string("ptLow"),        quantity = cms.untracked.string("ptByIso(0)")),
        cms.PSet( tag = cms.untracked.string("ptHigh"),       quantity = cms.untracked.string("ptByIso(1)")),
        cms.PSet( tag = cms.untracked.string("etaLow"),       quantity = cms.untracked.string("etaByIso(0)")),
        cms.PSet( tag = cms.untracked.string("etaHigh"),      quantity = cms.untracked.string("etaByIso(1)")),
        cms.PSet( tag = cms.untracked.string("phiLow"),       quantity = cms.untracked.string("phiByIso(0)")),
        cms.PSet( tag = cms.untracked.string("phiHigh"),      quantity = cms.untracked.string("phiByIso(1)")),
        cms.PSet( tag = cms.untracked.string("qLow"),         quantity = cms.untracked.string("qByIso(0)")),
        cms.PSet( tag = cms.untracked.string("qHigh"),        quantity = cms.untracked.string("qByIso(1)")),
        cms.PSet( tag = cms.untracked.string("isEBLow"),      quantity = cms.untracked.string("isEBByIso(0)")),
        cms.PSet( tag = cms.untracked.string("isEBHigh"),     quantity = cms.untracked.string("isEBByIso(1)")),
        cms.PSet( tag = cms.untracked.string("tkPtLow"),      quantity = cms.untracked.string("tkPtByIso(0)")),
        cms.PSet( tag = cms.untracked.string("tkPtHigh"),     quantity = cms.untracked.string("tkPtByIso(1)")),

        cms.PSet( tag = cms.untracked.string("tk1"),          quantity = cms.untracked.string("tkIsoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("tk2"),          quantity = cms.untracked.string("tkIsoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("ecal1"),        quantity = cms.untracked.string("ecalIsoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("ecal2"),        quantity = cms.untracked.string("ecalIsoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("hcal1"),        quantity = cms.untracked.string("hcalIsoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("hcal2"),        quantity = cms.untracked.string("hcalIsoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("iso1"),         quantity = cms.untracked.string("allIsoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("iso2"),         quantity = cms.untracked.string("allIsoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("tkVeto1"),      quantity = cms.untracked.string("tkVetoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("tkVeto2"),      quantity = cms.untracked.string("tkVetoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("ecalVeto1"),    quantity = cms.untracked.string("ecalVetoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("ecalVeto2"),    quantity = cms.untracked.string("ecalVetoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("hcalVeto1"),    quantity = cms.untracked.string("hcalVetoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("hcalVeto2"),    quantity = cms.untracked.string("hcalVetoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("isoVeto1"),     quantity = cms.untracked.string("allVetoByPt(0)")),
        cms.PSet( tag = cms.untracked.string("isoVeto2"),     quantity = cms.untracked.string("allVetoByPt(1)")),
        cms.PSet( tag = cms.untracked.string("pdgid1"),       quantity = cms.untracked.string("pdgIdByPt(0)")),
        cms.PSet( tag = cms.untracked.string("pdgid2"),       quantity = cms.untracked.string("pdgIdByPt(1)")),
        cms.PSet( tag = cms.untracked.string("pt1"),          quantity = cms.untracked.string("ptByPt(0)")),
        cms.PSet( tag = cms.untracked.string("pt2"),          quantity = cms.untracked.string("ptByPt(1)")),
        cms.PSet( tag = cms.untracked.string("eta1"),         quantity = cms.untracked.string("etaByPt(0)")),
        cms.PSet( tag = cms.untracked.string("eta2"),         quantity = cms.untracked.string("etaByPt(1)")),
        cms.PSet( tag = cms.untracked.string("phi1"),         quantity = cms.untracked.string("phiByPt(0)")),
        cms.PSet( tag = cms.untracked.string("phi2"),         quantity = cms.untracked.string("phiByPt(1)")),
        cms.PSet( tag = cms.untracked.string("q1"),           quantity = cms.untracked.string("qByPt(0)")),
        cms.PSet( tag = cms.untracked.string("q2"),           quantity = cms.untracked.string("qByPt(1)")),
        cms.PSet( tag = cms.untracked.string("isEB1"),        quantity = cms.untracked.string("isEBByPt(0)")),
        cms.PSet( tag = cms.untracked.string("isEB2"),        quantity = cms.untracked.string("isEBByPt(1)")),
        cms.PSet( tag = cms.untracked.string("tkPt1"),        quantity = cms.untracked.string("tkPtByPt(0)")),
        cms.PSet( tag = cms.untracked.string("tkPt2"),        quantity = cms.untracked.string("tkPtByPt(1)")),
        cms.PSet( tag = cms.untracked.string("tkVeto"),       quantity = cms.untracked.string("tkVeto()")),
        cms.PSet( tag = cms.untracked.string("tkPt"),         quantity = cms.untracked.string("tkPt()")),
        cms.PSet( tag = cms.untracked.string("pt"),           quantity = cms.untracked.string("pt()")),
    ),
#all optional:
#     weight = cms.untracked.double(RMMEWEIGHT),
#     weight = cms.untracked.double(1),
#     vtxWeights = cms.untracked.PSet(),
#     vtxLabel = cms.InputTag('goodPrimaryVertices'),
)
process.wwelmutup = process.wweleltup.clone( src = 'wwelmu' )
process.wwmumutup = process.wweleltup.clone( src = 'wwmumu' )

# process.TFileService = cms.Service("TFileService",fileName = cms.string("data/REPLACEME.tuple.root"))
# process.TFileService = cms.Service("TFileService",fileName = cms.string("test.tuple.root"))
# process.TFileService = cms.Service("TFileService",fileName = cms.string("wwEvents.tuple.root"))
process.TFileService = cms.Service("TFileService",fileName = cms.string("wEvents.tuple.root"))

process.p = cms.Path(process.wweleltup + process.wwelmutup + process.wwmumutup)

