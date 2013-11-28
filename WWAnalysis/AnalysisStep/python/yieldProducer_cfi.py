import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.cutPSets_cfi import *

#default yield params setup for wwelel
FWLiteParams = cms.PSet(

    sampleName  = cms.string('test'),
#     fileOutName = cms.string('output.root'),
    puWeights = cms.vdouble(),
# uncomment if needed
#     puLabel = cms.InputTag("addPileupInfo"),
#     ptWeight = cms.InputTag("higgsPt")
    doNMinus1 = cms.bool(True),
    doByCuts  = cms.bool(True),
    printSummary  = cms.bool(False),

    hypotheses = cms.PSet( 
        wwelel0 = cms.PSet(
            src = cms.InputTag("wwelel0"),
            cuts = cloneVPSet(defaultWW)
        ),
    ),

    histParams = cms.PSet(
        mll = cms.PSet(
            variable = cms.string("mll"),
            nbins = cms.uint32(300),
            low = cms.double(0),
            high = cms.double(300),
            xtitle = cms.string("m_{ll} [GeV]")
        ),
        ptMax = cms.PSet(
            variable = cms.string("ptMax"),
            nbins = cms.uint32(200),
            low = cms.double(0),
            high = cms.double(200),
            xtitle = cms.string("p_{T} [GeV]")
        ),
        ptMin = cms.PSet(
            variable = cms.string("ptMin"),
            nbins = cms.uint32(200),
            low = cms.double(0),
            high = cms.double(200),
            xtitle = cms.string("p_{T} [GeV]")
        ),
        dPhill = cms.PSet(
            variable = cms.string("dPhillInDegrees"),
            nbins = cms.uint32(180),
            low = cms.double(0),
            high = cms.double(180),
            xtitle = cms.string("#Delta#phi_{ll} [deg]")
        ),
        mTHiggs = cms.PSet(
            variable = cms.string("mTHiggs"),
            nbins = cms.uint32(300),
            low = cms.double(0),
            high = cms.double(300),
            xtitle = cms.string("m_{T} [GeV]")
        ),
        nJets = cms.PSet(
            variable = cms.string("nCentralJets(30.0,5.0,1)"),
            nbins = cms.uint32(10),
            low = cms.double(-0.5),
            high = cms.double(9.5),
            xtitle = cms.string("N_{Jets}")
        ),
#         tcMET = cms.PSet(
#             variable = cms.string("tcMet"),
#             nbins = cms.uint32(200),
#             low = cms.double(0.0),
#             high = cms.double(200),
#             xtitle = cms.string("tc#slash{E}_{T} [GeV]")
#         ),
        pfMET = cms.PSet(
            variable = cms.string("pfMet"),
            nbins = cms.uint32(200),
            low = cms.double(0.0),
            high = cms.double(200),
            xtitle = cms.string("pf#slash{E}_{T} [GeV]")
        ),
#         ptcMET = cms.PSet(
#             variable = cms.string("projTcMet"),
#             nbins = cms.uint32(200),
#             low = cms.double(0.0),
#             high = cms.double(200),
#             xtitle = cms.string("tc#slash{E}_{T} [GeV]")
#         ),
        ppfMET = cms.PSet(
            variable = cms.string("projPfMet"),
            nbins = cms.uint32(200),
            low = cms.double(0.0),
            high = cms.double(200),
            xtitle = cms.string("pf#slash{E}_{T} [GeV]")
        ),
        nVtx = cms.PSet(
            variable = cms.string("min(nGoodVertices,20.5)"),
            nbins = cms.uint32(20),
            low = cms.double(0.5),
            high = cms.double(20.5),
            xtitle = cms.string("N_{vtx}^{good}")
        ),
    )
)





