import FWCore.ParameterSet.Config as cms

process = cms.Process("skimEventCombination")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_1_1_Us9.root",
#"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_2_1_l6b.root",
#"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_3_1_vMP.root",
#"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_4_1_J8Y.root",
"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_5_1_WTG.root",
"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_6_1_dEd.root",
"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_7_1_SGZ.root",
"file:/nfs/bluearc/group/trees/hww/test/101160.mergeGluGlu160Vtx//ggToH160toWWto2L2Nu_8_1_gMd.root",
    )
)

### ID
process.crossID = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwelmuIDMerge"),
  inputB = cms.InputTag("wwmuelIDMerge"),
)

process.sfID = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwmumuIDMerge"),
  inputB = cms.InputTag("wwelelIDMerge"),
)

process.allID = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("crossID"),
  inputB = cms.InputTag("sfID"),
)


### ISO
process.crossISO = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwelmuISOMerge"),
  inputB = cms.InputTag("wwmuelISOMerge"),
)

process.sfISO = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwmumuISOMerge"),
  inputB = cms.InputTag("wwelelISOMerge"),
)

process.allISO = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("crossISO"),
  inputB = cms.InputTag("sfISO"),
)


### CONV
process.crossCONV = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwelmuCONVMerge"),
  inputB = cms.InputTag("wwmuelCONVMerge"),
)

process.sfCONV = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwmumuCONVMerge"),
  inputB = cms.InputTag("wwelelCONVMerge"),
)

process.allCONV = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("crossCONV"),
  inputB = cms.InputTag("sfCONV"),
)


### IP
process.crossIP = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwelmuIPMerge"),
  inputB = cms.InputTag("wwmuelIPMerge"),
)

process.sfIP = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("wwmumuIPMerge"),
  inputB = cms.InputTag("wwelelIPMerge"),
)

process.allIP = cms.EDProducer('SkimEventCombiner',
  inputA = cms.InputTag("crossIP"),
  inputB = cms.InputTag("sfIP"),
)


### OUTPUT
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('mergedSkimEvents2.root'),
    outputCommands = cms.untracked.vstring(
        'keep *_*_*_*',
        'drop *_sf*_*_*',
        )
)
  


process.p = cms.Path(process.crossID*process.sfID*process.allID*
                     process.crossISO*process.sfISO*process.allISO*
                     process.crossCONV*process.sfCONV*process.allCONV*
                     process.crossIP*process.sfIP*process.allIP
)

process.e = cms.EndPath(process.out)
