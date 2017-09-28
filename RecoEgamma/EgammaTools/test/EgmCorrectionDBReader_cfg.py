ONLY_GET=False

import FWCore.ParameterSet.Config as cms
process = cms.Process("elescaletxt")
process.load('Configuration.StandardSequences.Services_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

process.source = cms.Source("EmptySource")
process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = 'sqlite_file:Egm17_ele_prompt_Moriond17.db'

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDB,
    DumpStat=cms.untracked.bool(False),
    toGet = cms.VPSet(cms.PSet(
            record = cms.string('EgmCorrectorParametersRcd'),
            tag = cms.string('EgmCorrectorParameters_Ele_Prompt_Moriond17'),
            )),
)               

# this is just to test the retrieval of the payload in the record
process.get = cms.EDAnalyzer("EventSetupRecordDataGetter",
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('EgmCorrectorParametersRcd'),
        data = cms.vstring('EgmCorrectorParameters')
    )),
    verbose = cms.untracked.bool(True)
)


process.readEleScale    = cms.EDAnalyzer('EgmCorrectorDBReader',  
      # below is the communication to the database 
      record    = cms.untracked.string('EgmCorrectorParameters'),
      # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
      # but it is recommended to use the GT name that you retrieved the files from.
      globalTag      = cms.untracked.string('EgmTest'),  
      printScreen    = cms.untracked.bool(False),
      createTextFile = cms.untracked.bool(True)
)

if ONLY_GET:
    process.p = cms.Path(process.get)
else:
    process.p = cms.Path(process.readEleScale)
