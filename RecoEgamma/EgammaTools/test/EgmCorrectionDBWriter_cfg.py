import FWCore.ParameterSet.Config as cms 
process = cms.Process('egmscaledb') 

process.source = cms.Source('EmptySource') 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1)) 

process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = 'sqlite_file:Egm17_ele_prompt_Moriond17.db'

process.PoolDBOutputService = cms.Service('PoolDBOutputService', 
   process.CondDB,
   timetype = cms.untracked.string('runnumber'),
   toPut = cms.VPSet( 
      cms.PSet(
         record = cms.string('EgmCorrectorParametersRcd'),
         tag    = cms.string('EgmCorrectorParameters_Ele_Prompt_Moriond17'), 
      ),
   ) 
) 

process.dbWriterEleScale = cms.EDAnalyzer('EgmCorrectorDBWriter', 
   object   = cms.untracked.string("electrons"),
   name   = cms.untracked.string("scale"),
   reco   = cms.untracked.string('promptreco'), 
   version   = cms.untracked.string('Moriond17'),
   path   = cms.untracked.string('RecoEgamma/EgammaTools/data/'),
   record = cms.untracked.string('EgmCorrectorParametersRcd'),
) 

process.p = cms.Path( process.dbWriterEleScale ) 
