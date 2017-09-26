import FWCore.ParameterSet.Config as cms 
import sys
process = cms.Process('egmscaledb') 

FIRSTRUN=int(sys.argv[2])
TXTFILE=str(sys.argv[3])
DBFILE=str(sys.argv[4])
TAG=str(sys.argv[5])

process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    lastValue = cms.uint64(1),
    interval = cms.uint64(1)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1)) 

process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = 'sqlite_file:'+DBFILE

process.PoolDBOutputService = cms.Service('PoolDBOutputService', 
   process.CondDB,
   timetype = cms.untracked.string('runnumber'),
   toPut = cms.VPSet( 
      cms.PSet(
         record = cms.string('EgmCorrectorParametersRcd'),
         tag    = cms.string(TAG), 
      ),
   ) 
) 

process.dbWriterEleScale = cms.EDAnalyzer('EgmCorrectorDBWriter', 
   txtfile   = cms.untracked.string(TXTFILE),
   payloadTag = cms.untracked.string(TAG),
   record = cms.untracked.string('EgmCorrectorParametersRcd'),
   startIOVRun = cms.untracked.uint32(FIRSTRUN),
) 

process.p = cms.Path( process.dbWriterEleScale ) 
