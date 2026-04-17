import sys
import os.path
import FWCore.ParameterSet.Config as cms

POPULATE_MC = False
FIRST_RUN_DATA = '2'

if POPULATE_MC: suffix = "mc-_timeShifted_1.000000ns"
else: 
    suffix = sys.argv[2]
    suffixorig = sys.argv[3]


print "reading txt file with suffix ",suffixorig
print "writing into ",suffix

process = cms.Process("ProcessOne")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = 'sqlite_file:ecaltemplates_popcon_'+suffix+'.db'
process.CondDBCommon.DBParameters.authenticationPath = '.'
process.CondDBCommon.DBParameters.messageLevel=cms.untracked.int32(1)

process.MessageLogger = cms.Service("MessageLogger",
                                    debugModules = cms.untracked.vstring('*'),
                                    destinations = cms.untracked.vstring('cout')
                                    )

process.source = cms.Source("EmptyIOVSource",
                            firstValue = cms.uint64(1),
                            lastValue = cms.uint64(1),
                            timetype = cms.string('runnumber'),
                            interval = cms.uint64(1)
                            )

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    logconnect = cms.untracked.string('sqlite_file:logecaltemplates_popcon_'+suffix+'.db'),
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('EcalPh2CubicPulseShapesRcd'),
        tag = cms.string('EcalPh2CubicPulseShapes_data')
    ))
)

txtfile = "atemplate_histograms_ECAL_"+suffixorig+".txt"
if os.path.isfile(txtfile)==False:
    print "WARNING: file ",txtfile," does not exist. Exiting... "
    exit

process.Test1 = cms.EDAnalyzer("ExTestEcalPh2CubicPulseShapesAnalyzer",
    SinceAppendMode = cms.bool(True),
    record = cms.string('EcalPh2CubicPulseShapesRcd'),
    loggingOn = cms.untracked.bool(True),
    Source = cms.PSet(
        firstRun = cms.string('1' if POPULATE_MC else FIRST_RUN_DATA),
        inputFileName = cms.string(txtfile),
        EBCubicPulseShapeTemplate = cms.vdouble (
            0, 0.736756, 1, 0.891306, 0.685984, 0.489531, 0.333569, 0.220315, 0.142279, 0.0903356, 0.0565993, 0.0350856 # 1.0 ns
            ) ,
        EECubicPulseShapeTemplate = cms.vdouble (
            0.077564, 0.774131, 1, 0.916252, 0.723221, 0.524956, 0.361084, 0.239229, 0.154177, 0.097279, 0.0603566, 0.0369412 # 1.0 ns
            )
        )
)


process.p = cms.Path(process.Test1)
