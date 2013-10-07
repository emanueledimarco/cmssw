
import FWCore.ParameterSet.Config as cms


#==============================================================================
# regreesion energy gsf electrons
#==============================================================================

regressionEnergyPatElectrons = cms.EDProducer("RegressionEnergyPatElectronProducerCalifornia",
                                              
    # input collections
    inputPatElectronsTag = cms.InputTag("boostedElectrons"),

    # set to True to get debugging printout   
    debug = cms.bool(True),
                                              
    #regression energy type
    energyRegressionType = cms.untracked.uint32(1),

    #regression weight file                                              
    regressionInputFile = cms.string("weightFile_V11.root")
)


