
#ifndef RegressionEnergyPatElectronProducerCalifornia_h
#define RegressionEnergyPatElectronProducerCalifornia_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEnergyRegressionEvaluate.h"

class RegressionEnergyPatElectronProducerCalifornia: public edm::EDProducer 
 {
  public:

    explicit RegressionEnergyPatElectronProducerCalifornia( const edm::ParameterSet & ) ;
    virtual ~RegressionEnergyPatElectronProducerCalifornia();
    virtual void produce( edm::Event &, const edm::EventSetup & ) ;

  private:

    edm::InputTag inputPatElectrons ;
    edm::InputTag rhoInputTag_ ;
    uint32_t energyRegressionType ;
    std::string regressionInputFile;
    bool debug ;
    ElectronEnergyRegressionEvaluate *regressionEvaluator;

 } ;

#endif
