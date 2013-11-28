
#ifndef CalibratedPatElectronProducerCalifornia_h
#define CalibratedPatElectronProducerCalifornia_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

class CalibratedPatElectronProducerCalifornia: public edm::EDProducer 
 {
  public:

    //static void fillDescriptions( edm::ConfigurationDescriptions & ) ;

    explicit CalibratedPatElectronProducerCalifornia( const edm::ParameterSet & ) ;
    virtual ~CalibratedPatElectronProducerCalifornia();
    virtual void produce( edm::Event &, const edm::EventSetup & ) ;

  private:

    edm::InputTag inputPatElectrons ;
    std::string dataset ;
    bool isAOD ;
    bool isMC ;
    bool updateEnergyError ;
    double smearingRatio;
    bool debug ;
    uint energyMeasurementType;
    
 } ;

#endif
