#ifndef ECAL_CUBIC_PULSESHAPES_HANDLER_H
#define ECAL_CUBIC_PULSESHAPES_HANDLER_H

#include <vector>
#include <typeinfo>
#include <string>
#include <map>
#include <iostream>
#include <ctime>

#include "CondCore/PopCon/interface/PopConSourceHandler.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"

#include "CondFormats/EcalObjects/interface/EcalCubicPulseShapeT.h"
#include "CondFormats/DataRecord/interface/EcalCubicPulseShapesRcd.h"

#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/Provenance/interface/Timestamp.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}  // namespace edm

namespace popcon {

template <class P>
  class EcalCubicPulseShapesHandler : public popcon::PopConSourceHandler<P> {
  public:
    EcalCubicPulseShapesHandler(edm::ParameterSet const&);
    ~EcalCubicPulseShapesHandler() override;
    bool checkPulseShape(EcalCubicPulseShapeT<P>::Item* item);
    void fillSimPulseShape(EcalCubicPulseShapeT<P>::Item* item, bool isbarrel);
    void getNewObjects() override;
    std::string id() const override { return m_name; }

  private:
    const EcalCubicPulseShapeT<P>* mypulseshapes;

    unsigned int m_firstRun;
    unsigned int m_lastRun;

    std::string m_gentag;
    std::string m_filename;
    std::string m_name;
    std::vector<double> m_EBPulseShapeTemplate, m_EEPulseShapeTemplate;
  };
}  // namespace popcon
#endif
