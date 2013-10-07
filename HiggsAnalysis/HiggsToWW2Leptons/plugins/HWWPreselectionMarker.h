#ifndef HWWPRESELECTIONMARKER
#define HWWPRESELECTIONMARKER

#include <memory>
#include <math.h>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class HWWPreselectionMarker : public edm::EDProducer {
public:
  explicit HWWPreselectionMarker(const edm::ParameterSet&);
  ~HWWPreselectionMarker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  edm::InputTag muonslabel_;
  edm::InputTag electronslabel_;
  edm::InputTag metlabel_;
  edm::InputTag jetslabel_;

  double leptonPtMinMin_;
  double leptonPtMaxMin_;
  double leptonEtaMax_;
  double leptonChargeCombination_;
  double metMin_;
  double invMassMin_;  
  int selectedEvents[7];
};

#endif
