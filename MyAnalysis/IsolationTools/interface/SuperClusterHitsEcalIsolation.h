#ifndef SuperClusterHitsEcalIsolation_h
#define SuperClusterHitsEcalIsolation_h

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include <vector>

class SuperClusterHitsEcalIsolation {
 public:

  //! ctor
  SuperClusterHitsEcalIsolation() { };
  SuperClusterHitsEcalIsolation(const EcalRecHitCollection *ebRecHits,
                                const EcalRecHitCollection *eeRecHits);

  //! dtor
  ~SuperClusterHitsEcalIsolation() { }

  //! the radius of the external cone 
  void setExtRadius (float extRadius) { m_extRadius = extRadius; }

  //! exclude a halo of 1 additional crystal strip around the super cluster from energy sum
  void excludeHalo(bool what) { m_excludeHalo = what; }

  //! get the sum of the energies in a SuperCluster
  float getSum(const edm::Event & iEvent, const edm::EventSetup & iSetup, const reco::SuperCluster *scluster);

  //! get the sum of the energies in a BasicCluster
  float getSum(const edm::Event & iEvent, const edm::EventSetup & iSetup, const reco::BasicCluster *bcluster);

 private:

  float collect( const GlobalPoint &caloPosition, const CaloSubdetectorGeometry* subdet,
                 const  std::vector< std::pair<DetId, float> > scHitsByDetId, const EcalRecHitCollection &hits);

  std::vector<DetId> get3x3(const DetId *id);

  const EcalRecHitCollection *m_ebRecHits;
  const EcalRecHitCollection *m_eeRecHits;

  float m_extRadius;
  bool m_excludeHalo;

};

#endif
