#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "MyAnalysis/IsolationTools/interface/SuperClusterHitsEcalIsolation.h"

SuperClusterHitsEcalIsolation::SuperClusterHitsEcalIsolation(const EcalRecHitCollection *ebRecHits,
                                                             const EcalRecHitCollection *eeRecHits) {

  m_ebRecHits = ebRecHits;
  m_eeRecHits = eeRecHits;
  m_excludeHalo = false;

}

float SuperClusterHitsEcalIsolation::getSum(const edm::Event & iEvent, const edm::EventSetup & iSetup, const reco::SuperCluster *scluster) {

  math::XYZPoint caloPosition = scluster->position();
  const GlobalPoint point(caloPosition.x(), caloPosition.y() , caloPosition.z());

  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* caloGeom = pG.product(); 
  const CaloSubdetectorGeometry* barrelgeom = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  const CaloSubdetectorGeometry* endcapgeom = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

  float sum = 0.0;
  // sum all the rechits in the circle
  // not considering rechits belonging to sc
  std::vector< std::pair<DetId, float> > scHitsByDetId = scluster->hitsAndFractions();
  //  const  std::vector<DetId> scHitsByDetId = scluster->getHitsByDetId();
  if( fabs(scluster->eta()) < 1.479 ) {
    sum = collect(point, barrelgeom, scHitsByDetId, *m_ebRecHits);
  } else {
    sum = collect(point, endcapgeom, scHitsByDetId, *m_eeRecHits);
  }

  return sum;

}

float SuperClusterHitsEcalIsolation::getSum(const edm::Event & iEvent, const edm::EventSetup & iSetup, const reco::BasicCluster *bcluster) {

  math::XYZPoint caloPosition = bcluster->position();
  const GlobalPoint point(caloPosition.x(), caloPosition.y() , caloPosition.z());

  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* caloGeom = pG.product(); 
  const CaloSubdetectorGeometry* barrelgeom = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  const CaloSubdetectorGeometry* endcapgeom = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

  float sum = 0.0;
  // sum all the rechits in the circle
  // not considering rechits belonging to bc
  std::vector< std::pair<DetId, float> > bcHitsByDetId = bcluster->hitsAndFractions();
  if( fabs(bcluster->eta()) < 1.479 ) {
    sum = collect(point, barrelgeom, bcHitsByDetId, *m_ebRecHits);
  } else {
    sum = collect(point, endcapgeom, bcHitsByDetId, *m_eeRecHits);
  }

  return sum;

}

float SuperClusterHitsEcalIsolation::collect( const GlobalPoint &caloPosition, const CaloSubdetectorGeometry* subdet, 
                                              std::vector< std::pair<DetId, float> > scHitsByDetId, const EcalRecHitCollection &hits) {
  
  float energySum = 0.0;
  
  CaloSubdetectorGeometry::DetIdSet chosen = subdet->getCells(caloPosition,m_extRadius);
  EcalRecHitCollection::const_iterator j=hits.end();
  
  for (CaloSubdetectorGeometry::DetIdSet::const_iterator i = chosen.begin(), end = chosen.end() ; i != end;  ++i) {  
    j=hits.find(*i);
    
    bool usedInSuperCluster = false;
    bool scNeighbour = false;

    std::vector<DetId>::iterator scDetId;
    for(std::vector< std::pair<DetId, float> >::iterator scDetId = scHitsByDetId.begin(); scDetId != scHitsByDetId.end(); ++scDetId){
      //    for(scDetId=scHitsByDetId.begin(); scDetId!=scHitsByDetId.end(); ++scDetId) {
      if((*i)==(*scDetId).first) {
        usedInSuperCluster=true;
        break;
      } else {
        std::vector<DetId> neighbours = get3x3(&(*scDetId).first);
        std::vector<DetId>::const_iterator neighbourId;
        for(neighbourId=neighbours.begin(); neighbourId!=neighbours.end(); neighbourId++) {
          if((*i)==(*neighbourId)) {
            scNeighbour=true;
            break;
          }
        }
      }
    }


    if(j != hits.end() && !usedInSuperCluster){

      if( m_excludeHalo && scNeighbour ) {}
      else {
        double energy = j->energy();
        energySum += energy;
      }

    }
  }

  return energySum;

}

std::vector<DetId> SuperClusterHitsEcalIsolation::get3x3(const DetId *id) {

  std::vector<DetId> neighbours;
  for(int icry=0; icry<9; ++icry) {
    unsigned int row    = icry/3;
    unsigned int column = icry%3;
    if(id->subdetId()==EcalBarrel) {
      EBDetId ebid(*id);
      int icryEta = ebid.ieta()+column-1;
      int icryPhi = ebid.iphi()+row-1;
      if ( EBDetId::validDetId(icryEta, icryPhi) ) {
        EBDetId id3x3 = EBDetId(icryEta, icryPhi, EBDetId::ETAPHIMODE);
        neighbours.push_back(id3x3);
      }
    } else if(id->subdetId()==EcalEndcap) {
      EEDetId eeid(*id);
      int icryX = eeid.ix()+column-1;
      int icryY = eeid.iy()+row-1;
      int zside = eeid.zside();
      if ( EEDetId::validDetId(icryX, icryY, zside) ) {
        EEDetId id3x3 = EEDetId(icryX, icryY, zside, EEDetId::XYMODE);
        neighbours.push_back(id3x3);
      }
    }
  }
  
  return neighbours;

}
