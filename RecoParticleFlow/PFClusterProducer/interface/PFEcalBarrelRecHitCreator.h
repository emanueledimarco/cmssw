#ifndef RecoParticleFlow_PFClusterProducer_PFEcalRecHitCreator_h
#define RecoParticleFlow_PFClusterProducer_PFEcalRecHitCreator_h

#include "RecoParticleFlow/PFClusterProducer/interface/PFRecHitCreatorBase.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"

#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

template <typename Geometry,PFLayer::Layer Layer,int Detector>
  class PFEcalRecHitCreator :  public  PFRecHitCreatorBase {

 public:  
  PFEcalRecHitCreator(const edm::ParameterSet& iConfig,edm::ConsumesCollector& iC):
    PFRecHitCreatorBase(iConfig,iC)
    {
      recHitToken_ = iC.consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("src"));
    }

    void importRecHits(std::auto_ptr<reco::PFRecHitCollection>&out,std::auto_ptr<reco::PFRecHitCollection>& cleaned ,const edm::Event& iEvent,const edm::EventSetup& iSetup) {

      beginEvent(iEvent,iSetup);
      initSRPMap(iSetup);

      edm::Handle<EcalRecHitCollection> recHitHandle;

      edm::ESHandle<CaloGeometry> geoHandle;
      iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
      // get the ecal geometry
      const CaloSubdetectorGeometry *gTmp = 
	geoHandle->getSubdetectorGeometry(DetId::Ecal, Detector);

      const Geometry *ecalGeo =dynamic_cast< const Geometry* > (gTmp);

      iEvent.getByToken(recHitToken_,recHitHandle);
      for(const auto& erh : *recHitHandle ) {      
	const DetId& detid = erh.detid();
	auto energy = erh.energy();
	auto time = erh.time();
        bool hi = isHighInterest(detid);

	const CaloCellGeometry * thisCell= ecalGeo->getGeometry(detid);
  
	// find rechit geometry
	if(!thisCell) {
	  edm::LogError("PFEcalRecHitCreator")
	    <<"warning detid "<<detid.rawId()
	    <<" not found in geometry"<<std::endl;
	  continue;
	}

	out->emplace_back(thisCell, detid.rawId(),Layer,
			   energy); 

        auto & rh = out->back();
	
	bool rcleaned = false;
	bool keep=true;

	//Apply Q tests
	for( const auto& qtest : qualityTests_ ) {
	  if (!qtest->test(rh,erh,rcleaned,hi)) {
	    keep = false;	    
	  }
	}
	  
	if(keep) {
	  rh.setTime(time);
	  rh.setDepth(1);
	} 
        else {
	  if (rcleaned) 
	    cleaned->push_back(std::move(out->back()));
          out->pop_back();
        }
      }
    }



 protected:

  void initSRPMap(const edm::EventSetup &es);
  bool isHighInterest(const EBDetId& detid);
  bool srFullReadOut(EcalTrigTowerDetId& towid);

  edm::EDGetTokenT<EcalRecHitCollection> recHitToken_;

  const EcalTrigTowerConstituentsMap* eTTmap_;  
  // Array of the DetIds
  std::vector<EcalTrigTowerDetId> theTTDetIds_;
  // neighboring TT DetIds
  std::vector<std::vector<int> > neighboringTTs_;
  // the crystals in a given TT 
  std::vector<std::vector<int> > crystalsinTT_;
  // the towers which have been looked at 
  std::vector<int> theTTofHighInterest_;
  // the status of the towers. A tower is of high interest if it or one of its neighbour is above the threshold
  std::vector<int> TTHighInterest_;

};


typedef PFEcalRecHitCreator<EcalBarrelGeometry,PFLayer::ECAL_BARREL,EcalBarrel> PFEBRecHitCreator;
typedef PFEcalRecHitCreator<EcalEndcapGeometry,PFLayer::ECAL_ENDCAP,EcalEndcap> PFEERecHitCreator;

#endif
