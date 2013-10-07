#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "CommonTools/RecoAlgos/interface/GsfElectronSelector.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>

#include <RecoEgamma/EgammaTools/interface/ConversionFinder.h>
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"

template<typename Object,typename OutputObject>
struct ConvRejectionSelector {
  typedef std::vector<Object> collection;

  //typedef std::vector<const reco::GsfElectron *> container;
  typedef std::vector<const OutputObject *> container;

  typedef typename container::const_iterator const_iterator;

  ConvRejectionSelector ( const edm::ParameterSet & cfg ) :
    inputTracks_(cfg.getParameter<edm::InputTag>("srcTracks")),
    distMax_(cfg.getParameter<double>("distMax")),
    dCotanThetaMax_(cfg.getParameter<double>("dCotanThetaMax")) { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }


  void select( const edm::Handle<collection> & c, 
	       const edm::Event & event, 
	       const edm::EventSetup& eventSetup) 
  {
    using namespace std;
    using namespace edm;
    using namespace reco;
    
    Handle<TrackCollection> tracks;
    event.getByLabel(inputTracks_,tracks);
    

    
    selected_.clear();

    for( typename collection::const_iterator ele = c->begin(); 
         ele != c->end(); ++ ele ){
      ConversionFinder convFinder;
      ConversionInfo convInfo = convFinder.getConversionInfo(*ele, tracks, 3.8112);
      //if( convFinder.isFromConversion()  )
      if( fabs(convInfo.dist()) < distMax_ && fabs(convInfo.dcot()) < dCotanThetaMax_  )
	{
	  //the electron is likely coming from a conversion
	}else
	selected_.push_back( & * ele );
    }
	
      
  }
  size_t size() const { return selected_.size(); }

private:
  edm::InputTag inputTracks_;
  double distMax_,dCotanThetaMax_;

  container selected_;
};






// define your producer name
typedef ObjectSelector<ConvRejectionSelector<reco::GsfElectron,reco::GsfElectron> > ConvRejectionSelection;
typedef ObjectSelector<ConvRejectionSelector<pat::Electron,pat::Electron> > ConvRejectionSelectionPAT;

// declare the module as plugin
DEFINE_FWK_MODULE( ConvRejectionSelection );
DEFINE_FWK_MODULE( ConvRejectionSelectionPAT );
