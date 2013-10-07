#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//#include "CommonTools/RecoAlgos/interface/GsfElectronSelector.h"
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


struct GsfEleIPSelector {
  typedef reco::GsfElectronCollection collection;

  typedef std::vector<const reco::GsfElectron *> container;

  typedef container::const_iterator const_iterator;

  GsfEleIPSelector ( const edm::ParameterSet & cfg ) :
    d0Cut_( cfg.getParameter<double>( "d0Cut" ) )
  { }

  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }


  void select( const edm::Handle<collection> & c, 
	       const edm::Event & event, 
	       const edm::EventSetup& eventSetup) 
  {
    using namespace std;
    using namespace edm;
    using namespace reco;
    
    Handle<VertexCollection> vertices;
    event.getByLabel("offlinePrimaryVertices",vertices);

    Vertex::Point pvPos = vertices->front().position();

    selected_.clear();
    for( reco::GsfElectronCollection::const_iterator ele = c->begin(); 
         ele != c->end(); ++ ele )
      if ( fabs( ele->gsfTrack()->dxy(pvPos) ) < d0Cut_
	   ) selected_.push_back( & * ele );
	
      
  }
  size_t size() const { return selected_.size(); }

private:
  container selected_;
  double d0Cut_;
};






// define your producer name
typedef ObjectSelector<GsfEleIPSelector> GsfEleIPSelection;

// declare the module as plugin
DEFINE_FWK_MODULE( GsfEleIPSelection );
