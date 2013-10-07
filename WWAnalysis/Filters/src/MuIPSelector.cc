#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

//#include "CommonTools/RecoAlgos/interface/MuonSelector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


struct MuIPSelector {
  typedef reco::MuonCollection collection;

  typedef std::vector<const reco::Muon *> container;

  typedef container::const_iterator const_iterator;

  MuIPSelector ( const edm::ParameterSet & cfg ) :
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

    //if(vertices->size()>1) cout << "WARNING: there is more than 1 primary vertex" << endl;

    Vertex::Point pvPos = vertices->front().position();

    selected_.clear();
    for( reco::MuonCollection::const_iterator mu = c->begin(); 
         mu != c->end(); ++ mu )
      if ( mu->type() != 8         &&
	   fabs( mu->innerTrack()->dxy(pvPos) ) < d0Cut_
	   ) selected_.push_back( & * mu );
	
      
  }
  size_t size() const { return selected_.size(); }

private:
  container selected_;
  double d0Cut_;
};






// define your producer name
typedef ObjectSelector<MuIPSelector> MuIPSelection;

// declare the module as plugin
DEFINE_FWK_MODULE( MuIPSelection );
