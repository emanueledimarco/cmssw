#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//#include "CommonTools/RecoAlgos/interface/GsfElectronSelector.h"
//#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>
#include <DataFormats/PatCandidates/interface/Electron.h>


#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

//BDT ElectronID
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"


struct PatEleBDTSelector {
  typedef pat::ElectronCollection collection;

  typedef std::vector<const pat::Electron *> container;

  typedef container::const_iterator const_iterator;

  PatEleBDTSelector ( const edm::ParameterSet & cfg ) 
  {
    eleMVA  = new ElectronIDMVA();
    const char *base=getenv("CMSSW_BASE");
    std::string baseFolder(base);
    baseFolder += "/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/"; 
    eleMVA->Initialize("BDTG method",
		       (baseFolder+"Subdet0LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet1LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet2LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet0HighPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet1HighPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet2HighPt_WithIPInfo_BDTG.weights.xml").c_str(),                
		       ElectronIDMVA::kWithIPInfo);
  }

  ~PatEleBDTSelector ( ) 
  {
    delete eleMVA;
  }

  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }


  void select( const edm::Handle<collection> & c, 
	       const edm::Event & event, 
	       const edm::EventSetup& eventSetup) 
  {
    using namespace std;
    using namespace edm;
    using namespace reco;
    

    selected_.clear();
    for( collection::const_iterator ele = c->begin(); 
         ele != c->end(); ++ ele ){
      bool selection = false;
      double EleOneOverEMinusOneOverP = (1.0/(ele->superCluster()->energy())) - 1.0 / ele->gsfTrack()->p(); 
      double xieSign  = ( (-ele->userFloat("dxyPV")) >=0 )  ? 1: -1;
      double mvaValue = eleMVA->MVAValue(ele->pt() , 
					 ele->superCluster()->eta(), 
					 ele->sigmaIetaIeta(), 
					 ele->deltaEtaSuperClusterTrackAtVtx(),
					 ele->deltaPhiSuperClusterTrackAtVtx(),
					 ele->hcalOverEcal(),
					 -ele->userFloat("dxyPV"),
					 ele->fbrem(), 
					 ele->eSuperClusterOverP(), 
					 ele->eSeedClusterOverPout(),
					 ele->userFloat("sigmaIphiIphi"),
					 ele->basicClustersSize() - 1,
					 EleOneOverEMinusOneOverP,
					 ele->eSeedClusterOverP(),
					 xieSign*ele->userFloat("ip"),
					 xieSign*ele->userFloat("ip")/ele->userFloat("ipErr"));
      //here I have to define the cut on the mvaValue;
      float pt = ele->pt();
      float eta = ele->superCluster()->eta();
      if(pt<20. && fabs(eta)<=1.0) selection = mvaValue > 0.139;
      else if(pt<20. && fabs(eta) > 1.0 && fabs(eta) <= 1.479) selection = mvaValue > 0.525;
      else if(pt<20. && fabs(eta) > 1.479 && fabs(eta) <= 2.5) selection = mvaValue > 0.543;
      else if(pt>=20. && fabs(eta) > 0.0 && fabs(eta) <= 1.0) selection = mvaValue > 0.947;
      else if(pt>=20. && fabs(eta) > 1.0 && fabs(eta) <= 1.479) selection = mvaValue > 0.950;
      else if(pt>=20. && fabs(eta) > 1.479 && fabs(eta) <= 2.5) selection = mvaValue > 0.884;
      else selection = false;

      if ( selection )  selected_.push_back( & * ele );
    }	
      
  }
  size_t size() const { return selected_.size(); }

private:
  container selected_;
  ElectronIDMVA* eleMVA;
};






// define your producer name
  typedef ObjectSelector<PatEleBDTSelector> PatEleBDTSelection;
  
// declare the module as plugin
DEFINE_FWK_MODULE( PatEleBDTSelection );
