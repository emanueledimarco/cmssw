#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/EleTrackerIsolationAlgo.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWElectronIsolator.h"


HWWElectronIsolator::HWWElectronIsolator(const edm::ParameterSet& iConfig)
{
  theTrackIsolCut_	     = iConfig.getParameter<double>("trackIsolCut");
  absolute_                  = iConfig.getParameter<bool>("absolute");
}


HWWElectronIsolator::~HWWElectronIsolator()
{
}

void HWWElectronIsolator::select(edm::Handle<reco::GsfElectronCollection> electrons,
                                 const edm::Event& iEvent,
		                 const edm::EventSetup& iEventSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  if ( !absolute_ && theTrackIsolCut_ > 1.0 ) {
    LogWarning("HWWElectronIsolator") << "Requested relative electron tracker isolation with cut: " 
                                      << theTrackIsolCut_ << " > 1.0. Possible misconfiguration...";
  }

  selected_.clear();
  
  for(unsigned i =0; i<electrons->size(); i++) {
    
    Ref<reco::GsfElectronCollection> electronRef(electrons,i);
    double pt = electronRef->pt();
    double sumPt = electronRef->dr04TkSumPt();
    
    if ( !absolute_ && pt != 0 ) sumPt = sumPt/pt;
    
    if ( sumPt < theTrackIsolCut_ ) selected_.push_back(electronRef);

  }

}
