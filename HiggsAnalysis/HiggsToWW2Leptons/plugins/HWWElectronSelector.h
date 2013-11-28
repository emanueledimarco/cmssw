#ifndef HWWELECTRONSELECTOR
#define HWWELECTRONSELECTOR

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

class HWWElectronSelector  {
 public:
  explicit HWWElectronSelector(const edm::ParameterSet&);
  ~HWWElectronSelector();

  // Collections to be selected
  typedef reco::GsfElectronCollection collection;
  typedef std::vector<reco::GsfElectronRef> ::const_iterator const_iterator;

  //define iterators with above typedef
  const_iterator begin () const { return selected_.begin () ; }
  const_iterator end () const { return  selected_.end () ; }
 
  void select (edm::Handle<reco::GsfElectronCollection>,
               const edm::Event&, 
               const edm::EventSetup&) ;

  // ----------member data ---------------------------
 private:

  double hOverECut_;
  double sigmaiEtaiEtaEBCut_;
  double sigmaiEtaiEtaEECut_;
  double deltaphiinCut_;
  double deltaetainCut_;

  std::vector<reco::GsfElectronRef> selected_ ;

};

#endif
