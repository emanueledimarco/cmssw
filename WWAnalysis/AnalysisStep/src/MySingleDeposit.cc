#include "WWAnalysis/AnalysisStep/interface/MySingleDeposit.h"

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"



#include <boost/regex.hpp>


bool MySingleDeposit::isNumber(const std::string &str) {
   static boost::regex re("^[+-]?(\\d+\\.?|\\d*\\.\\d*)$");
   return regex_match(str.c_str(), re);
}
double MySingleDeposit::toNumber(const std::string &str) {
    return atof(str.c_str());
}

MySingleDeposit::MySingleDeposit(const edm::ParameterSet &iConfig) :
  src_(iConfig.getParameter<edm::InputTag>("src")),
  deltaR_(iConfig.getParameter<double>("deltaR")),
  weightExpr_(iConfig.getParameter<std::string>("weight")),
  skipDefaultVeto_(iConfig.getParameter<bool>("skipDefaultVeto")),
  labelName_(iConfig.getParameter<std::string>("label"))
						      //,vetos_(new AbsVetos())
{
  std::string mode = iConfig.getParameter<std::string>("mode");
  if (mode == "sum") mode_ = Sum; 
  else if (mode == "sumRelative") mode_ = SumRelative; 
  else if (mode == "sum2") mode_ = Sum2;                  
  else if (mode == "sum2Relative") mode_ = Sum2Relative;
  else if (mode == "max") mode_ = Max;                  
  else if (mode == "maxRelative") mode_ = MaxRelative;
  else if (mode == "nearestDR") mode_ = NearestDR;
  else if (mode == "count") mode_ = Count;
  else throw cms::Exception("Not Implemented") << "Mode '" << mode << "' not implemented. " <<
    "Supported modes are 'sum', 'sumRelative', 'count'." << 
    //"Supported modes are 'sum', 'sumRelative', 'max', 'maxRelative', 'count'." << // TODO: on request only
    "New methods can be easily implemented if requested.";
  typedef std::vector<std::string> vstring;
  vstring vetos = iConfig.getParameter< vstring >("vetos");
  reco::isodeposit::EventDependentAbsVeto *evdep; 
  for (vstring::const_iterator it = vetos.begin(), ed = vetos.end(); it != ed; ++it) {
    vetos_.push_back(IsoDepositVetoFactory::make(it->c_str(), evdep));
    if (evdep) evdepVetos_.push_back(evdep);
  }
  std::string weight = iConfig.getParameter<std::string>("weight");
  if (isNumber(weight)) {
    //std::cout << "Weight is a simple number, " << toNumber(weight) << std::endl;
    weight_ = toNumber(weight);
    usesFunction_ = false;
  } else {
    usesFunction_ = true;
    //std::cout << "Weight is a function, this might slow you down... " << std::endl;
  }
  //std::cout << "CandIsolatorFromDeposits::MySingleDeposit::MySingleDeposit: Total of " << vetos_.size() << " vetos" << std::endl;
}
void MySingleDeposit::cleanup() {
    for (reco::IsoDeposit::AbsVetos::iterator it = vetos_.begin(), ed = vetos_.end(); it != ed; ++it) {
        delete *it;
    }
    vetos_.clear();
    // NOTE: we DON'T have to delete the evdepVetos_, they have already been deleted above. We just clear the vectors
    evdepVetos_.clear();
}
void MySingleDeposit::open(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    iEvent.getByLabel(src_, hDeps_);
    for (reco::isodeposit::EventDependentAbsVetos::iterator it = evdepVetos_.begin(), ed = evdepVetos_.end(); it != ed; ++it) {
        (*it)->setEvent(iEvent,iSetup);
    }
}

double MySingleDeposit::compute(const reco::CandidateBaseRef &cand) {
    const reco::IsoDeposit &dep = (*hDeps_)[cand];
    double eta = dep.eta(), phi = dep.phi(); // better to center on the deposit direction
                                             // that could be, e.g., the impact point at calo
    for (reco::IsoDeposit::AbsVetos::iterator it = vetos_.begin(), ed = vetos_.end(); it != ed; ++it) {
        (*it)->centerOn(eta, phi);
    }
    double weight = (usesFunction_ ? weightExpr_(*cand) : weight_);
    switch (mode_) {
        case Count:        return weight * dep.countWithin(deltaR_, vetos_, skipDefaultVeto_);
        case Sum:          return weight * dep.sumWithin(deltaR_, vetos_, skipDefaultVeto_);
        case SumRelative:  return weight * dep.sumWithin(deltaR_, vetos_, skipDefaultVeto_) / dep.candEnergy() ;
        case Sum2:         return weight * dep.sum2Within(deltaR_, vetos_, skipDefaultVeto_);
        case Sum2Relative: return weight * dep.sum2Within(deltaR_, vetos_, skipDefaultVeto_) / (dep.candEnergy() * dep.candEnergy()) ;
        case Max:          return weight * dep.maxWithin(deltaR_, vetos_, skipDefaultVeto_);
//         case NearestDR:    return weight * dep.nearestDR(deltaR_, vetos_, skipDefaultVeto_);
        case NearestDR:    return -9999;
        case MaxRelative:  return weight * dep.maxWithin(deltaR_, vetos_, skipDefaultVeto_) / dep.candEnergy() ;
    }
    throw cms::Exception("Logic error") << "Should not happen at " << __FILE__ << ", line " << __LINE__; // avoid gcc warning
}

