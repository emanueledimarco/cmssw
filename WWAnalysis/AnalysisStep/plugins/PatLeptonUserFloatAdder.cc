// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include <DataFormats/Common/interface/ValueMap.h>
#include <DataFormats/Common/interface/View.h>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <vector>

template<typename T>
class PatLeptonUserFloatAdder : public edm::EDProducer {
    public:
        explicit PatLeptonUserFloatAdder(const edm::ParameterSet&);
        ~PatLeptonUserFloatAdder();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag leptonTag_;
        std::vector<std::pair<std::string,edm::InputTag> > mapTags_;

        typedef StringObjectFunction<reco::Candidate> Func;
        typedef StringCutObjectSelector<reco::Candidate> Cut;
        std::vector<std::pair<std::string,Func> > mapFuncs_;
        std::vector<std::pair<std::string,Cut> > mapCuts_;

};

template<typename T>
PatLeptonUserFloatAdder<T>::PatLeptonUserFloatAdder(const edm::ParameterSet& iConfig) :
    leptonTag_(iConfig.getParameter<edm::InputTag>("src")) 
{
    if (iConfig.existsAs<edm::ParameterSet>("valueMaps")) {
        edm::ParameterSet pset = iConfig.getParameter<edm::ParameterSet>("valueMaps");
        std::vector<std::string> names = pset.getParameterNamesForType<edm::InputTag>();
        mapTags_.reserve(names.size());
        for (unsigned int i = 0, n = names.size(); i < n; ++i) {
            mapTags_.push_back(std::make_pair(names[i], pset.getParameter<edm::InputTag>(names[i])));
        }
    }
    if (iConfig.existsAs<edm::ParameterSet>("flags")) {
        edm::ParameterSet pset = iConfig.getParameter<edm::ParameterSet>("flags");
        std::vector<std::string> names = pset.getParameterNamesForType<std::string>();
        mapCuts_.reserve(names.size());
        for (unsigned int i = 0, n = names.size(); i < n; ++i) {
            mapCuts_.push_back(std::make_pair(names[i], Cut(pset.getParameter<std::string>(names[i]), true)));
        }
    }
    if (iConfig.existsAs<edm::ParameterSet>("variables")) {
        edm::ParameterSet pset = iConfig.getParameter<edm::ParameterSet>("variables");
        std::vector<std::string> names = pset.getParameterNamesForType<std::string>();
        mapFuncs_.reserve(names.size());
        for (unsigned int i = 0, n = names.size(); i < n; ++i) {
            mapFuncs_.push_back(std::make_pair(names[i], Func(pset.getParameter<std::string>(names[i]), true)));
        }
    }
   
    produces<std::vector<T> >();  
}


template<typename T>
PatLeptonUserFloatAdder<T>::~PatLeptonUserFloatAdder() {}


// ------------ method called to produce the data  ------------
template<typename T>
void PatLeptonUserFloatAdder<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<T> > leptons;
    iEvent.getByLabel(leptonTag_,leptons);

    std::vector<edm::Handle<edm::ValueMap<float> > > maps_(mapTags_.size());
    for (unsigned int i = 0, n = mapTags_.size(); i < n; ++i) {
        iEvent.getByLabel(mapTags_[i].second, maps_[i]);
    }

    std::auto_ptr<std::vector<T> > out(new std::vector<T>());
    out->reserve(leptons->size());
    for (unsigned int j = 0, m = leptons->size(); j < m; ++j) {
        edm::Ptr<T> lref(leptons->ptrAt(j));
        out->push_back((*leptons)[j]); 
        T &lep = out->back();
        for (unsigned int i = 0, n = mapTags_.size(); i < n; ++i) {
            lep.addUserFloat( mapTags_[i].first, (*maps_[i])[lref] );
        }
        for (unsigned int i = 0, n = mapFuncs_.size(); i < n; ++i) {
            lep.addUserFloat( mapFuncs_[i].first, mapFuncs_[i].second(lep) );
        }
        for (unsigned int i = 0, n = mapCuts_.size(); i < n; ++i) {
            lep.addUserInt( mapCuts_[i].first, mapCuts_[i].second(lep) );
        }
    }
    iEvent.put(out);
}

template<typename T>
void PatLeptonUserFloatAdder<T>::beginJob() { }
template<typename T>
void PatLeptonUserFloatAdder<T>::endJob() { } 

typedef PatLeptonUserFloatAdder<pat::Muon>     PatMuonUserFloatAdder;
typedef PatLeptonUserFloatAdder<pat::Electron> PatElectronUserFloatAdder;
typedef PatLeptonUserFloatAdder<pat::PFParticle> PatPFParticleUserFloatAdder;
DEFINE_FWK_MODULE(PatMuonUserFloatAdder);
DEFINE_FWK_MODULE(PatElectronUserFloatAdder);
DEFINE_FWK_MODULE(PatPFParticleUserFloatAdder);
