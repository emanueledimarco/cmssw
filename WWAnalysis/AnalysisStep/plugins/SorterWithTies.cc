#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

template<typename Object>
class SorterWithTies : public edm::EDProducer {
    public:
        explicit SorterWithTies(const edm::ParameterSet&);
        ~SorterWithTies() {}

        typedef StringObjectFunction<Object,true> Func;

    private:
        edm::InputTag src_;
        Func sortBy_;
        bool sortAscending_;
        uint32_t maxNumber_;        
        bool keepTies_;

        virtual void produce(edm::Event&, const edm::EventSetup&);
};



template<typename Object>
SorterWithTies<Object>::SorterWithTies(const edm::ParameterSet& iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    sortBy_(iConfig.getParameter<std::string>("sortBy")),
    sortAscending_(iConfig.getParameter<bool>("sortAscending")),
    maxNumber_(iConfig.existsAs<uint32_t>("maxNumber") ? iConfig.getParameter<uint32_t>("maxNumber") : 999999),
    keepTies_(iConfig.existsAs<bool>("keepTies") ? iConfig.getParameter<bool>("keepTies") : false)
{
    produces<std::vector<Object> >();
}

template<typename Object>
void SorterWithTies<Object>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<edm::View<Object> > src;
    iEvent.getByLabel(src_,src);

    typedef std::pair<double, const Object *> ValObj;
    std::vector<ValObj> pairs;
    
    for (typename edm::View<Object>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {    
        double val = sortBy_(*it);
        if (!sortAscending_) val = -val;
        pairs.push_back(std::make_pair(val,&*it));
    }

    std::sort(pairs.begin(), pairs.end());

    std::auto_ptr<std::vector<Object> > out(new std::vector<Object>());
  
    // figure out how many to copy
    size_t n = std::min<size_t>(pairs.size(), maxNumber_);
    if (keepTies_ && n > 0) {
        while (n < pairs.size() && pairs[n].first == pairs[n-1].first) {
           ++n;
        }
    }

    // now I must copy the first n pairs to the output
    out->reserve(n);
    for (size_t i = 0; i < n; ++i) {
        out->push_back(*pairs[i].second);
    }

    iEvent.put(out);
}

#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"
typedef SorterWithTies<reco::SkimEvent4L> SkimEvent4LSorterWithTies;
typedef SorterWithTies<reco::SkimEvent2L> SkimEvent2LSorterWithTies;
DEFINE_FWK_MODULE(SkimEvent2LSorterWithTies);
DEFINE_FWK_MODULE(SkimEvent4LSorterWithTies);
