#ifndef WWAnalysis_AnalysisStep_StringComparator_h
#define WWAnalysis_AnalysisStep_StringComparator_h


#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/Utils/src/ExpressionPtr.h"
#include "CommonTools/Utils/src/ExpressionBase.h"
#include "CommonTools/Utils/interface/expressionParser.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "Reflex/Object.h"
#include <string>

template<typename T, bool DefaultLazyness=false>
struct StringComparator {
  StringComparator(std::string expr, bool sortAscending, bool lazy=DefaultLazyness) : 
    type_(Reflex::Type::ByTypeInfo(typeid(T))),
    sortAscending_(sortAscending) {
    if(! reco::parser::expressionParser<T>(expr, expr_, lazy)) {
      throw edm::Exception(edm::errors::Configuration,
                           "failed to parse \"" + expr + "\"");
    }
  }
  StringComparator(const reco::parser::ExpressionPtr & expr) : 
    expr_(expr),
    type_(Reflex::Type::ByTypeInfo(typeid(T))) {
  }
  double operator()(const T & t1, const T & t2) const {
    using namespace Reflex;
    Object o1(type_, const_cast<T *>(& t1));
    Object o2(type_, const_cast<T *>(& t2));
    if (sortAscending_) return expr_->value(o1) < expr_->value(o2);
    else return expr_->value(o1) > expr_->value(o2);          
  }

private:
  reco::parser::ExpressionPtr expr_;
  Reflex::Type type_;
  bool sortAscending_;
};

namespace reco {
  namespace modules {
    
    template<typename T, bool Lazy>
    struct ParameterAdapter<StringComparator<T, Lazy> > {
      static StringComparator<T, Lazy> make( const edm::ParameterSet & cfg ) {
        return StringComparator<T, Lazy>( cfg.template getParameter<std::string>("sortBy"), cfg.template getParameter<bool>("sortAscending") );
      }
    };
    
  }
}

#endif
