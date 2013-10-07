#ifndef WWAnalysis_AnalysisStep_SortSelector_h
#define WWAnalysis_AnalysisStep_SortSelector_h

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StoreContainerTrait.h"
#include "CommonTools/UtilAlgos/interface/SelectionAdderTrait.h"
#include "WWAnalysis/AnalysisStep/interface/StringSortCollectionSelector.h"

template<typename InputCollection, typename Comparator, 
         typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type,
         typename StoreContainer = typename ::helper::StoreContainerTrait<OutputCollection>::type,
         typename PostProcessor = ::helper::NullPostProcessor<OutputCollection>,
         typename StoreManager = typename ::helper::StoreManagerTrait<OutputCollection>::type,
         typename Base = typename ::helper::StoreManagerTrait<OutputCollection>::base,
         typename RefAdder = typename ::helper::SelectionAdderTrait<InputCollection, StoreContainer>::type>
class SortSelector : 
  public ObjectSelector<StringSortCollectionSelector<InputCollection, Comparator, OutputCollection, StoreContainer, RefAdder>, 
                        OutputCollection, NonNullNumberSelector, PostProcessor, StoreManager, Base> {
public:
  explicit SortSelector( const edm::ParameterSet & cfg ) :
    ObjectSelector<StringSortCollectionSelector<InputCollection, Comparator, OutputCollection, StoreContainer, RefAdder>, 
                   OutputCollection, NonNullNumberSelector, PostProcessor>( cfg ) { }
  virtual ~SortSelector() { }
};

#endif


