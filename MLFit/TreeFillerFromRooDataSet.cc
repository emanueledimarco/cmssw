#include "TTree.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include <vector>
#include <string>

#include "MLFit/TreeFillerFromRooDataSet.hh"

using namespace std;

TTree* TreeFillerFromRooDataSet::getTree() {
  _tree = new TTree(_dataset->GetName(), "tree with selected variables");
  const RooArgSet *fullSet = _dataset->get();

  int size = _variables.size();
  vector<float> varsVal;
  vector<RooRealVar*> vars;
  varsVal.reserve(size);
  vars.reserve(size);

  // make the necessary branches 
  for(int i=0;i<size;++i) {
    RooRealVar *var = (RooRealVar*)fullSet->find(_variables[i].c_str());
    if(var) {
      cout << "Adding the branch for variable " << var->GetName() << endl;
      std::string typecol(_variables[i]);
      typecol+="/F";
      _tree->Branch(_variables[i].c_str(), &varsVal[i], typecol.c_str());
      vars.push_back(var);
    }
  }

  // loop and fill the tree
  for(int ievt=0; ievt<_dataset->numEntries();++ievt) {
    _dataset->get(ievt);
    for(int ivar=0;ivar<(int)vars.size(); ++ivar) {
      varsVal[ivar] = vars[ivar]->getVal();
    }
    _tree->Fill();
  }
  return _tree;
}
