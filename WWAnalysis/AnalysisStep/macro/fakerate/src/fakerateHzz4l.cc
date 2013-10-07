#include "estimateMuonFakeRateHzz4lTree.C"
#include "estimateElectronFakeRateHzz4lTree.C"

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int main(int argc, char* argv[]) {

  char outname[300];

  if( argc < 2) {
    cout << "missing argument: fakerate <basename>" << endl;
    return 1;
  }
  strcpy(outname,argv[1]);

  cout << "run the FR for non triggering muons..." << endl;
  TFile *file = TFile::Open("results_data/hzzTree_data.root");
  TTree *tree = (TTree*)file->Get("zllmtree/probe_tree");
  estimateMuonFakeRateHzz4lTree analyzer(tree);
  TString outfileUnbias(outname);
  outfileUnbias += TString("_zllm");
  analyzer.Loop(outfileUnbias);
  cout << "DONE unbiased mu." << endl;

  cout << "run the FR for non triggering muons..." << endl;
  TFile *file2 = TFile::Open("results_data/hzzTree_data.root");
  TTree *tree2 = (TTree*)file2->Get("zlletree/probe_tree");
  estimateElectronFakeRateHzz4lTree analyzer2(tree2);
  TString outfileUnbias2(outname);
  outfileUnbias2 += TString("_zlle");
  analyzer2.Loop(outfileUnbias2);
  cout << "DONE unbiased electrons." << endl;

  return 0;
}
