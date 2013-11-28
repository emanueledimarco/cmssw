#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"


#include "WWAnalysis/TreeModifiers/interface/BranchAdderFactory.h"


#include "TFile.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TClass.h"
#include "TCollection.h"
#include "TKey.h"
#include "TTree.h"

#include <sstream>
#include <vector>





int main(int argc,char* argv[]) {


  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  // ----------- READING THE CONFIGURATION FROM THE INPUT PYTHON FILE ---------------------
  PythonProcessDesc builder(argv[1]);
  boost::shared_ptr<edm::ParameterSet > pPset = builder.processDesc()->getProcessPSet();
  const edm::ParameterSet&  genPSet = pPset->getParameter<edm::ParameterSet>("General");
  
  using namespace std;  
  vector<string> inputFiles( genPSet.getParameter<vector<string> >("inputFiles") );
  vector<string> branchAdderNames( genPSet.getParameter<vector<string> >("branchAdders") );



  // ----------- INIZIALIZING THE BRANCHADDERS -----------
  vector<BranchAdder*> branchAdders;
  BranchAdderFactory factory;
  for(unsigned int i=0; i<branchAdderNames.size(); ++i){
    const edm::ParameterSet& specificPSet = pPset->getParameter<edm::ParameterSet>(branchAdderNames[i]);
    branchAdders.push_back( factory.getBranchAdder(branchAdderNames[i],specificPSet) );
  }




  // ----------- LOOPING OVER FILES AND TREES ------------

  //loop over files
  for(unsigned int i=0; i!=inputFiles.size(); ++i){
    string fileName = inputFiles[i];
    //extracting the sample ID from the file name
    int start, length;
    start = fileName.find("_id")+3;
    length = fileName.find(".root") - start;
    int id = atoi(fileName.substr(start,length).c_str());

    cout << "processing " << fileName << " that has sample ID = " << id << endl;
    TFile* file = new TFile(fileName.c_str());
    
    if (file->IsZombie()) {
      continue;
    }

    file->ReOpen("update");

    //loop over branchAdders
    for(unsigned int i=0; i<branchAdderNames.size(); ++i){
      BranchAdder* branchAdder  = branchAdders[i];
      cout << "  using branchAdder " << branchAdderNames[i] << endl; 

      //here I loop over all the directories in the file
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TDirectory")) continue;
	TDirectory *dir = (TDirectory*)key->ReadObj();
	string dirName = dir->GetName();

	//check if this directory is supposed to be processed
	vector<string> treesToProcess = branchAdder->treesToProcess();
	if(treesToProcess.empty() ||  
	   std::find(treesToProcess.begin(),treesToProcess.end(),dirName) != treesToProcess.end()){      
	  file->cd(dirName.c_str());
	  cout << "    processing directory " << dirName << endl;
	  //dir->ls();
	  TTree* tree = (TTree*) dir->Get("probe_tree");
	  if(tree) {
	    //cout << "   found probe_tree for directory " << dirName << endl;
	    branchAdder->addBranch(*tree,id);	  
	  }else{
	    cout << "    WARNING: prob_tree not found in current directory" << endl;
	  }
	}else{
	  //cout << "skipped directory " << dirName << endl;
	}

	dir->Write(0,TObject::kWriteDelete);
      }//end loop over file keys
    }//end loop over branchAdders

    file->Close();
  }// end loop files


  return 0;
}
