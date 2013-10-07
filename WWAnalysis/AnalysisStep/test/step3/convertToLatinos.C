#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include<iostream>

int convertToLatinos(TString inDir,TString treeName) {

    const char *chans[4] = { "elel", "elmu", "muel", "mumu" };
    for (int c = 0; c < 4; ++c) {
        TFile *fIn = TFile::Open(inDir+"/"+treeName);
        TTree *tIn  = (TTree *) fIn->Get(TString::Format("%sTree/probe_tree",chans[c]));
        TFile *fOut = TFile::Open(TString("./")+chans[c]+"_"+treeName,"RECREATE");
        TTree *tOut = tIn->CloneTree();
        tOut->SetTitle("latino");
        tOut->SetName("latino");
        tOut->Write("latino",TObject::kOverwrite);
        fOut->Close();
    }
    return 0;
}
