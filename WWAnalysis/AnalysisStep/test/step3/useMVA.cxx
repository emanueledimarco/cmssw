TMVA::Reader * reader_ = 0;
TString algo_; double mvaCut_;

void use(TString dsIn, TString name, int mass=160, TString vars="mll:pt1:pt2:dphill:gammaMRStar", TString mvaAlgo="BDT") {
    TString algo = TString::Format("%s_%s_%d",mvaAlgo.Data(), name.Data(), mass);

    TFile *fIn = TFile::Open(TString::Format("tree_%s.root", dsIn.Data()));    
    TString fileOut = TString::Format("tree_%s.%s.root", dsIn.Data(), algo.Data());
    TFile *fOut = TFile::Open(fileOut, "RECREATE");    

    TMVA::Reader *reader = new TMVA::Reader("V:Color"); 

    float vals[255]; int nvars = 0; char *varname[255];
    for (char *vptr = strtok((char*)vars.Data(), ":"); vptr != 0; vptr = strtok((char*)0,":"), ++nvars) {
        varname[nvars] = vptr;
        std::cout << "  adding variable: " << vptr << std::endl;
        reader->AddVariable(vptr, &vals[nvars]);
    }
    varname[nvars] = NULL;

    reader->BookMVA(algo, "weights/TMVA_"+algo+".weights.xml");

    const char *channel[4] = { "mumu", "elmu", "muel", "elel" };
    for (int c = 0; c < 4; ++c) {
        float mva;
        TTree *tIn  = (TTree *) fIn->Get(TString::Format("%sTree/probe_tree", channel[c]));
        TTree *tOut = new TTree(algo, "");
        tOut->Branch("mva", &mva);
        for (int v = 0; v < nvars; ++v) {
            tIn->SetBranchAddress(varname[v], &vals[v]);
        }
        for (size_t i = 0, n = tIn->GetEntries(); i < n; ++i) {
            tIn->GetEntry(i);
            mva = reader->EvaluateMVA(algo);
            tOut->Fill();
        }
        TDirectory *dOut = fOut->mkdir(TString::Format("%sTree", channel[c]));
        dOut->WriteTObject(tOut);
    }

    fOut->Close();
    fIn->Close();
}
