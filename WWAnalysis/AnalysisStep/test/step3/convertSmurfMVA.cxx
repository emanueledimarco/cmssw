// for x in 115 `seq 120 10 200` `seq 250 50 600`; do for y in Av2 A; do for z in SingleMuon DoubleMuon DoubleElectron MuEG; do  root -b -l -q "useMVA.cxx(\"${z}2011$y\",\"smooth\",$x)" > out/$x.$y.$z.out 2>&1 & done; done; done

TMVA::Reader * reader_ = 0;
TString algo_; double mvaCut_;

float convertType(float &lat) {
// smurf:  mumu = 0, elel, elmu, muel 
// latino: elel = 3, elmu, muel, mumu
    return (((int)lat)+2)%4;
//     switch(lat) {
//         case 3: return 1; break;
//         case 4: return 2; break;
//         case 5: return 3; break;
//         case 6: return 0; break;
//     }
//     return -1;
}

int useMVA(TString dsIn, TString name, int mass=160, TString mvaAlgo="BDTG") {

    TString algo = TString::Format("%s_%s_%d",mvaAlgo.Data(), name.Data(), mass);

    TFile *fIn = TFile::Open(TString::Format("../yields/data/tree_%s.root", dsIn.Data()));    
    TString fileOut = TString::Format("trees/tree_%s.%s.root", dsIn.Data(), algo.Data());
    TFile *fOut = TFile::Open(fileOut, "RECREATE");    

    TMVA::Reader *r0j = new TMVA::Reader("V:Color"); 
    TMVA::Reader *r1j = new TMVA::Reader("V:Color"); 

    float smurfType, latinoType;
    float lep1pt, lep2pt, dPhi, dR, dilmass, mt;
    float dPhiDiLepMET, dPhiDiLepJet1;
    float njet;

    r0j->AddVariable("lep1pt",  &lep1pt);   
    r0j->AddVariable("lep2pt",  &lep2pt);   
    r0j->AddVariable("dPhi",    &dPhi);     
    r0j->AddVariable("dR",      &dR);       
    r0j->AddVariable("dilmass", &dilmass);  
    r0j->AddVariable("type",    &smurfType);
    r0j->AddVariable("mt",      &mt);       
    r1j->AddVariable("lep1pt",  &lep1pt);   
    r1j->AddVariable("lep2pt",  &lep2pt);   
    r1j->AddVariable("dPhi",    &dPhi);     
    r1j->AddVariable("dR",      &dR);       
    r1j->AddVariable("dilmass", &dilmass);  
    r1j->AddVariable("type",    &smurfType);
    r1j->AddVariable("mt",      &mt);       
    r1j->AddVariable("dPhiDiLepMET",      &dPhiDiLepMET); 
    r1j->AddVariable("dPhiDiLepJet1",     &dPhiDiLepJet1);

    r0j->BookMVA(algo, TString::Format("weights/ntuples_%dtrain_0jets_%s.weights.xml",mass,mvaAlgo.Data()));
    r1j->BookMVA(algo, TString::Format("weights/ntuples_%dtrain_1jets_%s.weights.xml",mass,mvaAlgo.Data()));

    const char *channel[4] = { "mumu", "elmu", "muel", "elel" };
    for (int c = 0; c < 4; ++c) {
        float mva;
        TTree *tIn  = (TTree *) fIn->Get(TString::Format("%sTree/probe_tree", channel[c]));
        TTree *tOut = new TTree(algo, "");
        tOut->Branch("mva", &mva);
        tIn->SetBranchAddress("njet",          &njet);         
        tIn->SetBranchAddress("pt1",           &lep1pt);         
        tIn->SetBranchAddress("pt2",           &lep2pt);         
        tIn->SetBranchAddress("dphill",        &dPhi);           
        tIn->SetBranchAddress("drll",          &dR);             
        tIn->SetBranchAddress("mll",           &dilmass);        
        tIn->SetBranchAddress("hypo",          &latinoType);      
        tIn->SetBranchAddress("mth",           &mt);             
        tIn->SetBranchAddress("dphillmet",     &dPhiDiLepMET);  
        tIn->SetBranchAddress("dphilljet",     &dPhiDiLepJet1); 

        for (size_t i = 0, n = tIn->GetEntries(); i < n; ++i) {
            tIn->GetEntry(i);
            smurfType = convertType(latinoType);
            if      (njet == 0) mva = r0j->EvaluateMVA(algo);
            else if (njet == 1) mva = r1j->EvaluateMVA(algo);
            else                mva = -9999;
            tOut->Fill();
        }
        TDirectory *dOut = fOut->mkdir(TString::Format("%sTree", channel[c]));
        dOut->WriteTObject(tOut);
    }

    fOut->Close();
    fIn->Close();
}
