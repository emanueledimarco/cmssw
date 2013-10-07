#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "scales2.h"
#include "WWAnalysis/TreeModifiers/interface/HiggsMassWeightProvider.h"

float getMassCut(float cardmass, bool low) {
    std::map<float, float> higgswidth;
    higgswidth[200.0] = 1.43E+00;
    higgswidth[202.0] = 1.51E+00;
    higgswidth[204.0] = 1.59E+00;
    higgswidth[206.0] = 1.68E+00;
    higgswidth[208.0] = 1.76E+00;
    higgswidth[210.0] = 1.85E+00;
    higgswidth[212.0] = 1.93E+00;
    higgswidth[214.0] = 2.02E+00;
    higgswidth[216.0] = 2.12E+00;
    higgswidth[218.0] = 2.21E+00;
    higgswidth[220.0] = 2.31E+00;
    higgswidth[222.0] = 2.40E+00;
    higgswidth[224.0] = 2.50E+00;
    higgswidth[226.0] = 2.61E+00;
    higgswidth[228.0] = 2.71E+00;
    higgswidth[230.0] = 2.82E+00;
    higgswidth[232.0] = 2.93E+00;
    higgswidth[234.0] = 3.04E+00;
    higgswidth[236.0] = 3.16E+00;
    higgswidth[238.0] = 3.27E+00;
    higgswidth[240.0] = 3.40E+00;
    higgswidth[242.0] = 3.52E+00;
    higgswidth[244.0] = 3.64E+00;
    higgswidth[246.0] = 3.77E+00;
    higgswidth[248.0] = 3.91E+00;
    higgswidth[250.0] = 4.04E+00;
    higgswidth[252.0] = 4.18E+00;
    higgswidth[254.0] = 4.32E+00;
    higgswidth[256.0] = 4.46E+00;
    higgswidth[258.0] = 4.61E+00;
    higgswidth[260.0] = 4.76E+00;
    higgswidth[262.0] = 4.91E+00;
    higgswidth[264.0] = 5.07E+00;
    higgswidth[266.0] = 5.23E+00;
    higgswidth[268.0] = 5.39E+00;
    higgswidth[270.0] = 5.55E+00;
    higgswidth[272.0] = 5.72E+00;
    higgswidth[274.0] = 5.89E+00;
    higgswidth[276.0] = 6.07E+00;
    higgswidth[278.0] = 6.25E+00;
    higgswidth[280.0] = 6.43E+00;
    higgswidth[282.0] = 6.61E+00;
    higgswidth[284.0] = 6.80E+00;
    higgswidth[286.0] = 6.99E+00;
    higgswidth[288.0] = 7.19E+00;
    higgswidth[290.0] = 7.39E+00;
    higgswidth[295.0] = 7.90E+00;
    higgswidth[300.0] = 8.43E+00;
    higgswidth[305.0] = 8.99E+00;
    higgswidth[310.0] = 9.57E+00;
    higgswidth[315.0] = 1.02E+01;
    higgswidth[320.0] = 1.08E+01;
    higgswidth[325.0] = 1.14E+01;
    higgswidth[330.0] = 1.21E+01;
    higgswidth[335.0] = 1.28E+01;
    higgswidth[340.0] = 1.35E+01;
    higgswidth[345.0] = 1.42E+01;
    higgswidth[350.0] = 1.52E+01;
    higgswidth[360.0] = 1.76E+01;
    higgswidth[370.0] = 2.02E+01;
    higgswidth[380.0] = 2.31E+01;
    higgswidth[390.0] = 2.61E+01;
    higgswidth[400.0] = 2.92E+01;
    higgswidth[410.0] = 3.25E+01;
    higgswidth[420.0] = 3.59E+01;
    higgswidth[430.0] = 3.94E+01;
    higgswidth[440.0] = 4.30E+01;
    higgswidth[450.0] = 4.68E+01;
    higgswidth[460.0] = 5.08E+01;
    higgswidth[470.0] = 5.49E+01;
    higgswidth[480.0] = 5.91E+01;
    higgswidth[490.0] = 6.35E+01;
    higgswidth[500.0] = 6.80E+01;
    higgswidth[510.0] = 7.27E+01;
    higgswidth[520.0] = 7.75E+01;
    higgswidth[530.0] = 8.25E+01;
    higgswidth[540.0] = 8.77E+01;
    higgswidth[550.0] = 9.30E+01;
    higgswidth[560.0] = 9.86E+01;
    higgswidth[570.0] = 1.04E+02;
    higgswidth[580.0] = 1.10E+02;
    higgswidth[590.0] = 1.16E+02;
    higgswidth[600.0] = 1.23E+02;
    higgswidth[610.0] = 1.29E+02;
    higgswidth[620.0] = 1.36E+02;
    higgswidth[630.0] = 1.43E+02;
    higgswidth[640.0] = 1.50E+02;
    higgswidth[650.0] = 1.58E+02;
    higgswidth[660.0] = 1.65E+02;
    higgswidth[670.0] = 1.73E+02;
    higgswidth[680.0] = 1.82E+02;
    higgswidth[690.0] = 1.90E+02;
    higgswidth[700.0] = 1.99E+02;
    higgswidth[710.0] = 2.08E+02;
    higgswidth[720.0] = 2.17E+02;
    higgswidth[730.0] = 2.27E+02;
    higgswidth[740.0] = 2.37E+02;
    higgswidth[750.0] = 2.47E+02;
    higgswidth[760.0] = 2.58E+02;
    higgswidth[770.0] = 2.69E+02;
    higgswidth[780.0] = 2.80E+02;
    higgswidth[790.0] = 2.92E+02;
    higgswidth[800.0] = 3.04E+02;
    higgswidth[810.0] = 3.17E+02;
    higgswidth[820.0] = 3.30E+02;
    higgswidth[830.0] = 3.43E+02;
    higgswidth[840.0] = 3.57E+02;
    higgswidth[850.0] = 3.71E+02;
    higgswidth[860.0] = 3.86E+02;
    higgswidth[870.0] = 4.01E+02;
    higgswidth[880.0] = 4.16E+02;
    higgswidth[890.0] = 4.32E+02;
    higgswidth[900.0] = 4.49E+02;
    higgswidth[910.0] = 4.66E+02;
    higgswidth[920.0] = 4.84E+02;
    higgswidth[930.0] = 5.02E+02;
    higgswidth[940.0] = 5.21E+02;
    higgswidth[950.0] = 5.40E+02;
    higgswidth[960.0] = 5.60E+02;
    higgswidth[970.0] = 5.81E+02;
    higgswidth[980.0] = 6.02E+02;
    higgswidth[990.0] = 6.24E+02;
    higgswidth[1000.0] = 6.47E+02;

    double windowVal = max(higgswidth[cardmass], float(1.));
    double lowside = (cardmass >= 275) ? 180. : 100.;
    //if (low) return std::max((cardmass - 20.*windowVal), lowside);
    //else return std::min((cardmass + 15.*windowVal), 2000.);


    if (low) return 100.;
    else return 1600.;
}


std::pair<float, std::pair<float, float> > getYield (std::string filename, float hmass, int ch, bool is7, int id) {

    TFile file(filename.c_str());
    TTree* tree = (TTree*)file.Get("zz4lTree/probe_tree");
   
    TBranch *bmass      = tree->GetBranch("genhiggsmass"); 
    TBranch *bnumsim    = tree->GetBranch("numTrueInteractions");
    TBranch *bl1pt      = tree->GetBranch("l1pt");
    TBranch *bl1eta     = tree->GetBranch("l1eta");
    TBranch *bl1pdgId   = tree->GetBranch("l1pdgId");
    TBranch *bl2pt      = tree->GetBranch("l2pt");
    TBranch *bl2eta     = tree->GetBranch("l2eta");
    TBranch *bl2pdgId   = tree->GetBranch("l2pdgId");
    TBranch *bl3pt      = tree->GetBranch("l3pt");
    TBranch *bl3eta     = tree->GetBranch("l3eta");
    TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
    TBranch *bl4pt      = tree->GetBranch("l4pt");
    TBranch *bl4eta     = tree->GetBranch("l4eta");
    TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");

    float vmass     = 0.0;
    float numsim    = 0.0;
    float l1pt      = 0.0;
    float l1eta     = 0.0;
    float l1pdgId   = 0.0;
    float l2pt      = 0.0;
    float l2eta     = 0.0;
    float l2pdgId   = 0.0;
    float l3pt      = 0.0;
    float l3eta     = 0.0;
    float l3pdgId   = 0.0;
    float l4pt      = 0.0;
    float l4eta     = 0.0;
    float l4pdgId   = 0.0;

    bmass   ->SetAddress(&vmass);
    bnumsim ->SetAddress(&numsim);
    bl1pt   ->SetAddress(&l1pt);
    bl1eta  ->SetAddress(&l1eta);
    bl1pdgId->SetAddress(&l1pdgId);
    bl2pt   ->SetAddress(&l2pt);
    bl2eta  ->SetAddress(&l2eta);
    bl2pdgId->SetAddress(&l2pdgId);
    bl3pt   ->SetAddress(&l3pt);
    bl3eta  ->SetAddress(&l3eta);
    bl3pdgId->SetAddress(&l3pdgId);
    bl4pt   ->SetAddress(&l4pt);
    bl4eta  ->SetAddress(&l4eta);
    bl4pdgId->SetAddress(&l4pdgId);

    Float_t newweight;
    TBranch* newBranch = tree->Branch("newweight", &newweight, "newweight/F");
    Long64_t nentries = tree->GetEntries();

    cout << filename << " " << nentries << endl;

    std::stringstream weightss;
    if (hmass>=400 && (int(hmass))%50 == 0) {
        if (is7)weightss <<"/data1/sixie/releases/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs" << int(hmass) << "_7TeV_Lineshape+Interference.txt"; 
        else    weightss <<"/data1/sixie/releases/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs" << int(hmass) << "_8TeV_Lineshape+Interference.txt"; 
    }    
    else weightss << "";

    HiggsMassWeightProvider hmwp(weightss.str());
    for (unsigned int i = 0; i < nentries; i++){
        bmass   ->GetEvent(i);
        bnumsim ->GetEvent(i);
        bl1pt   ->GetEvent(i);
        bl1eta  ->GetEvent(i);
        bl1pdgId->GetEvent(i);
        bl2pt   ->GetEvent(i);
        bl2eta  ->GetEvent(i);
        bl2pdgId->GetEvent(i);
        bl3pt   ->GetEvent(i);
        bl3eta  ->GetEvent(i);
        bl3pdgId->GetEvent(i);
        bl4pt   ->GetEvent(i);
        bl4eta  ->GetEvent(i);
        bl4pdgId->GetEvent(i);
       
        newweight = 1.0;
        if (hmass>=400 && (int(hmass))%50 == 0) newweight *= hmwp.getWeight(vmass); 
        newweight *= getPUWeight(numsim); 
        newweight *= getSF(l1pt, l1eta, l1pdgId);
        newweight *= getSF(l2pt, l2eta, l2pdgId);
        newweight *= getSF(l3pt, l3eta, l3pdgId);
        newweight *= getSF(l4pt, l4eta, l4pdgId);

        newweight *= is7 ? (getXsecggH(hmass)/getXsecggHByChannel(hmass, ch))/evt_7TeV(id) : (getXsecggH(hmass)/getXsecggHByChannel(hmass, ch))/evt_8TeV(id);
        //newweight *= (getXsecggHByChannel(hmass, ch) + getXsecVBFByChannel(hmass, ch) + getXsecWHiByChannel(hmass, ch) + getXsecZHiByChannel(hmass, ch));

        newBranch->Fill();
    }

    stringstream cutss;
    cutss << "newweight * (z1mass>40 && z1mass<120 && z2mass>12 && z2mass<120 && melaLD > 0.5 && mass>";
    cutss << getMassCut(hmass, true ) << " && mass<";
    cutss << getMassCut(hmass, false);
    if (ch == 0 || ch == 1) cutss << " && " << "channel == " << ch << ")";
    else                    cutss << " && " << "(channel == 2 || channel == 3))";

    TH1F* hist = new TH1F("hist", "", 1, 100., 2000.);
    hist->Sumw2();
    

    tree->Draw("mass>>hist", cutss.str().c_str());

    return std::pair<float, std::pair<float, float> >(hmass, std::pair<float, float>(hist->Integral(), hist->GetBinError(1)));
}

void yieldFit() {

  int channel = 2;

  bool do7TeV = false;

    if (do7TeV == false) {
        init(do7TeV);
        
        std::string treeFolder = "/data1/sixie/ntuples/step2/hcp/2012/";
        
        std::vector<std::pair<float, std::pair<float, float> > > massyieldpairs;
        
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1115.root"  , 115. , channel, do7TeV, 1115 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1116.root"  , 116. , channel, do7TeV, 1116 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1117.root"  , 117. , channel, do7TeV, 1117 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1118.root"  , 118. , channel, do7TeV, 1118 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1119.root"  , 119. , channel, do7TeV, 1119 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1121.root"  , 121. , channel, do7TeV, 1121 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1122.root"  , 122. , channel, do7TeV, 1122 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1123.root"  , 123. , channel, do7TeV, 1123 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1124.root"  , 124. , channel, do7TeV, 1124 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1125.root"  , 125. , channel, do7TeV, 1125 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1126.root"  , 126. , channel, do7TeV, 1126 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1127.root"  , 127. , channel, do7TeV, 1127 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1128.root"  , 128. , channel, do7TeV, 1128 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1129.root"  , 129. , channel, do7TeV, 1129 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1130.root"  , 130. , channel, do7TeV, 1130 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1140.root"  , 140. , channel, do7TeV, 1140 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1150.root"  , 150. , channel, do7TeV, 1150 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1160.root"  , 160. , channel, do7TeV, 1160 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1190.root"  , 190. , channel, do7TeV, 1190 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1200.root"  , 200. , channel, do7TeV, 1200 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1220.root"  , 220. , channel, do7TeV, 1220 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1250.root"  , 250. , channel, do7TeV, 1250 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1300.root"  , 300. , channel, do7TeV, 1300 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1325.root"  , 325. , channel, do7TeV, 1325 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1350.root"  , 350. , channel, do7TeV, 1350 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1400.root"  , 400. , channel, do7TeV, 1400 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1450.root"  , 450. , channel, do7TeV, 1450 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1500.root"  , 500. , channel, do7TeV, 1500 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1550.root"  , 550. , channel, do7TeV, 1550 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1600.root"  , 600. , channel, do7TeV, 1600 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1650.root"  , 650. , channel, do7TeV, 1650 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1700.root"  , 700. , channel, do7TeV, 1700 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1750.root"  , 750. , channel, do7TeV, 1750 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1800.root"  , 800. , channel, do7TeV, 1800 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1850.root"  , 850. , channel, do7TeV, 1850 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1900.root"  , 900. , channel, do7TeV, 1900 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1950.root"  , 950. , channel, do7TeV, 1950 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id11000.root" , 1000., channel, do7TeV, 11000));
        

        Float_t x[18];
        Float_t y[18];
        Float_t xerr[18];
        Float_t yerr[18];
        
        for (int i = 0; i < 18; i++) {
            x[i] = massyieldpairs[i].first;
            y[i] = massyieldpairs[i].second.first;
            xerr[i] = 0.0;
            yerr[i] = massyieldpairs[i].second.second;

            std::cout << x[i] << " " << y[i] << std::endl;
        }
       
        TF1* f72e2mu = new TF1("f72e2mu", "(-4.41814+4.62134*TMath::Erf((x+66.2786)/115.288))*(2.12426+0.00221063*x-1.12005e-06*x*x)", 100, 1000.);
        f72e2mu->SetLineColor(kBlue);
        TF1* f74e    = new TF1("f74e", "(-4.42814+4.61116*TMath::Erf((x+68.4113)/116.562))*(1.8253+0.00230693*x-1.16985e-06*x*x)", 100, 1000.);
        f74e->SetLineColor(kBlue);
        TF1* f74mu    = new TF1("f74mu", "(-4.41479+4.62543*TMath::Erf((x+44.1078)/98.0653))*(2.43353+0.00245858*x-1.16657e-06*x*x)", 100, 1000.);
        f74mu->SetLineColor(kBlue);

        TFile* outfile = new TFile("eff_8TeV_2e2mu.root", "RECREATE");
        TGraphErrors* gr = new TGraphErrors(18, x, y, xerr, yerr);

        gr->Fit("pol3");
    
        gr->Write();           
    
        gr->Draw("AP");

        //f72e2mu->Draw("SAME");
        //f74e->Draw("SAME");
        //f74mu->Draw("SAME");

        TF1* f = (TF1*)gr->GetFunction("pol3");

        std::stringstream ss;
        
        for (int i = 0; i < f->GetNumberFreeParameters(); i++) {
            if (i != 0) ss << " + ";
            ss << "(" << f->GetParameter(i) << ")";
            for (int j = 0; j < i; j++) {
                ss << "*@0";
            }
        }

        std::cout << ss.str() << std::endl;
        outfile->Close();
    }

    else {
        init(do7TeV);

        std::string treeFolder = "/data1/sixie/ntuples/step2/hcp/2011/";

        std::vector<std::pair<float, std::pair<float, float> > > massyieldpairs;

        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1115.root"  , 115. , channel, do7TeV, 1115 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1120.root"  , 120. , channel, do7TeV, 1120 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1130.root"  , 130. , channel, do7TeV, 1130 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1140.root"  , 140. , channel, do7TeV, 1140 ));
        massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1160.root"  , 160. , channel, do7TeV, 1160 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1170.root"  , 170. , channel, do7TeV, 1170 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1180.root"  , 180. , channel, do7TeV, 1180 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1210.root"  , 210. , channel, do7TeV, 1210 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1220.root"  , 220. , channel, do7TeV, 1220 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1230.root"  , 230. , channel, do7TeV, 1230 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1250.root"  , 250. , channel, do7TeV, 1250 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1300.root"  , 300. , channel, do7TeV, 1300 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1325.root"  , 325. , channel, do7TeV, 1325 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1350.root"  , 350. , channel, do7TeV, 1350 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1400.root"  , 400. , channel, do7TeV, 1400 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1450.root"  , 450. , channel, do7TeV, 1450 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1500.root"  , 500. , channel, do7TeV, 1500 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1550.root"  , 550. , channel, do7TeV, 1550 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1600.root"  , 600. , channel, do7TeV, 1600 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1650.root"  , 650. , channel, do7TeV, 1650 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1700.root"  , 700. , channel, do7TeV, 1700 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1750.root"  , 750. , channel, do7TeV, 1750 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1800.root"  , 800. , channel, do7TeV, 1800 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1900.root"  , 900. , channel, do7TeV, 1900 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id1950.root"  , 950. , channel, do7TeV, 1950 ));
//         massyieldpairs.push_back(getYield(treeFolder+"hzzTree_id11000.root" , 1000., channel, do7TeV, 11000));

        Float_t x[5];
        Float_t y[5];
        Float_t xerr[5];
        Float_t yerr[5];

        for (int i = 0; i < 5; i++) {
            x[i] = massyieldpairs[i].first;
            y[i] = massyieldpairs[i].second.first;
            xerr[i] = 0.0;
            yerr[i] = massyieldpairs[i].second.second;

            std::cout << x[i] << " " << y[i] << std::endl;
        }

        TF1* f72e2mu = new TF1("f72e2mu", "(-4.41691+4.62292*TMath::Erf((x+36.5501)/96.1332))*(2.11887+0.00295207*x-1.7409e-06*x*x)", 100, 1000.);
        f72e2mu->SetLineColor(kBlue);
        TF1* f74e    = new TF1("f74e", "(-4.4277+4.6121*TMath::Erf((x+34.5267)/94.973))*(1.95611+0.00282739*x-1.43709e-06*x*x)", 100, 1000.);
        f74e->SetLineColor(kBlue);
        TF1* f74mu    = new TF1("f74mu", "(-4.40824+4.63167*TMath::Erf((x+29.5154)/89.0166))*(2.47423+0.00218329*x-9.59231e-07*x*x)", 100, 1000.);
        f74mu->SetLineColor(kBlue);

        TGraphErrors* gr = new TGraphErrors(5, x, y, xerr, yerr);
        TFile* outfile = new TFile("eff_7TeV_4mu_no.root", "RECREATE");

        gr->Fit("pol3");

        gr->Write();
        
        gr->Draw("AP");

        //f72e2mu->Draw("SAME");
        //f74e->Draw("SAME");
        //f74mu->Draw("SAME");

        TF1* f = (TF1*)gr->GetFunction("pol3");

        std::stringstream ss;

        for (int i = 0; i < f->GetNumberFreeParameters(); i++) {
            if (i != 0) ss << " + ";
            ss << "(" << f->GetParameter(i) << ")";
            for (int j = 0; j < i; j++) {
                ss << "*@0";
            }
        }

        std::cout << ss.str() << std::endl;

        outfile->Close();
    }
    

}   

