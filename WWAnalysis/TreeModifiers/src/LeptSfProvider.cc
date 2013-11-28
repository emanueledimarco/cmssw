#include "WWAnalysis/TreeModifiers/interface/LeptSfProvider.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include <cmath>
#include <cstdlib>
using namespace std;


float LeptSfProvider::getMuonTkSF(float eta) {

    int bin = -1;

    if (eta>-2.4 && eta<=-2.1) bin = 0;
    if (eta>-2.1 && eta<=-1.6) bin = 1;
    if (eta>-1.6 && eta<=-1.2) bin = 2;
    if (eta>-1.2 && eta<=-0.9) bin = 3;
    if (eta>-0.9 && eta<=-0.6) bin = 4;
    if (eta>-0.6 && eta<=-0.3) bin = 5;
    if (eta>-0.3 && eta<=-0.2) bin = 6;
    if (eta>-0.2 && eta<= 0.2) bin = 7;
    if (eta> 0.2 && eta<= 0.3) bin = 8;
    if (eta> 0.3 && eta<= 0.6) bin = 9;
    if (eta> 0.6 && eta<= 0.9) bin = 10;
    if (eta> 0.9 && eta<= 1.2) bin = 11;
    if (eta> 1.2 && eta<= 1.6) bin = 12;
    if (eta> 1.6 && eta<= 2.1) bin = 13;
    if (eta> 2.1 && eta<= 2.4) bin = 14;

    if (bin >= 0) return mutkscalefactors[bin];

    return 0.0;

}


float LeptSfProvider::getMuonIDSF(float pt, float eta) {

    int bin = -1;

    if (pt > 20) {

        if (eta>-2.4 && eta<=-2.1) bin = 0;
        if (eta>-2.1 && eta<=-1.6) bin = 1;
        if (eta>-1.6 && eta<=-1.2) bin = 2;
        if (eta>-1.2 && eta<=-0.9) bin = 3;
        if (eta>-0.9 && eta<=-0.6) bin = 4;
        if (eta>-0.6 && eta<=-0.3) bin = 5;
        if (eta>-0.3 && eta<=-0.2) bin = 6;
        if (eta>-0.2 && eta<= 0.2) bin = 7;
        if (eta> 0.2 && eta<= 0.3) bin = 8;
        if (eta> 0.3 && eta<= 0.6) bin = 9;
        if (eta> 0.6 && eta<= 0.9) bin = 10;
        if (eta> 0.9 && eta<= 1.2) bin = 11;
        if (eta> 1.2 && eta<= 1.6) bin = 12;
        if (eta> 1.6 && eta<= 2.1) bin = 13;
        if (eta> 2.1 && eta<= 2.4) bin = 14;

        if (bin >= 0) {
            if (is2011_) return 0.42*ID_eta_2011A[bin] + 0.58*ID_eta_2011B[bin];
            else return ID_eta[bin];
        }

    }

    else {
                
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_barrel[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_endcaps[bin];  
        }

    }

    return 0.0;
}

float LeptSfProvider::getMuonIPSF(float pt, float eta) {

    int bin = -1;

    if (pt > 20) {

        if (eta>-2.4 && eta<=-2.1) bin = 0;
        if (eta>-2.1 && eta<=-1.6) bin = 1;
        if (eta>-1.6 && eta<=-1.2) bin = 2;
        if (eta>-1.2 && eta<=-0.9) bin = 3;
        if (eta>-0.9 && eta<=-0.6) bin = 4;
        if (eta>-0.6 && eta<=-0.3) bin = 5;
        if (eta>-0.3 && eta<=-0.2) bin = 6;
        if (eta>-0.2 && eta<= 0.2) bin = 7;
        if (eta> 0.2 && eta<= 0.3) bin = 8;
        if (eta> 0.3 && eta<= 0.6) bin = 9;
        if (eta> 0.6 && eta<= 0.9) bin = 10;
        if (eta> 0.9 && eta<= 1.2) bin = 11;
        if (eta> 1.2 && eta<= 1.6) bin = 12;
        if (eta> 1.6 && eta<= 2.1) bin = 13;
        if (eta> 2.1 && eta<= 2.4) bin = 14;

        if (bin >= 0) return IP_eta[bin];

    }

    else {
                     
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_barrel[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_endcaps[bin];
        }
        
    }

    return 0.0;

}

float LeptSfProvider::getMuonIsoSF(float pt, float eta) {

    int bin = -1;

    if (pt > 20) {

        if (eta>-2.4 && eta<=-2.1) bin = 0;
        if (eta>-2.1 && eta<=-1.6) bin = 1;
        if (eta>-1.6 && eta<=-1.2) bin = 2;
        if (eta>-1.2 && eta<=-0.9) bin = 3;
        if (eta>-0.9 && eta<=-0.6) bin = 4;
        if (eta>-0.6 && eta<=-0.3) bin = 5;
        if (eta>-0.3 && eta<=-0.2) bin = 6;
        if (eta>-0.2 && eta<= 0.2) bin = 7;
        if (eta> 0.2 && eta<= 0.3) bin = 8;
        if (eta> 0.3 && eta<= 0.6) bin = 9;
        if (eta> 0.6 && eta<= 0.9) bin = 10;
        if (eta> 0.9 && eta<= 1.2) bin = 11;
        if (eta> 1.2 && eta<= 1.6) bin = 12;
        if (eta> 1.6 && eta<= 2.1) bin = 13;
        if (eta> 2.1 && eta<= 2.4) bin = 14;

        if (bin >= 0) return ISO_eta[bin];

    }

    else {

        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_barrel[bin];                               
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_endcaps[bin];
        }

    }

    return 0.0;
}

float LeptSfProvider::getMuonSF(float pt, float eta) {
    if (is2011_) return getMuonIDSF(pt, eta) * getMuonIsoSF(pt, eta) * getMuonIPSF(pt, eta);
    else         return getMuonIDSF(pt, eta) * getMuonIsoSF(pt, eta) * getMuonIPSF(pt, eta) * getMuonTkSF(eta);
}

float LeptSfProvider::getElectronSF(float pt, float eta) {
    int etabin = -1;
    int ptbin = -1;

    float abseta = fabs(eta);

    if (pt> 7 && pt<=10) ptbin = 0;
    if (pt>10 && pt<=15) ptbin = 1;
    if (pt>15 && pt<=20) ptbin = 2;
    if (pt>20 && pt<=30) ptbin = 3;
    if (pt>30 && pt<=40) ptbin = 4;
    if (pt>40 && pt<=50) ptbin = 5;
    if (pt>50 && pt<=80) ptbin = 6;
    if (pt>80)           ptbin = 6;

    if (abseta<=0.8)                     etabin = 0;
    if (abseta>0.8    && abseta<=1.4442) etabin = 1;
    if (abseta>1.4442 && abseta<=1.566)  etabin = 2;
    if (abseta>1.566  && abseta<=2.0)    etabin = 3;
    if (abseta>2.0    && abseta<=2.5)    etabin = 4;

    if (ptbin>=0 && etabin>=0) return elesf[ptbin][etabin];
    else return 0.0;
}

float LeptSfProvider::getPR(float pt, float eta, float id) {

    if (id == 13 || id == -13) return getMuonIDSF(pt, eta) * getMuonIsoSF(pt, eta);
    else return getElectronSF(pt, eta);

}

float LeptSfProvider::getSF(float pt, float eta, float id) {

    if (id == 13 || id == -13) return getMuonSF(pt, eta);

    else return getElectronSF(pt, eta);

}



void LeptSfProvider::initMu(bool is2011) {

    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + "/src/WWAnalysis/TreeModifiers/data/mu_sf.root";

    TFile file(filepath.c_str());

    if (is2011) {
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("IP_eta_2011")), IP_eta );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("IP_pt_barrel_2011")), IP_pt_barrel );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("IP_pt_endcaps_2011")), IP_pt_endcaps );
        
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ISO_eta_2011")), ISO_eta );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ISO_pt_barrel_2011")), ISO_pt_barrel );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ISO_pt_endcaps_2011")), ISO_pt_endcaps );
      
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ID_eta_2011A")), ID_eta_2011A );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ID_eta_2011B")), ID_eta_2011B );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ID_pt_barrel_2011")), ID_pt_barrel );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ID_pt_endcaps_2011")), ID_pt_endcaps );
    }

    else {
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("IP_eta_2012")), IP_eta );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("IP_pt_barrel_2012")), IP_pt_barrel );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("IP_pt_endcaps_2012")), IP_pt_endcaps );
      
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ISO_eta_2012")), ISO_eta );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ISO_pt_barrel_2012")), ISO_pt_barrel );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ISO_pt_endcaps_2012")), ISO_pt_endcaps );
      
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ID_eta_2012")), ID_eta );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ID_pt_barrel_2012")), ID_pt_barrel );
      setArrayFromGraphY( dynamic_cast<TGraphAsymmErrors*>(file.Get("ID_pt_endcaps_2012")), ID_pt_endcaps );

    }
        
    TH1* mutkscalefactorshist = (TH1*)file.Get("TH1D_TK_2012");    

    for (int i = 0; i < mutkscalefactorshist->GetNbinsX(); i++) {
        mutkscalefactors.push_back( mutkscalefactorshist->GetBinContent(i+1));
    }
}

void LeptSfProvider::initEl(bool is2011) {
    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + (is2011 ? "/src/WWAnalysis/TreeModifiers/data/el_sf_2011.root" : "/src/WWAnalysis/TreeModifiers/data/el_sf_2012.root");

    TFile file(filepath.c_str());


    TH2* electronscalefactors = (TH2*)file.Get("heff_electron_selection");

    if (elesf.size() == 0) {
        for (std::size_t i = 0; i < 7; i++) {
            elesf.push_back(std::vector<float>(5, 0.0));
        }
    }
    for (std::size_t i = 0; i < 7; i++) {
        for (std::size_t j = 0; j < 5; j++) {
            elesf[i][j] = electronscalefactors->GetBinContent(i+1, j+1);
        }
    }
}




void LeptSfProvider::setArrayFromGraphY(TGraph* graph,std::vector<double>& vect){
  double* values = graph->GetY();
  vect.clear();
  for(int i=0; i< graph->GetN(); ++i){
    vect.push_back(values[i]);
  }
}
