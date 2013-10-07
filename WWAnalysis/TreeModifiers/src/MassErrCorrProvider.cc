#include "WWAnalysis/TreeModifiers/interface/MassErrCorrProvider.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
using namespace std;

void MassErrCorrProvider::init(bool is2011, bool ismc) {
  
  is2011_ = is2011;
  ismc_ = ismc;

  std::string baseFolder(getenv("CMSSW_BASE"));
  std::string filepath = baseFolder + "/src/WWAnalysis/TreeModifiers/data/ebeOverallCorrections.HCP2012.v1.root";
  
  TFile *file = TFile::Open(filepath.c_str());

  stringstream namemu, nameel;
  namemu << "mu_" << (ismc_ ? "mc" : "reco") << (is2011_ ? "42x" : "53x");
  nameel << "el_" << (ismc_ ? "mc" : "reco") << (is2011_ ? "42x" : "53x");

  mapMu_ = new TH2F;
  mapEl_ = new TH2F;

  mapMu_ = (TH2F*)file->Get(namemu.str().c_str())->Clone("mapMu");
  mapEl_ = (TH2F*)file->Get(nameel.str().c_str())->Clone("mapEl");

}


float MassErrCorrProvider::getMassErr(float pt[4], float eta[4], float lid[4], float lmassErr[4], float phoMassErr[2]) {
  float sumErr = pow(phoMassErr[0],2) + pow(phoMassErr[1],2);
  for(int i=1; i<5; ++i) {
    sumErr += pow(lmassErr[i],2);
  }
  return sqrt(sumErr);
}

float MassErrCorrProvider::getMassErrCorr(float pt[4], float eta[4], float lid[4], float lmassErr[4], float phoMassErr[2]) {
  
  float sumErr = pow(phoMassErr[0],2) + pow(phoMassErr[1],2);
  for(int i=1; i<5; ++i) {
    TH2F *map = abs(lid[i]==13) ? mapMu_ : mapEl_;
    int ptBin  = min(max(1,map->GetXaxis()->FindBin(pt[i])), map->GetNbinsX());
    int etaBin  = min(max(1,map->GetYaxis()->FindBin(fabs(eta[i]))), map->GetNbinsY());
    sumErr += pow(lmassErr[i] * map->GetBinContent(ptBin,etaBin),2);
  }
  return sqrt(sumErr);
}
