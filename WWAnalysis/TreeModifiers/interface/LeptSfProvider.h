#ifndef TreeModifiers_LeptSfProvider_h
#define TreeModifiers_LeptSfProvider_h

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include <string>
#include <iostream>

/** Class cache SF values from 2D histograms and returns them through convenient methods.
 */

class LeptSfProvider {
 public:
 LeptSfProvider() : is2011_(false)  {    }
  virtual ~LeptSfProvider(){
  };

  float getMuonTkSF(float eta);

  float getMuonIDSF(float pt, float eta);

  float getMuonIPSF(float pt, float eta);

  float getMuonIsoSF(float pt, float eta);

  float getMuonSF(float pt, float eta);

  float getElectronSF(float pt, float eta);
  
  float getPR(float pt, float eta, float id);

  float getSF(float pt, float eta, float id);

  void init(bool is2011){
    is2011_ = is2011;
    initMu(is2011);
    initEl(is2011);
  }

  void initMu(bool is2011); 

  void initEl(bool is2011); 
  
 private:
  void setArrayFromGraphY(TGraph* graph, std::vector<double>& vect);


 private:  
  bool is2011_;

  std::vector<double> IP_eta; 
  std::vector<double> IP_pt_barrel; 
  std::vector<double> IP_pt_endcaps; 

  std::vector<double> ISO_eta; 
  std::vector<double> ISO_pt_barrel; 
  std::vector<double> ISO_pt_endcaps; 

  std::vector<double> ID_eta_2011A; 
  std::vector<double> ID_eta_2011B; 
  std::vector<double> ID_eta; 
  std::vector<double> ID_pt_barrel; 
  std::vector<double> ID_pt_endcaps; 


  std::vector<float> mutkscalefactors;    
  std::vector<std::vector<float> > elesf;


};

#endif
