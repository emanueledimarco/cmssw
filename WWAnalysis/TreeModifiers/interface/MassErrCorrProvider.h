#ifndef TreeModifiers_MassErrCorrProvider_h
#define TreeModifiers_MassErrCorrProvider_h

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include <string>
#include <iostream>

/** Class cache lepton momentum error corrections values from 2D histograms and returns them through convenient methods.
 */

class MassErrCorrProvider {
 public:
 MassErrCorrProvider(): 
   is2011_(false),
   ismc_(false)
     { 
       mapMu_ = 0;
       mapEl_ = 0;
     }
  virtual ~MassErrCorrProvider(){
  };

  float getMassErr(float pt[4], float eta[4], float lid[4], float lmassErr[4], float phoMassErr[2]);

  float getMassErrCorr(float pt[4], float eta[4], float lid[4], float lmassErr[4], float phoMassErr[2]);

  void init(bool is2011, bool ismc);

 private:  

  bool is2011_;
  bool ismc_;

  TH2F *mapMu_, *mapEl_;

};

#endif
