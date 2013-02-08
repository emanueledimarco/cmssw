#ifndef QGLikelihoodCalculator_h
#define QGLikelihoodCalculator_h

#include <string>

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"

#include <map>
#include <vector>

class QGLikelihoodCalculator{
 public:
   QGLikelihoodCalculator(const TString dataDir, Bool_t chs = false);
   ~QGLikelihoodCalculator();
   Float_t QGvalue(std::map<TString, Float_t>);

 private:
   std::map<TString,JetCorrectorParameters*> JCP;
   std::map<TString,SimpleJetCorrector*> SJC;
   std::vector<TString> names;
};

#endif
