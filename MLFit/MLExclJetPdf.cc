#include <iostream>
#include <math.h>

#include "MLFit/MLExclJetPdf.hh"
#include "RooRealVar.h"
#include "RooRealConstant.h"

ClassImp(MLExclJetPdf)

MLExclJetPdf::MLExclJetPdf(const char *name, const char *title,
			     RooAbsReal& _m, RooAbsReal& _njets)
  :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  njets("njets", "njets", this, _njets)
{
}

MLExclJetPdf::MLExclJetPdf(const MLExclJetPdf& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), njets("njets", this, other.njets)
{
}

Double_t MLExclJetPdf::evaluate() const 
{
  double val = ( (m-njets<=0.5 && m-njets>-0.5) ? 1. : 0.);
  return val;
}


Int_t MLExclJetPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  if (matchArgs(allVars, analVars,m)) return 1;
  return 0;
}

Double_t MLExclJetPdf::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  Double_t sum(1.0) ;
  return sum;    
}
