#include <iostream>
#include <math.h>

#include "MLFit/MLFitNonResPdf.hh"
#include "RooRealVar.h"
#include "RooRealConstant.h"

ClassImp(MLFitNonResPdf)

MLFitNonResPdf::MLFitNonResPdf(const char *name, const char *title,
			       RooAbsReal& _m, RooAbsReal& _sqrtS, 
			       RooAbsReal& _p0, RooAbsReal& _p1, 
			       RooAbsReal& _p2, RooAbsReal& _p3) 
  :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  sqrtS("sqrtS", "sqrtS", this, _sqrtS),
  p0("p0", "p0", this, _p0),
  p1("p1", "p1", this, _p1),
  p2("p2", "p2", this, _p2),
  p3("p3", "p3", this, _p3)
{
}

MLFitNonResPdf::MLFitNonResPdf(const MLFitNonResPdf& other, const char* name) :
  RooAbsPdf(other, name), 
  m("m", this, other.m), 
  sqrtS("sqrtS", this, other.sqrtS),
  p0("p0", this, other.p0),
  p1("p1", this, other.p1),
  p2("p2", this, other.p2),
  p3("p3", this, other.p3)
{
}

Double_t MLFitNonResPdf::evaluate() const 
{
  double val = p0*pow(1-m/sqrtS,p1)/pow(m/sqrtS,p2+p3*log(m/sqrtS));

  return val;
}
