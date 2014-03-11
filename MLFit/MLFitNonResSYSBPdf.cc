#include <iostream>
#include <math.h>

#include "MLFit/MLFitNonResSYSBPdf.hh"
#include "RooRealVar.h"
#include "RooRealConstant.h"

ClassImp(MLFitNonResSYSBPdf)

MLFitNonResSYSBPdf::MLFitNonResSYSBPdf(const char *name, const char *title,
					     RooAbsReal& _m, RooAbsReal& _sqrtS, RooAbsReal& _p0, 
					     RooAbsReal& _p1, RooAbsReal& _p2)  
  :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  sqrtS("sqrtS", "sqrtS", this, _sqrtS),
  p0("p0", "p0", this, _p0),
  p1("p1", "p1", this, _p1),
  p2("p2", "p2", this, _p2)
{
}

MLFitNonResSYSBPdf::MLFitNonResSYSBPdf(const MLFitNonResSYSBPdf& other, const char* name) :
  RooAbsPdf(other, name), 
  m("m", this, other.m), 
  sqrtS("sqrtS", this, other.sqrtS),
  p0("p0", this, other.p0),
  p1("p1", this, other.p1),
  p2("p2", this, other.p2)
{
}

Double_t MLFitNonResSYSBPdf::evaluate() const 
{
  double val = pow(1-m/sqrtS+p2*pow(m/sqrtS,2.),p0)/pow(m,p1);

  return val;
}
