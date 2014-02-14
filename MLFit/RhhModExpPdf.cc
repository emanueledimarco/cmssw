#include <iostream>
#include <math.h>

#include "MLFit/RhhModExpPdf.hh"
#include "RooRealVar.h"
#include "RooRealConstant.h"

ClassImp(RhhModExpPdf)

RhhModExpPdf::RhhModExpPdf(const char *name, const char *title,RooAbsReal& _m,
	     RooAbsReal& _slope, RooAbsReal& _shape)
  :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  slope("slope", "slope", this, _slope),
  shape("shape", "shape", this, _shape)
{
}

RhhModExpPdf::RhhModExpPdf(const RhhModExpPdf& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), slope("slope", this, other.slope),
  shape("shape", this, other.shape)
{
}

Double_t RhhModExpPdf::evaluate() const 
{
  return exp( (m-m.min())/(slope + slope*shape*(m-m.min())) ) ;
}
