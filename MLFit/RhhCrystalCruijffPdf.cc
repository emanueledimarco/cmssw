#include <iostream>
#include <math.h>

#include "MLFit/RhhCrystalCruijffPdf.hh"
#include "RooRealVar.h"
#include "RooRealConstant.h"

ClassImp(RhhCrystalCruijffPdf)

RhhCrystalCruijffPdf::RhhCrystalCruijffPdf(const char *name, const char *title,
					   RooAbsReal& _m, RooAbsReal& _m0, 
					   RooAbsReal& _sigma, RooAbsReal& _alpha, 
					   RooAbsReal& _alphaCB, RooAbsReal& _nCB)
  :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "SigmaL", this, _sigma),
  alpha("alpha", "Alpha", this, _alpha),
  alphaCB("alphaCB", "AlphaCB", this, _alphaCB),
  nCB("nCB", "nCB", this, _nCB)
{
}

RhhCrystalCruijffPdf::RhhCrystalCruijffPdf(const RhhCrystalCruijffPdf& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), alpha("alpha", this, other.alpha), 
  alphaCB("alphaCB", this, other.alphaCB), nCB("nCB", this, other.nCB) 
{
}

Double_t RhhCrystalCruijffPdf::evaluate() const 
{
  double dx = (m-m0) ;
  double val;
  if(dx<0) { // Cruijff 
    double f = 2*sigma*sigma + alpha*dx*dx ;
    val = exp(-dx*dx/f);
  } else { // CryBall
    if(dx/sigma < alphaCB) {
      val = exp(-dx*dx/2./sigma/sigma);
    } else {
      double A = pow(nCB/fabs(alphaCB),nCB)*exp(-pow(alphaCB,2.)/2.);
      double B = nCB/fabs(alphaCB)-fabs(alphaCB);
      val = A*pow(B+dx/sigma,-nCB);
    }
  }
  return val;
}
