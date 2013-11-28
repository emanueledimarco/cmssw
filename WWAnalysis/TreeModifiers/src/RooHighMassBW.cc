#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "TMath.h"

#include "RooHighMassBW.h"

ClassImp(RooHighMassBW)

RooHighMassBW::RooHighMassBW(const char *name, const char *title, RooAbsReal& _m, RooAbsReal& _mH, RooAbsReal& _gamma) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  mH("mH", "MH", this, _mH),
  gamma("gamma", "Gamma", this, _gamma)
{
}

RooHighMassBW::RooHighMassBW(const RooHighMassBW& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), mH("mH", this, other.mH), gamma("sigmaL", this, other.gamma)
{
}

Double_t RooHighMassBW::evaluate() const {

  return m / ((m*m - mH*mH)*(m*m - mH*mH) + m*m*gamma*gamma);

}

