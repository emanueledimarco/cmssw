// Code from Pierini

#ifndef RooCruijff_H
#define RooCruijff_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"

using namespace RooFit;

class RooCruijff : public RooAbsPdf {
public:
  RooCruijff() {} ;
  RooCruijff(const char *name, const char *title, RooAbsReal& _m,
        RooAbsReal& _m0, 
        RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
        RooAbsReal& _alphaL, RooAbsReal& _alphaR) ;
  
  RooCruijff(const RooCruijff& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new RooCruijff(*this,newname); }

  inline virtual ~RooCruijff() { }

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;
  RooRealProxy alphaL;
  RooRealProxy alphaR;

  Double_t evaluate() const;

private:
  
  ClassDef(RooCruijff,1)
};



#endif
