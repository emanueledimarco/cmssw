#ifndef RooHighMassBW_H
#define RooHighMassBW_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"

using namespace RooFit;

class RooHighMassBW : public RooAbsPdf {
public:
  RooHighMassBW() {} ;
  RooHighMassBW(const char *name, const char *title, RooAbsReal& _m,
        RooAbsReal& _m0, RooAbsReal& _gamma) ; 
  RooHighMassBW(const RooHighMassBW& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new RooHighMassBW(*this,newname); }

  inline virtual ~RooHighMassBW() { }

protected:

  RooRealProxy m;
  RooRealProxy mH;
  RooRealProxy gamma;

  Double_t evaluate() const;

private:
  
  ClassDef(RooHighMassBW,1)
};



#endif
