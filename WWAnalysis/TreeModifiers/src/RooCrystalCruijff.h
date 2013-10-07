// Code from Pierini


#ifndef RooCrystalCruijff_h
#define RooCrystalCruijff_h

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class RooCrystalCruijff : public RooAbsPdf {
public:
  RooCrystalCruijff(const char *name, const char *title, RooAbsReal& _m,
        RooAbsReal& _m0, 
        RooAbsReal& _sigma, RooAbsReal& _alpha,
        RooAbsReal& _alphaCR, RooAbsReal& _nCR) ;
  
  RooCrystalCruijff(const RooCrystalCruijff& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new RooCrystalCruijff(*this,newname); }

  inline virtual ~RooCrystalCruijff() { }

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alpha;
  RooRealProxy alphaCB;
  RooRealProxy nCB;

  Double_t evaluate() const;

private:
  
  ClassDef(RooCrystalCruijff,1)
};

#endif

