 /*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitBabar                                                      *
 *    File: $Id: RhhBinnedPdf.cc,v 1.1 2010/06/23 09:42:50 mpierini Exp $
 * Authors:                                                                  *
 *    Aaron Roodman, Stanford Linear Accelerator Center, Stanford University *
 *    Adapted by Wouter                                                      *
 *                                                                           *
 * Copyright (c) 2004, Stanford University. All rights reserved.        *
 *           
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
//

// This is a reimplementation of the
// RooModels/RooParametricStepFunction. In the RooModels implementation
// the coefficient of bin n represents the density in bin n. In this
// implementation, the coefficient for bin n is the integral of the
// pdf in bin 'n', divided by the integral over bins [n,N]. The
// advantage of this approach is that as long as each coefficient is
// limited to values [0,1], there will never be a problem with the
// normalization.
//
// Analytical expressions for calculating each bin integral for a given set of
// parameters:
//    b1 = N p1
//    b2 = N p2 (1-p1)
//    b3 = N p3 (1-p1) (1-p2)
//    b4 = N p4 (1-p1) (1-p2) (1-p3)
//    ... 
//
// An example of usage is:
//
// Int_t nbins(10);
// TArrayD limits(nbins+1);
// limits[0] = 0.0; //etc...
// RooArgList* list = new RooArgList("list");
// RooRealVar* binHeight0 = new RooRealVar("binHeight0","bin 0 Value",0.1,0.0,1.0);
// list->add(binHeight0); // up to binHeight8, ie. 9 parameters
//
// RhhBinnedPdf  aPdf = ("aPdf","PSF",*x,*list,limits);
//

#include <math.h>

#include "MLFit/RhhBinnedPdf.hh"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

ClassImp(RhhBinnedPdf);

RhhBinnedPdf::RhhBinnedPdf(const char* name, const char* title, 
			     RooAbsReal& x, const RooArgList& coefList, const TArrayD& limits) :
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefList","List of coefficients",this),
  _nBins(limits.GetSize()-1)
{
  // Check lowest order
  if (_nBins<0) {
    cout << "RhhBinnedPdf::ctor(" << GetName() 
	 << ") WARNING: nBins must be >=0, setting value to 0" << endl ;
    _nBins=0 ;
  }

  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while(coef = (RooAbsArg*)coefIter->Next()) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RhhBinnedPdf::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter ;

  // Bin limits  
  limits.Copy(_limits);
}


RhhBinnedPdf::RhhBinnedPdf(const RhhBinnedPdf& other, const char* name) :
  RooAbsPdf(other, name), 
  _x("x", this, other._x), 
  _coefList("coefList",this,other._coefList),
  _nBins(other._nBins)
{
  // Copy constructor
  (other._limits).Copy(_limits);
}



RhhBinnedPdf::~RhhBinnedPdf()
{
}


Int_t RhhBinnedPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}



Double_t RhhBinnedPdf::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;

  Double_t sum(1.0) ;
  return sum;  
  
}


Double_t RhhBinnedPdf::evaluate() const 
{
  Double_t xval = _x;
  Double_t value(0);
  if (xval >= _limits[0] && xval < _limits[_nBins]){
    double sum(0),binval(0);
    int i=0 ;
    for (i=0; i<_nBins-1 && xval>=_limits[i];++i) {
      binval = (1-sum)*static_cast<RooRealVar*>(_coefList.at(i))->getVal() ;
      //binval = static_cast<RooRealVar*>(_coefList.at(i))->getVal() ;
      sum    += binval ;
    }
    if( xval>=_limits[_nBins-1] ) { // the last bin
      binval = 1-sum ;
      i = _nBins ;
    }
    double binwidth = _limits[i] - _limits[i-1];
    value = binval/binwidth ;
    if (value<0){
      cout << "RhhBinnedPdf: sum of values gt 1.0 -- beware!!" 
	   << value << " " << binval << " " << sum << " " << i << " " << xval << endl;
      value = 0.000000001;
    }
  }
  return value;
}

Int_t RhhBinnedPdf::getnBins(){
  return _nBins;
}

Double_t* RhhBinnedPdf::getLimits(){
  Double_t* limoutput = _limits.GetArray();
  return limoutput;
}
