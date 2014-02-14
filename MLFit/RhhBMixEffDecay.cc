/*********************************************************************************
 * RhhBMixEffDecay
 * Created Sep 2002 Nick Danielson
 *********************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// The mixing PDF using S and C (no lambda) and using
// reco flavor instead of mixing state.  Uses mu's
// to take into account different B0, B0bar efficiencies.
// Based in large part on RooBMixDecay. 

#include "MLFit/RhhBMixEffDecay.hh"

#include "RooRealIntegral.h"

ClassImp(RhhBMixEffDecay) 
;


RhhBMixEffDecay::RhhBMixEffDecay(const char *name, const char *title, 
				 RooRealVar& t, RooAbsCategory& recoFlav,
				 RooAbsCategory& tagFlav,
				 RooAbsReal& tau, RooAbsReal& dm,			   
				 RooAbsReal& mistag, RooAbsReal& delMistag,
				 RooAbsReal& mu, RooAbsReal& nu,
				 const RooResolutionModel& model, 
				 DecayType type) :
  RooAbsAnaConvPdf(name,title,model,t), 
  _mistag("mistag","Mistag rate",this,mistag),
  _recoFlav("recoFlav","Reco Flavor",this,recoFlav),
  _delMistag("delMistag","Delta mistag rate",this,delMistag),
  _tagFlav("tagFlav","Flavour of tagged B0",this,tagFlav),
  _mu("mu", "Tagging efficiency asymmetry",this,mu),
  _nu("nu", "Reconstruction efficiency asymmetry",this,nu),
  _type(type),
  _tau("tau","Mixing life time",this,tau),
  _dm("dm","Mixing frequency",this,dm),
  _t("_t","time",this,t), _genMixFrac(0)
{
  // Constructor
  switch(type) {
  case SingleSided:
    _basisExp = declareBasis("exp(-@0/@1)",RooArgList(tau,dm)) ;
    _basisCos = declareBasis("exp(-@0/@1)*cos(@0*@2)",RooArgList(tau,dm)) ;
    break ;
  case Flipped:
    _basisExp = declareBasis("exp(@0)/@1)",RooArgList(tau,dm)) ;
    _basisCos = declareBasis("exp(@0/@1)*cos(@0*@2)",RooArgList(tau,dm)) ;
    break ;
  case DoubleSided:
    _basisExp = declareBasis("exp(-abs(@0)/@1)",RooArgList(tau,dm)) ;
    _basisCos = declareBasis("exp(-abs(@0)/@1)*cos(@0*@2)",RooArgList(tau,dm)) ;
    break ;
  }

}


RhhBMixEffDecay::RhhBMixEffDecay(const RhhBMixEffDecay& other, const char* name) : 
  RooAbsAnaConvPdf(other,name),
  _mistag("mistag",this,other._mistag),
  _recoFlav("recoFlav",this,other._recoFlav),
  _delMistag("delMistag",this,other._delMistag),
  _tagFlav("tagFlav",this,other._tagFlav),
  _mu("mu",this,other._mu),
  _nu("nu",this,other._nu),
  _type(other._type),
  _tau("tau",this,other._tau),
  _dm("dm",this,other._dm),
  _t("t",this,other._t),
  _basisExp(other._basisExp),
  _basisCos(other._basisCos),
  _genMixFrac(other._genMixFrac),
  _genFlavFrac(other._genFlavFrac),
  _genFlavFracMix(other._genFlavFracMix),
  _genFlavFracUnmix(other._genFlavFracUnmix)
{
  // Copy constructor
}



RhhBMixEffDecay::~RhhBMixEffDecay()
{
  // Destructor
}


Double_t RhhBMixEffDecay::coefficient(Int_t basisIndex) const 
{
  // Comp with tFit MC: must be (1 - tagFlav*...)
  Double_t D  = 1-2*_mistag;
  Double_t dD = -2*_delMistag;

  if (basisIndex==_basisExp) {
    //return (1 - _tagFlav*_delMistag) ;
    //cout << "exp part:" << (1 + _recoFlav*_nu)*(1 + _tagFlav*dD/2 + _tagFlav*_mu*D) << endl;
    //cout << "while recof=" << 1*_recoFlav << endl;
    //_recoFlav.Print();
    //_recoFlav.arg().Print();
    //_tagFlav.arg().Print();
    
    
    return (1 + _recoFlav*_nu)*(1 + _tagFlav*dD/2 + _tagFlav*_mu*D);
    
  }

  if (basisIndex==_basisCos) {
    return -_recoFlav*(1 + _recoFlav*_nu)*_tagFlav*(D + _tagFlav*_mu*(1 + _tagFlav*dD/2));
  }
  
  return 0 ;
}



Int_t RhhBMixEffDecay::getCoefAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars) const 
{
//   cout << "RhhBMixEffDecay::getCoefAI " ; allVars.Print("1") ;

  if (matchArgs(allVars,analVars,_recoFlav,_tagFlav)) return 3 ;
  if (matchArgs(allVars,analVars,_recoFlav)) return 2 ;
  if (matchArgs(allVars,analVars,_tagFlav)) return 1 ;
  return 0 ;
}



Double_t RhhBMixEffDecay::coefAnalyticalIntegral(Int_t basisIndex, Int_t code, const char* rangeName) const 
{  

  //Double_t chi=1/(1+_dm*_dm*_tau*_tau);

  // BIG HACK: Multiply everything by 1-mu*chi to get agreement with LMinuit!
  //Double_t fudge=1-_recoFlav*_tagFlav*_mu*chi;
  Double_t fudge = 1;
  
  switch(code) {
    // No integration
  case 0: return coefficient(basisIndex) ;

    
    // Integration over 'recoFlav' and 'tagFlav' 
  case 3:
    
   if (basisIndex==_basisExp) {
      return 4.0 ;
      
    }    
    if (basisIndex==_basisCos) {
      //      return 0.0 ;
      return -4*_mu*_nu ;
    }

    // Integration over 'recoFlav'
  case 2:
    
    if (basisIndex==_basisExp) {
      //      return 2.0*coefficient(basisIndex) ;
      Double_t D  = 1-2*_mistag;
      Double_t dD = -2*_delMistag;
      return 2*(1 + _tagFlav*dD/2 + _tagFlav*_mu*D);
    }    
    if (basisIndex==_basisCos) {
      return 0.0 ;
    }

    // Integration over 'tagFlav'
  case 1:
    
    if (basisIndex==_basisExp) {
      //return 2.0/fudge ;
      return (2.0/fudge)*(1 + _recoFlav*_nu);
    }    
    if (basisIndex==_basisCos) {
      //return -2*_recoFlav*_mu/fudge;
      return -(2*_recoFlav*_mu/fudge)*(1 + _recoFlav*_nu);
      
    }
  default:
    assert(0) ;
  }
    
  return 0 ;
}


Int_t RhhBMixEffDecay::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
  if (staticInitOK) {
    if (matchArgs(directVars,generateVars,_t,_recoFlav,_tagFlav)) return 4 ;  
    if (matchArgs(directVars,generateVars,_t,_recoFlav)) return 3 ;  
    if (matchArgs(directVars,generateVars,_t,_tagFlav)) return 2 ;  
  }

  if (matchArgs(directVars,generateVars,_t)) return 1 ;  
  return 0 ;
}



void RhhBMixEffDecay::initGenerator(Int_t code)
{
  switch (code) {
  case 2:
    {
      // Calculate the fraction of B0bar events to generate
      Double_t sumInt = RooRealIntegral("sumInt","sum integral",*this,RooArgSet(_t.arg(),_tagFlav.arg())).getVal() ;
      _tagFlav = 1 ; // B0 
      Double_t flavInt = RooRealIntegral("flavInt","flav integral",*this,RooArgSet(_t.arg())).getVal() ;
      _genFlavFrac = flavInt/sumInt ;
      break ;
    }  
  case 3:
    {
      // Calculate the fraction of mixed events to generate
      Double_t sumInt = RooRealIntegral("sumInt","sum integral",*this,RooArgSet(_t.arg(),_recoFlav.arg())).getVal() ;
      _recoFlav = -1 ; // mixed
      Double_t mixInt = RooRealIntegral("mixInt","mix integral",*this,RooArgSet(_t.arg())).getVal() ;
      _genMixFrac = mixInt/sumInt ;
      break ;
    }  
  case 4:
    {
      // Calculate the fraction of mixed events to generate
      Double_t sumInt = RooRealIntegral("sumInt","sum integral",*this,RooArgSet(_t.arg(),_recoFlav.arg(),_tagFlav.arg())).getVal() ;
      _recoFlav = -1 ; // mixed
      Double_t mixInt = RooRealIntegral("mixInt","mix integral",*this,RooArgSet(_t.arg(),_tagFlav.arg())).getVal() ;
      _genMixFrac = mixInt/sumInt ;
      
      // Calculate the fractio of B0bar tags for mixed and unmixed
      RooRealIntegral dtInt("mixInt","mix integral",*this,RooArgSet(_t.arg())) ;
      _recoFlav = -1 ; // Mixed
      _tagFlav  =  1 ; // B0
      _genFlavFracMix   = dtInt.getVal() / mixInt ;
      _recoFlav =  1 ; // Unmixed
      _tagFlav  =  1 ; // B0
      _genFlavFracUnmix = dtInt.getVal() / (sumInt - mixInt) ;
      break ;
    }
  }
}




void RhhBMixEffDecay::generateEvent(Int_t code)
{
  // Generate mix-state dependent
  switch(code) {
  case 2:
    {
      Double_t rand = RooRandom::uniform() ;
      _tagFlav = (Int_t) ((rand<=_genFlavFrac) ?  1 : -1) ;
      break ;
    }
  case 3:
    {
      Double_t rand = RooRandom::uniform() ;
      _recoFlav = (Int_t) ((rand<=_genMixFrac) ? -1 : 1) ;
      break ;
    }
  case 4:
    {
      Double_t rand = RooRandom::uniform() ;
      _recoFlav = (Int_t) ((rand<=_genMixFrac) ? -1 : 1) ;

      rand = RooRandom::uniform() ;
      Double_t genFlavFrac = (_recoFlav==-1) ? _genFlavFracMix : _genFlavFracUnmix ;
      _tagFlav = (Int_t) ((rand<=genFlavFrac) ?  1 : -1) ;
      break ;
    }
  }

  // Generate delta-t dependent
  while(1) {
    Double_t rand = RooRandom::uniform() ;
    Double_t tval(0) ;

    switch(_type) {
    case SingleSided:
      tval = -_tau*log(rand);
      break ;
    case Flipped:
      tval= +_tau*log(rand);
      break ;
    case DoubleSided:
      tval = (rand<=0.5) ? -_tau*log(2*rand) : +_tau*log(2*(rand-0.5)) ;
      break ;
    }

    // Accept event if T is in generated range
    Double_t maxDil = 1.0 ;
    // 2 in next line is conservative and inefficient - allows for delMistag=1!
    Double_t maxAcceptProb = 2 + maxDil;        
    Double_t acceptProb    = (1 + _recoFlav*_nu)*((1-_tagFlav*_delMistag + _mu*_tagFlav*(1. - 2.*_mistag)) 
                           - (_tagFlav*_recoFlav*(1-2*_mistag) + _recoFlav*_mu*(1. - _tagFlav*_delMistag))*cos(_dm*tval));

    Bool_t accept = maxAcceptProb*RooRandom::uniform() < acceptProb ? kTRUE : kFALSE ;
    
    if (tval<_t.max() && tval>_t.min() && accept) {
      _t = tval ;
      break ;
    }

  }

}

