/*******************************************************
 * MLPdfFactory
 * Nick Danielson / Amir Farbin
 *******************************************************/

//#include "TH1F.h"

#include "TClass.h"
#include "RooAddModel.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooResolutionModel.h"
#include "RooStringVar.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooArgusBG.h"
#include "RooBCPGenDecay.h"
#include "RooBifurGauss.h"
#include "RooBMixDecay.h"
#include "RooDecay.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooCBShape.h"
#include "RooLandau.h"
#include "RooPolynomial.h"
#include "RooParametricStepFunction.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "MLFit/RhhBMixEffDecay.hh"
#include "MLFit/RhhCruijffPdf.hh"
#include "MLFit/MLPdfFactory.hh"
#include "MLFit/RhhCruijffPdf.hh"
#include "MLFit/RhhCrystalCruijffPdf.hh"
#include "MLFit/MLExclJetPdf.hh"
#include "MLFit/RhhBinnedPdf.hh"
#include "MLFit/RhhModExpPdf.hh"
#include "MLFit/VecbosBtagPdf.hh"
#include "Roo2DKeysPdf.h"
#include "RooKeysPdf.h"

#include "MLFit/MLFitNonResPdf.hh"
#include "MLFit/MLFitNonResSYSAPdf.hh"
#include "MLFit/MLFitNonResSYSBPdf.hh"
#include "MLFit/MLFitNonResSYSCPdf.hh"

ClassImp(MLPdfFactory);

RooArgList MLPdfFactory::makeParameters(TString pdftype, TString basename, const TList &args)
{
  if (pdftype == "Argus")              return makeArgusParams(basename);
  if (pdftype == "BCPGenDecay")        return makeBCPGenDecayParams(basename);
  if (pdftype == "B0BarCPGenDecay")    return makeBCPGenDecayParams(basename);
  if (pdftype == "BMixDecay")          return makeBMixDecayParams(basename);
  if (pdftype == "BMixEffDecay")       return makeBMixEffDecayParams(basename);
  if (pdftype == "BMixEffDecayMistag") return makeBMixEffDecayMistagParams(basename);
  if (pdftype == "B0MixEffDecay")      return makeBMixEffDecayParams(basename);
  if (pdftype == "B0BarMixEffDecay")   return makeBMixEffDecayParams(basename);
  if (pdftype == "BifurGauss")         return makeBifurGaussParams(basename);
  if (pdftype == "DoubleGaussian")     return makeDoubleGaussianParams(basename);
  if (pdftype == "TripleGaussian")     return makeTripleGaussianParams(basename);
  if (pdftype == "DoubleLandau")       return makeDoubleLandauParams(basename);
  if (pdftype == "Pulse")              return makePulseParams(basename);
  if (pdftype == "CryBall")            return makeCryBallParams(basename);
  if (pdftype == "Voigtian")           return makeVoigtianParams(basename);
  if (pdftype == "DoubleVoigtian")     return makeDoubleVoigtianParams(basename);
  if (pdftype == "EffDecay3G")         return makeEffDecay3GParams(basename);
  if (pdftype == "EffDecay")           return makeEffDecayParams(basename);
  if (pdftype == "Gaussian")           return makeGaussianParams(basename);
  if (pdftype == "hhDeGaussian")       return makeHhDeGaussianParams(basename);
  if (pdftype == "hhDe2Gaussian")      return makeHhDe2GaussianParams(basename);
  if (pdftype == "hhDe3Gaussian")      return makeHhDe3GaussianParams(basename);
  if (pdftype == "hhDeCryBall")        return makeHhDeCryBallParams(basename);
  if (pdftype == "hhDeCryBallG")       return makeHhDeCryBallGParams(basename);
  if (pdftype == "hhDeCruijff")        return makeHhDeCruijffParams(basename);
  if (pdftype == "DeDoubleGaussian")   return makeDeDoubleGaussianParams(basename);
  if (pdftype == "hh0DeCryBall")       return makeHh0DeCryBallParams(basename);
  if (pdftype == "InvArgus")           return makeArgusParams(basename);
  if (pdftype == "Landau")             return makeLandauParams(basename);
  if (pdftype == "Line")               return makeLineParams(basename);
  if (pdftype == "Poly2")              return makePoly2Params(basename);
  if (pdftype == "Totti")              return makeTottiParams(basename);
  if (pdftype == "Cruijff")            return makeCruijffParams(basename);
  if (pdftype == "DoubleCruijff")      return makeDoubleCruijffParams(basename);
  if (pdftype == "CrystalCruijff")     return makeCrystalCruijffParams(basename);
  if (pdftype == "BadSigDTPdf")        return makeBadSigDTParams(basename);
  if (pdftype == "BadBkgDTPdf")        return makeBadBkgDTParams(basename) ;
  if (pdftype == "NoPdf")              return RooArgList();
  if (pdftype == "StepFunction")       return makeStepFunctionParams(basename, args);
  if (pdftype == "BinnedPdf")          return makeBinnedPdfParams(basename, args);
  if (pdftype == "BreitWigner")        return makeBreitWignerParams(basename);
  if (pdftype == "BreitWignerAndPoly") return makeBreitWignerAndPolyParams(basename);
  if (pdftype == "BreitWignerAndExpo") return makeBreitWignerAndExpoParams(basename);
  if (pdftype == "Keys")               return makeKeysParams(basename, args);
  if (pdftype == "2DKeys")             return make2DKeysParams(basename, args);
  if (pdftype == "Hist")               return makeHistParams(basename, args);
  if (pdftype == "Expo")               return makeExpoParams(basename, args);
  if (pdftype == "ModExp")             return makeModExpParams(basename, args);
  if (pdftype == "ExclJet")            return makeExclJetParams(basename);
  if (pdftype == "BtagPdf")            return makeBtagParams(basename);

  if (pdftype == "NonRes")             return makeNonResParams(basename);
  if (pdftype == "NonResSYSA")         return makeNonResParamsSYSA(basename);
  if (pdftype == "NonResSYSB")         return makeNonResParamsSYSB(basename);
  if (pdftype == "NonResSYSC")         return makeNonResParamsSYSC(basename);

  // RFs...
  if (pdftype == "Scaled2G")          return makeScaled2GParams(basename);
  if (pdftype == "Scaled3G")          return makeScaled3GParams(basename);
 
  

  // If we got here, pdftype ain't in the list!
  cout << "MLPdfFactory couldn't find pdf named: " << pdftype << endl;
  return RooArgList(); 
}

RooAbsPdf* MLPdfFactory::makePdf(TString pdftype, TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  
  args.Print();

  if (pdftype == "Argus")              return makeArgusPdf(pdfname, obs, parameters, args);
  if (pdftype == "BifurGauss")         return makeBifurGaussPdf(pdfname, obs, parameters, args);
  if (pdftype == "DoubleGaussian")     return makeDoubleGaussianPdf(pdfname, obs, parameters, args);
  if (pdftype == "TripleGaussian")     return makeTripleGaussianPdf(pdfname, obs, parameters, args);
  if (pdftype == "DoubleLandau")       return makeDoubleLandauPdf(pdfname, obs, parameters, args);
  if (pdftype == "Pulse")              return makePulsePdf(pdfname, obs, parameters, args);
  if (pdftype == "CryBall")            return makeCryBallPdf(pdfname, obs, parameters, args);
  if (pdftype == "Voigtian")           return makeVoigtianPdf(pdfname, obs, parameters, args);
  if (pdftype == "DoubleVoigtian")     return makeDoubleVoigtianPdf(pdfname, obs, parameters, args);
  if (pdftype == "Gaussian")           return makeGaussianPdf(pdfname, obs, parameters, args);
  if (pdftype == "hhDeGaussian")       return makeHhDeGaussianPdf(pdfname, obs, parameters, args);
  if (pdftype == "hhDe2Gaussian")      return makeHhDe2GaussianPdf(pdfname, obs, parameters, args);
  if (pdftype == "hhDe3Gaussian")      return makeHhDe3GaussianPdf(pdfname, obs, parameters, args);
  if (pdftype == "hhDeCryBall")        return makeHhDeCryBallPdf(pdfname, obs, parameters, args);
  if (pdftype == "hhDeCryBallG")       return makeHhDeCryBallGPdf(pdfname, obs, parameters, args);
  if (pdftype == "hhDeCruijff")        return makeHhDeCruijffPdf(pdfname, obs, parameters, args);
  if (pdftype == "DeDoubleGaussian")   return makeDeDoubleGaussianPdf(pdfname, obs, parameters, args);
  if (pdftype == "hh0DeCryBall")       return makeHh0DeCryBallPdf(pdfname, obs, parameters, args);
  if (pdftype == "InvArgus")           return makeInvArgusPdf(pdfname, obs, parameters, args);
  if (pdftype == "Landau")             return makeLandauPdf(pdfname, obs, parameters, args);
  if (pdftype == "Line")               return makeLinePdf(pdfname, obs, parameters, args);
  if (pdftype == "Poly2")              return makePoly2Pdf(pdfname, obs, parameters, args);
  if (pdftype == "Cruijff")            return makeCruijffPdf(pdfname, obs, parameters, args);
  if (pdftype == "DoubleCruijff")      return makeDoubleCruijffPdf(pdfname, obs, parameters, args);
  if (pdftype == "CrystalCruijff")     return makeCrystalCruijffPdf(pdfname, obs, parameters, args);
  if (pdftype == "BadSigDTPdf")        return makeBadSigDTPdf(pdfname, obs, parameters, args);
  if (pdftype == "BadBkgDTPdf")        return makeBadBkgDTPdf(pdfname,obs,parameters,args);
  if (pdftype == "NoPdf")              return makeNoPdf(pdfname, obs, parameters, args);
  if (pdftype == "Prebuilt")           return makePrebuiltPdf(args);
  if (pdftype == "Totti")              return makeTottiPdf(pdfname, obs, parameters, args);
  if (pdftype == "Cruijff")            return makeCruijffPdf(pdfname, obs, parameters, args);
  if (pdftype == "StepFunction")       return makeStepFunctionPdf(pdfname, obs, parameters, args);
  if (pdftype == "BinnedPdf")          return makeBinnedPdf(pdfname, obs, parameters, args);
  if (pdftype == "BreitWigner")        return makeBreitWignerPdf(pdfname, obs, parameters, args);
  if (pdftype == "BreitWignerAndPoly") return makeBreitWignerAndPolyPdf(pdfname, obs, parameters, args);
  if (pdftype == "Keys")               return makeKeysPdf(pdfname, obs, parameters, args);
  if (pdftype == "BreitWignerAndExpo") return makeBreitWignerAndExpoPdf(pdfname, obs, parameters, args);
  if (pdftype == "2DKeys")             return make2DKeysPdf(pdfname, obs, parameters, args);
  if (pdftype == "Hist")               return makeHistPdf(pdfname, obs, parameters, args);
  if (pdftype == "Expo")               return makeExpoPdf(pdfname, obs, parameters, args);
  if (pdftype == "ModExp")             return makeModExpPdf(pdfname, obs, parameters, args);
  if (pdftype == "ExclJet")            return makeExclJetPdf(pdfname, obs, parameters, args);
  if (pdftype == "BtagPdf")            return makeBtagPdf(pdfname, obs, parameters, args);

  if (pdftype == "NonRes")          return makeNonResPdf(pdfname, obs, parameters, args);
  if (pdftype == "NonResSYSA")      return makeNonResPdfSYSA(pdfname, obs, parameters, args);
  if (pdftype == "NonResSYSB")      return makeNonResPdfSYSB(pdfname, obs, parameters, args);
  if (pdftype == "NonResSYSC")      return makeNonResPdfSYSC(pdfname, obs, parameters, args);

  // If we got here, pdftype ain't in the list!
  cout << "MLPdfFactory couldn't find pdf of type: " << pdftype << endl;
  return 0;
}

RooResolutionModel* MLPdfFactory::makeRF(TString pdftype, TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (pdftype == "Scaled2G")   return makeScaled2GRF(pdfname, obs, parameters, args);
  if (pdftype == "Scaled3G")   return makeScaled3GRF(pdfname, obs, parameters, args);

  cout << "MLPdfFactory couldn't find RF of type: " << pdftype << endl;
  return 0;
}


RooAbsPdf* MLPdfFactory::makeConvolutedPdf(TString pdftype, TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  
  if (pdftype == "BCPGenDecay")        return makeBCPGenDecayPdf(pdfname, obs, parameters, rf, args);
  if (pdftype == "B0BarCPGenDecay")    return makeB0BarCPGenDecayPdf(pdfname, obs, parameters, rf, args);
  if (pdftype == "BMixDecay")          return makeBMixDecayPdf(pdfname, obs, parameters, rf, args);
  if (pdftype == "BMixEffDecay")       return makeBMixEffDecayPdf(pdfname, obs, parameters, rf, args);
  if (pdftype == "BMixEffDecayMistag") return makeBMixEffDecayMistagPdf(pdfname, obs, parameters, rf, args);
  if (pdftype == "B0MixEffDecay")      return makeB0MixEffDecayPdf(pdfname, obs, parameters, rf, args);
  if (pdftype == "B0BarMixEffDecay")   return makeB0BarMixEffDecayPdf(pdfname, obs, parameters, rf, args);
  if (pdftype == "EffDecay")           return makeEffDecayPdf(pdfname, obs, parameters, rf, args);
  
  // If we got here, pdftype ain't in the list!
  cout << "MLPdfFactory couldn't find convoluted pdf of type: " << pdftype << endl;
  return 0;
}


Bool_t MLPdfFactory::checkArgLength(TString pdfname, const RooArgList &obs, Int_t nobs, RooArgList parameters, Int_t npar, const TList &args, Int_t narg) 
{
  // This takes the observable list, parameter list and argument list and makes sure they have the right lengths...
  if (obs.getSize() != nobs || parameters.getSize() != npar || args.GetSize() != narg) {
    cout << "MLPdfFactory Error in " << pdfname << ": wrong number of arguments." << endl;
    cout << "Expected " << nobs << " observables, got:" << endl;
    obs.Print();
    cout << "Expected " << npar << " parameters, got:" << endl;
    parameters.Print();
    cout << "Expected " << narg << " arguments, got:" << endl;
    args.Print();
    return kFALSE;
  }
  return kTRUE;
  
}


//===============================
// PDFs
//===============================


RooArgList MLPdfFactory::makeArgusParams(TString basename)
{
  RooRealVar *cutoff = new RooRealVar(basename+"_cutoff",basename+"_cutoff",5.295);
  RooRealVar *shape  = new RooRealVar(basename+"_shape",basename+"_shape",-20.524);

  return RooArgList(*cutoff, *shape);
}

RooAbsPdf* MLPdfFactory::makeArgusPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  RooAbsPdf* rc(0) ;
  if( checkArgLength(pdfname, obs, 1, parameters, 2, args, 0)) {
    RooRealVar *half = new RooRealVar("half", "half", 0.5);
    rc = new RooArgusBG(pdfname, pdfname, *(RooRealVar*)obs.at(0), 
			*(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1), *half);
  } else if( checkArgLength(pdfname, obs, 2, parameters, 2, args, 0)) {
    RooRealVar *half = new RooRealVar("half", "half", 0.5);
    rc = new RooArgusBG(pdfname, pdfname, *(RooRealVar*)obs.at(0), 
			*(RooRealVar*)obs.at(1), *(RooAbsReal*)parameters.at(1), *half);
  }
  return rc ;
}

RooArgList MLPdfFactory::makeBCPGenDecayParams(TString basename) 
{
  RooRealVar *tauB0 = new RooRealVar(basename+"_tauB0",basename+"_tauB0",1);
  RooRealVar *dm = new RooRealVar(basename+"_dm",basename+"_dm",1);
  RooRealVar *D = new RooRealVar(basename+"_D",basename+"_D",1);
  RooRealVar *dD = new RooRealVar(basename+"_dD",basename+"_dD",0);
  RooRealVar *S  = new RooRealVar(basename+"_S",basename+"_S",0);
  RooRealVar *C  = new RooRealVar(basename+"_C",basename+"_C",0);
  RooRealVar *mu = new RooRealVar(basename+"_mu",basename+"_mu",0);
  
  return RooArgList(*tauB0,*dm,*D,*dD,*S,*C,*mu);
}

RooAbsPdf* MLPdfFactory::makeBCPGenDecayPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 7, args, 0)) return 0;

  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(1);

  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *dm        = (RooRealVar*)parameters.at(1);
  RooRealVar *D         = (RooRealVar*)parameters.at(2);
  RooRealVar *dD        = (RooRealVar*)parameters.at(3);
  RooRealVar *S         = (RooRealVar*)parameters.at(4);
  RooRealVar *C         = (RooRealVar*)parameters.at(5);
  RooRealVar *mu        = (RooRealVar*)parameters.at(6);
  RooFormulaVar *w  = new RooFormulaVar(pdfname+"_w",pdfname+"_w","(1-@0)/2",*D);
  RooFormulaVar *dw = new RooFormulaVar(pdfname+"_dw",pdfname+"_dw","-@0/2",*dD);

  RooBCPGenDecay *pdf = new RooBCPGenDecay(pdfname,pdfname,*dt,*tag,*tau,*dm,*w,*C,*S,*dw,*mu,rf,RooBCPGenDecay::DoubleSided);
  return pdf;
}

RooAbsPdf* MLPdfFactory::makeB0BarCPGenDecayPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 7, args, 0)) return 0;

  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(1);

  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *dm        = (RooRealVar*)parameters.at(1);
  RooRealVar *D         = (RooRealVar*)parameters.at(2);
  RooRealVar *dD        = (RooRealVar*)parameters.at(3);
  RooRealVar *S         = (RooRealVar*)parameters.at(4);
  RooRealVar *C         = (RooRealVar*)parameters.at(5);
  RooRealVar *mu        = (RooRealVar*)parameters.at(6);
  RooFormulaVar *minusC = new RooFormulaVar(pdfname+"_-C",pdfname+"_-C","-@0",*C);
  RooFormulaVar *w  = new RooFormulaVar(pdfname+"_w",pdfname+"_w","(1-@0)/2",*D);
  RooFormulaVar *dw = new RooFormulaVar(pdfname+"_dw",pdfname+"_dw","-@0/2",*dD);

  RooBCPGenDecay *pdf = new RooBCPGenDecay(pdfname,pdfname,*dt,*tag,*tau,*dm,*w,*minusC,*S,*dw,*mu,rf,RooBCPGenDecay::DoubleSided);
  return pdf;
}

RooArgList MLPdfFactory::makeBMixDecayParams(TString basename)
{
  RooRealVar *tau = new RooRealVar(basename+"_tau",basename+"_tau",1);
  RooRealVar *dm = new RooRealVar(basename+"_dm",basename+"_dm",0);
  RooRealVar *D = new RooRealVar(basename+"_D",basename+"_D",0);
  RooRealVar *dD = new RooRealVar(basename+"_dD",basename+"_dD",0);
  RooArgList params(*tau,*dm,*D,*dD);
  return params; 
}

RooArgList MLPdfFactory::makeBMixEffDecayParams(TString basename)
{
  RooRealVar *tau = new RooRealVar(basename+"_tau",basename+"_tau",1);
  RooRealVar *dm = new RooRealVar(basename+"_dm",basename+"_dm",0);
  RooRealVar *D = new RooRealVar(basename+"_D",basename+"_D",0);
  RooRealVar *dD = new RooRealVar(basename+"_dD",basename+"_dD",0);
  RooRealVar *mu = new RooRealVar(basename+"_mu",basename+"_mu",0);
  RooRealVar *nu = new RooRealVar(basename+"_nu",basename+"_nu",0);
  RooArgList pars(*tau,*dm,*D,*dD,*mu,*nu);
  return pars;
}

RooAbsPdf* MLPdfFactory::makeBMixEffDecayPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 3, parameters, 6, args, 0)) return 0;

  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *mix = (RooAbsCategory*)obs.at(1);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(2);
  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *dm        = (RooRealVar*)parameters.at(1);
  RooRealVar *D         = (RooRealVar*)parameters.at(2);
  RooRealVar *dD        = (RooRealVar*)parameters.at(3);
  RooRealVar *mu        = (RooRealVar*)parameters.at(4);
  RooRealVar *nu        = (RooRealVar*)parameters.at(5);
  RooFormulaVar *w  = new RooFormulaVar(pdfname+"_w",pdfname+"_w","(1-@0)/2",*D);
  RooFormulaVar *dw = new RooFormulaVar(pdfname+"_dw",pdfname+"_dw","-@0/2",*dD);
  
  RhhBMixEffDecay *pdf = new RhhBMixEffDecay(pdfname,pdfname,*dt,*mix,*tag,*tau,*dm,*w,*dw,*mu,*nu,rf,RhhBMixEffDecay::DoubleSided);
  return pdf;
 
}

RooArgList MLPdfFactory::makeBMixEffDecayMistagParams(TString basename)
{
  RooRealVar *tau = new RooRealVar(basename+"_tau",basename+"_tau",1);
  RooRealVar *dm = new RooRealVar(basename+"_dm",basename+"_dm",0);
  RooRealVar *w = new RooRealVar(basename+"_w",basename+"_w",0.5, 0, 1);
  RooRealVar *dw = new RooRealVar(basename+"_dw",basename+"_dw",0);
  RooRealVar *mu = new RooRealVar(basename+"_mu",basename+"_mu",0);
  RooRealVar *nu = new RooRealVar(basename+"_nu",basename+"_nu",0);
  RooArgList pars(*tau,*dm,*w,*dw,*mu,*nu);
  return pars;
}

RooAbsPdf* MLPdfFactory::makeBMixEffDecayMistagPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 3, parameters, 6, args, 0)) return 0;

  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *mix = (RooAbsCategory*)obs.at(1);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(2);
  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *dm        = (RooRealVar*)parameters.at(1);
  RooRealVar *w         = (RooRealVar*)parameters.at(2);
  RooRealVar *dw        = (RooRealVar*)parameters.at(3);
  RooRealVar *mu        = (RooRealVar*)parameters.at(4);
  RooRealVar *nu        = (RooRealVar*)parameters.at(5);
  
  RhhBMixEffDecay *pdf = new RhhBMixEffDecay(pdfname,pdfname,*dt,*mix,*tag,*tau,*dm,*w,*dw,*mu,*nu,rf,RhhBMixEffDecay::DoubleSided);
  return pdf;
 
}


RooAbsPdf* MLPdfFactory::makeBMixDecayPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args) 
{
  if (!checkArgLength(pdfname, obs, 3, parameters, 4, args, 0)) return 0;

  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *mix = (RooAbsCategory*)obs.at(1);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(2);
  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *dm        = (RooRealVar*)parameters.at(1);
  RooRealVar *D         = (RooRealVar*)parameters.at(2);
  RooRealVar *dD        = (RooRealVar*)parameters.at(3);
  RooFormulaVar *w  = new RooFormulaVar(pdfname+"_w",pdfname+"_w","(1-@0)/2",*D);
  RooFormulaVar *dw = new RooFormulaVar(pdfname+"_dw",pdfname+"_dw","-@0/2",*dD);
  RooBMixDecay *pdf = new RooBMixDecay(pdfname,pdfname,*dt,*mix,*tag,*tau,*dm,*w,*dw,rf,RooBMixDecay::DoubleSided);
  return pdf;
}

RooAbsPdf* MLPdfFactory::makeB0MixEffDecayPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 6, args, 0)) return 0;

  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(1);
  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *dm        = (RooRealVar*)parameters.at(1);
  RooRealVar *D         = (RooRealVar*)parameters.at(2);
  RooRealVar *dD        = (RooRealVar*)parameters.at(3);
  RooRealVar *mu        = (RooRealVar*)parameters.at(4);
  RooRealVar *nu        = (RooRealVar*)parameters.at(5);
  RooCategory *reco      = new RooCategory(pdfname+"_recof",pdfname+"_recof");
  reco->defineType("B0", 1);
  reco->setConstant(kTRUE);
  RooFormulaVar *w  = new RooFormulaVar(pdfname+"_w",pdfname+"_w","(1-@0)/2",*D);
  RooFormulaVar *dw = new RooFormulaVar(pdfname+"_dw",pdfname+"_dw","-@0/2",*dD);
  
  RhhBMixEffDecay *pdf = new RhhBMixEffDecay(pdfname,pdfname,*dt,*reco,*tag,*tau,*dm,*w,*dw,*mu,*nu,rf,RhhBMixEffDecay::DoubleSided);
  return pdf;
 
}

RooAbsPdf* MLPdfFactory::makeB0BarMixEffDecayPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 6, args, 0)) return 0;

  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(1);
  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *dm        = (RooRealVar*)parameters.at(1);
  RooRealVar *D         = (RooRealVar*)parameters.at(2);
  RooRealVar *dD        = (RooRealVar*)parameters.at(3);
  RooRealVar *mu        = (RooRealVar*)parameters.at(4);
  RooRealVar *nu        = (RooRealVar*)parameters.at(5);
  RooCategory *reco      = new RooCategory(pdfname+"_recof",pdfname+"_recof");
  reco->defineType("B0Bar", -1);
  reco->setConstant(kTRUE);
  RooFormulaVar *w  = new RooFormulaVar(pdfname+"_w",pdfname+"_w","(1-@0)/2",*D);
  RooFormulaVar *dw = new RooFormulaVar(pdfname+"_dw",pdfname+"_dw","-@0/2",*dD);
  
  RhhBMixEffDecay *pdf = new RhhBMixEffDecay(pdfname,pdfname,*dt,*reco,*tag,*tau,*dm,*w,*dw,*mu,*nu,rf,RhhBMixEffDecay::DoubleSided);
  return pdf;
 
}

RooArgList MLPdfFactory::makeBifurGaussParams(TString basename)
{
  RooRealVar *mean = new RooRealVar(basename+"_mean", basename+"_mean", 0);
  RooRealVar *sigmaL = new RooRealVar(basename+"_sigmaL", basename+"sigmaL", 1);
  RooRealVar *sigmaR = new RooRealVar(basename+"_sigmaR", basename+"sigmaR", 1);
  return RooArgList(*mean,*sigmaL,*sigmaR);
}

RooAbsPdf* MLPdfFactory::makeBifurGaussPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 3, args, 0)) return 0;
  return new RooBifurGauss(pdfname, pdfname, *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1), *(RooAbsReal*)parameters.at(2));
}

TH1F* MLPdfFactory::makeHisto(const char* name, const char* title, const char* file, const char* xlabel, const char* ylabel)
{
  //
  // This returns a TH1F histogram created from the ascii file given.
  // The file format should be the same as makeHistoAsciiFile uses.
  // Used for reading in DIRC histograms.
  

  FILE *f = fopen(file, "r");
  if (f == 0) {
    cout << "MLPdfFactory::makeHisto: Error opening " << file << endl;
    return 0;
  }
  
  Int_t nbins;
  Double_t xmin, xmax;
  fscanf(f, "%d %lf %lf", &nbins, &xmin, &xmax);
  
  Double_t leftedge[nbins+1], value[nbins], error[nbins];
  for (Int_t i = 0; i < nbins; i++){
    Int_t ncols = fscanf(f, "%lf %lf %lf", &leftedge[i], &value[i], &error[i]);
    if (ncols < 0) {
      cout << "MLPdfFactory::makeHisto: ERROR: file " << file << " has the wrong number of lines or wrong format!" << endl;
      return 0;
    }
  }
  leftedge[nbins] = xmax;
  TH1F *h = new TH1F(name, title, nbins, leftedge);
  for (Int_t i = 0; i < nbins; i++){
    h->SetBinContent(i+1, value[i]);
    h->SetBinError(i+1, error[i]);
  }
  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);  
  return h;
  
}

RooArgList MLPdfFactory::makeDoubleGaussianParams(TString basename)
{
  RooRealVar *mean1 = new RooRealVar(basename+"_mean1", basename+"_mean1", 0);
  RooRealVar *mean2 = new RooRealVar(basename+"_mean2", basename+"_mean2", 0);
  RooRealVar *sigma1 = new RooRealVar(basename+"_sigma1", basename+"sigma1", 1);
  RooRealVar *sigma2 = new RooRealVar(basename+"_sigma2", basename+"sigma2", 1);
  RooRealVar *f1     = new RooRealVar(basename+"_f1",basename+"_f1",0.5);
  
  return RooArgList(*mean1,*mean2,*sigma1,*sigma2,*f1);
}

RooAbsPdf* MLPdfFactory::makeDoubleGaussianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 5, args, 0)) return 0;

  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *m1       = (RooRealVar*)parameters.at(0);
  RooRealVar *m2       = (RooRealVar*)parameters.at(1);
  RooRealVar *s1       = (RooRealVar*)parameters.at(2);
  RooRealVar *s2       = (RooRealVar*)parameters.at(3);
  RooRealVar *f1       = (RooRealVar*)parameters.at(4);
  
  RooGaussian *g1 = new RooGaussian(pdfname+"_g1",pdfname+"_g1", *theobs, *m1, *s1);
  RooGaussian *g2 = new RooGaussian(pdfname+"_g2",pdfname+"_g2", *theobs, *m2, *s2);
  RooAddPdf   *pdf = new RooAddPdf(pdfname, pdfname, RooArgList(*g1, *g2), *f1);
  return pdf;

}

RooArgList MLPdfFactory::makeTripleGaussianParams(TString basename)
{
  RooRealVar *mean1 = new RooRealVar(basename+"_mean1", basename+"_mean1", 0);
  RooRealVar *mean2 = new RooRealVar(basename+"_mean2", basename+"_mean2", 0);
  RooRealVar *mean3 = new RooRealVar(basename+"_mean3", basename+"_mean3", 0);
  RooRealVar *sigma1 = new RooRealVar(basename+"_sigma1", basename+"sigma1", 1);
  RooRealVar *sigma2 = new RooRealVar(basename+"_sigma2", basename+"sigma2", 1);
  RooRealVar *sigma3 = new RooRealVar(basename+"_sigma3", basename+"sigma3", 1);
  RooRealVar *f1     = new RooRealVar(basename+"_f1",basename+"_f1",0.3);
  RooRealVar *f2     = new RooRealVar(basename+"_f2",basename+"_f2",0.3);
  
  return RooArgList(*mean1,*mean2,*mean3,*sigma1,*sigma2,*sigma3,*f1,*f2);
}

RooAbsPdf* MLPdfFactory::makeTripleGaussianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 8, args, 0)) return 0;

  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *m1       = (RooRealVar*)parameters.at(0);
  RooRealVar *m2       = (RooRealVar*)parameters.at(1);
  RooRealVar *m3       = (RooRealVar*)parameters.at(2);
  RooRealVar *s1       = (RooRealVar*)parameters.at(3);
  RooRealVar *s2       = (RooRealVar*)parameters.at(4);
  RooRealVar *s3       = (RooRealVar*)parameters.at(5);
  RooRealVar *f1       = (RooRealVar*)parameters.at(6);
  RooRealVar *f2       = (RooRealVar*)parameters.at(7);
  
  RooGaussian *g1 = new RooGaussian(pdfname+"_g1",pdfname+"_g1", *theobs, *m1, *s1);
  RooGaussian *g2 = new RooGaussian(pdfname+"_g2",pdfname+"_g2", *theobs, *m2, *s2);
  RooGaussian *g3 = new RooGaussian(pdfname+"_g3",pdfname+"_g3", *theobs, *m3, *s3);
  RooAddPdf   *pdf = new RooAddPdf(pdfname, pdfname, RooArgList(*g1, *g2, *g3), RooArgList(*f1,*f2));
  return pdf;

}


RooArgList MLPdfFactory::makeDoubleLandauParams(TString basename)
{
  RooRealVar *mean1 = new RooRealVar(basename+"_mean1", basename+"_mean1", 0);
  RooRealVar *mean2 = new RooRealVar(basename+"_mean2", basename+"_mean2", 0);
  RooRealVar *sigma1 = new RooRealVar(basename+"_sigma1", basename+"sigma1", 1);
  RooRealVar *sigma2 = new RooRealVar(basename+"_sigma2", basename+"sigma2", 1);
  RooRealVar *f1     = new RooRealVar(basename+"_f1",basename+"_f1",0.5);
  
  return RooArgList(*mean1,*mean2,*sigma1,*sigma2,*f1);
}

RooAbsPdf* MLPdfFactory::makeDoubleLandauPdf(TString pdfname, 
					     const RooArgList &obs, 
					     RooArgList parameters, 
					     const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 5, args, 0)) return 0;

  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *m1       = (RooRealVar*)parameters.at(0);
  RooRealVar *m2       = (RooRealVar*)parameters.at(1);
  RooRealVar *s1       = (RooRealVar*)parameters.at(2);
  RooRealVar *s2       = (RooRealVar*)parameters.at(3);
  RooRealVar *f1       = (RooRealVar*)parameters.at(4);
  
  RooLandau *dau1 = new RooLandau(pdfname+"_dau1",pdfname+"_dau1", *theobs, *m1, *s1);
  RooLandau *dau2 = new RooLandau(pdfname+"_dau2",pdfname+"_dau2", *theobs, *m2, *s2);
  RooAddPdf   *pdf = new RooAddPdf(pdfname, pdfname, RooArgList(*dau1, *dau2), *f1);
  return pdf;

}

RooArgList MLPdfFactory::makePulseParams(TString basename)
{
  RooRealVar *mean = new RooRealVar(basename+"_mean", basename+"_mean", 0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma", basename+"sigma", 1);
  RooRealVar *pedSlope = new RooRealVar(basename+"_pedSlope", basename+"pedSlope", 1);
  RooRealVar *f1 = new RooRealVar(basename+"_f1", basename+"f1", 1);
  
  return RooArgList(*mean,*sigma,*pedSlope,*f1);
}

RooAbsPdf* MLPdfFactory::makePulsePdf(TString pdfname, 
				      const RooArgList &obs, 
				      RooArgList parameters, 
				      const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;

  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *m        = (RooRealVar*)parameters.at(0);
  RooRealVar *s        = (RooRealVar*)parameters.at(1);
  RooRealVar *pedSlope = (RooRealVar*)parameters.at(2);
  RooRealVar *f1       = (RooRealVar*)parameters.at(3);
  
  RooLandau *peak = new RooLandau(pdfname+"_peak",pdfname+"_peak", *theobs, *m, *s);
  RooPolynomial *ped = new RooPolynomial(pdfname+"_ped",pdfname+"_peak", *theobs, *pedSlope);
  RooAddPdf   *pdf = new RooAddPdf(pdfname, pdfname, RooArgList(*peak, *ped), *f1);
  return pdf;

}

RooArgList MLPdfFactory::makeVoigtianParams(TString basename)
{
  RooRealVar *mean  = new RooRealVar(basename+"_mean", basename+"_mean", 0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma", basename+"sigma", 1);
  RooRealVar *width = new RooRealVar(basename+"_width", basename+"width", 1);
  return RooArgList(*mean,*sigma,*width);
}

RooAbsPdf* MLPdfFactory::makeVoigtianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 3, args, 0)) return 0;
  return new RooVoigtian(pdfname, pdfname, *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1), *(RooAbsReal*)parameters.at(2));
}

RooArgList MLPdfFactory::makeDoubleVoigtianParams(TString basename)
{
  RooRealVar *mean  = new RooRealVar(basename+"_mean", basename+"_mean", 0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma", basename+"sigma", 1);
  RooRealVar *width1 = new RooRealVar(basename+"_width1", basename+"width1", 1);
  RooRealVar *width2 = new RooRealVar(basename+"_width2", basename+"width2", 1);
  RooRealVar *f1 = new RooRealVar(basename+"_f1", basename+"f1", 1);
  return RooArgList(*mean,*sigma,*width1, *width2, *f1);
}

RooAbsPdf* MLPdfFactory::makeDoubleVoigtianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 5, args, 0)) return 0;
  RooVoigtian* pdf1 = new RooVoigtian(pdfname+"_1", pdfname+"_1", *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1), *(RooAbsReal*)parameters.at(2));
  RooVoigtian* pdf2 = new RooVoigtian(pdfname+"_2", pdfname+"_2", *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1), *(RooAbsReal*)parameters.at(3));
  return new RooAddPdf(pdfname, pdfname, RooArgSet(*pdf1, *pdf2), RooArgSet( *(RooAbsReal*)parameters.at(4) ) );

  return new RooVoigtian(pdfname, pdfname, *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1), *(RooAbsReal*)parameters.at(2));
}

RooArgList MLPdfFactory::makeCryBallParams(TString basename)
{
  RooRealVar *mean  = new RooRealVar(basename+"_mean", basename+"_mean", 0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma", basename+"sigma", 1);
  RooRealVar *alpha = new RooRealVar(basename+"_alpha", basename+"alpha", 1);
  RooRealVar *N     = new RooRealVar(basename+"_N", basename+"N", 1);
  return RooArgList(*mean,*sigma,*alpha,*N);
}

RooAbsPdf* MLPdfFactory::makeCryBallPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;
  return new RooCBShape(pdfname, pdfname, *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1), *(RooAbsReal*)parameters.at(2), *(RooAbsReal*)parameters.at(3));
}

RooArgList MLPdfFactory::makeEffDecay3GParams(TString basename)
{
  RooArgList params(make3GRFParams(basename));
  params.add(makeEffDecayParams(basename));
  return params;
}

RooArgList MLPdfFactory::makeEffDecayParams(TString basename)
{
  RooRealVar *tau = new RooRealVar(basename+"_tau",basename+"_tau",1);
  RooRealVar *mu  = new RooRealVar(basename+"_mu",basename+"_mu",0);
  return RooArgList(*tau,*mu);
}


RooAbsPdf* MLPdfFactory::makeEffDecayPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const RooResolutionModel &rf, const TList &args)
{
  // (1 + tag*recof*mu)*(e^-|t/tau| conv. 3-gaussians)
  if (!checkArgLength(pdfname, obs, 2, parameters, 2, args, 1)) return 0;
  
  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(1);
  
  RooRealVar *tau       = (RooRealVar*)parameters.at(0);
  RooRealVar *mu        = (RooRealVar*)parameters.at(1);
  RooAbsCategory *recof = (RooAbsCategory*)args.At(0);
  
  RooDecay *noEffPdf = new RooDecay(pdfname+"_noEff",pdfname+"_noEff",*dt,*tau,rf, RooDecay::DoubleSided);
  RooGenericPdf *coef = new RooGenericPdf(pdfname+"_coef",pdfname+"_coef", "1+@0*@1*@2", RooArgList(*tag,*recof,*mu));
  RooProdPdf *pdf = new RooProdPdf(pdfname,pdfname,RooArgList(*coef,*noEffPdf));
  return pdf;
}


RooAbsPdf* MLPdfFactory::makeEffDecay3GPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  // (1 + tag*recof*mu)*(e^-|t/tau| conv. 3-gaussians)
  if (!checkArgLength(pdfname, obs, 3, parameters, 10, args, 1)) return 0;
  
  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooRealVar *dte = (RooRealVar*)obs.at(1);  
  RooAbsCategory *tag = (RooAbsCategory*)obs.at(2);  
  RooRealVar *corebias  = (RooRealVar*)parameters.at(0);
  RooRealVar *coresigma = (RooRealVar*)parameters.at(1);
  RooRealVar *tailbias  = (RooRealVar*)parameters.at(2);
  RooRealVar *tailsigma = (RooRealVar*)parameters.at(3);
  RooRealVar *outbias   = (RooRealVar*)parameters.at(4);
  RooRealVar *outsigma  = (RooRealVar*)parameters.at(5);
  RooRealVar *fout      = (RooRealVar*)parameters.at(6);
  RooRealVar *ftail     = (RooRealVar*)parameters.at(7);
  RooRealVar *tau       = (RooRealVar*)parameters.at(8);
  RooRealVar *mu        = (RooRealVar*)parameters.at(9);
  RooAbsCategory *recof = (RooAbsCategory*)args.At(0);
  
  RooArgList rfpars(*corebias,*coresigma,*tailbias,*tailsigma,*outbias,*outsigma,*fout,*ftail);
  RooArgList rfobs(*dt,*dte);
  RooResolutionModel *rf = make3GRF(pdfname, rfobs, rfpars, new TList());
  RooDecay *noEffPdf = new RooDecay(pdfname+"_noEff",pdfname+"_noEff",*dt,*tau,*rf, RooDecay::DoubleSided);
  RooGenericPdf *coef = new RooGenericPdf(pdfname+"_coef",pdfname+"_coef", "1+@0*@1*@2", RooArgList(*tag,*recof,*mu));
  RooProdPdf *pdf = new RooProdPdf(pdfname,pdfname,RooArgList(*coef,*noEffPdf));
  return pdf;
}


RooArgList MLPdfFactory::makeGaussianParams(TString basename)
{
  RooRealVar *mean = new RooRealVar(basename+"_mean",basename+"_mean",0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma",basename+"_sigma",1);
  return RooArgList(*mean,*sigma);
}

RooAbsPdf* MLPdfFactory::makeGaussianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 2, args, 0)) return 0;
  return new RooGaussian(pdfname, pdfname, *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1));
}

RooArgList MLPdfFactory::makeHhDeGaussianParams(TString basename)
{
  RooRealVar *offset = new RooRealVar(basename+"_offset",basename+"_offset",0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma",basename+"_sigma",1);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass1 = new RooRealVar(basename+"_mass1",basename+"_mass1",0);
  RooRealVar *mass2 = new RooRealVar(basename+"_mass2",basename+"_mass2",0);
  
  return RooArgList(*offset,*sigma, *recoMass, *mass1, *mass2);
}

RooAbsPdf* MLPdfFactory::makeHhDeGaussianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 8, parameters, 5, args, 0)) return 0;
   
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px1      = (RooRealVar*)obs.at(1);
  RooRealVar *py1      = (RooRealVar*)obs.at(2);
  RooRealVar *pz1      = (RooRealVar*)obs.at(3);
  RooRealVar *px2      = (RooRealVar*)obs.at(4);
  RooRealVar *py2      = (RooRealVar*)obs.at(5);
  RooRealVar *pz2      = (RooRealVar*)obs.at(6);
  RooRealVar *gamma    = (RooRealVar*)obs.at(7);
  RooRealVar *offset   = (RooRealVar*)parameters.at(0);
  RooRealVar *sigma    = (RooRealVar*)parameters.at(1);
  RooRealVar *massReco = (RooRealVar*)parameters.at(2);
  RooRealVar *mass1    = (RooRealVar*)parameters.at(3);
  RooRealVar *mass2    = (RooRealVar*)parameters.at(4);

  RooFormulaVar *p1    = new RooFormulaVar(pdfname+"_p1", pdfname+"_p1", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px1,*py1,*pz1));
  RooFormulaVar *p2    = new RooFormulaVar(pdfname+"_p2", pdfname+"_p2", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px2,*py2,*pz2));
  RooArgList vars(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset);
  RooFormulaVar *mean = new RooFormulaVar(pdfname+"_mean", pdfname+"_mean", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars);
  RooAbsPdf *pdf = new RooGaussian(pdfname, pdfname, *theobs, *mean, *sigma);
  return pdf;
}

RooArgList MLPdfFactory::makeHhDe2GaussianParams(TString basename)
{
  RooRealVar *offset1 = new RooRealVar(basename+"_offset1",basename+"_offset1",0);
  RooRealVar *sigma1 = new RooRealVar(basename+"_sigma1",basename+"_sigma1",1);
  RooRealVar *offset2 = new RooRealVar(basename+"_offset2",basename+"_offset2",0);
  RooRealVar *sigma2 = new RooRealVar(basename+"_sigma2",basename+"_sigma2",1);
  RooRealVar *frac1  = new RooRealVar(basename+"_frac1",basename+"_frac1",0.1);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass1 = new RooRealVar(basename+"_mass1",basename+"_mass1",0);
  RooRealVar *mass2 = new RooRealVar(basename+"_mass2",basename+"_mass2",0);

  RooArgList pars(*offset1,*sigma1, *offset2, *sigma2);
  pars.add(RooArgList(*frac1, *recoMass, *mass1, *mass2));
  return pars;
}

RooAbsPdf* MLPdfFactory::makeHhDe2GaussianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 8, parameters, 8, args, 0)) return 0;
   
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px1      = (RooRealVar*)obs.at(1);
  RooRealVar *py1      = (RooRealVar*)obs.at(2);
  RooRealVar *pz1      = (RooRealVar*)obs.at(3);
  RooRealVar *px2      = (RooRealVar*)obs.at(4);
  RooRealVar *py2      = (RooRealVar*)obs.at(5);
  RooRealVar *pz2      = (RooRealVar*)obs.at(6);
  RooRealVar *gamma    = (RooRealVar*)obs.at(7);
  RooRealVar *offset1   = (RooRealVar*)parameters.at(0);
  RooRealVar *sigma1    = (RooRealVar*)parameters.at(1);
  RooRealVar *offset2   = (RooRealVar*)parameters.at(2);
  RooRealVar *sigma2    = (RooRealVar*)parameters.at(3);
  RooRealVar *frac1     = (RooRealVar*)parameters.at(4);
  RooRealVar *massReco = (RooRealVar*)parameters.at(5);
  RooRealVar *mass1    = (RooRealVar*)parameters.at(6);
  RooRealVar *mass2    = (RooRealVar*)parameters.at(7);

  RooFormulaVar *p1    = new RooFormulaVar(pdfname+"_p1", pdfname+"_p1", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px1,*py1,*pz1));
  RooFormulaVar *p2    = new RooFormulaVar(pdfname+"_p2", pdfname+"_p2", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px2,*py2,*pz2));
  RooArgList vars1(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset1);
  RooFormulaVar *mean1 = new RooFormulaVar(pdfname+"_mean1", pdfname+"_mean1", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars1);
  RooArgList vars2(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset2);
  RooFormulaVar *mean2 = new RooFormulaVar(pdfname+"_mean2", pdfname+"_mean2", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars2);

  RooGaussian *pdf1 = new RooGaussian(pdfname+"_g1", pdfname+"_g1", *theobs, *mean1, *sigma1);
  RooGaussian *pdf2 = new RooGaussian(pdfname+"_g2", pdfname+"_g2", *theobs, *mean2, *sigma2);
  RooAbsPdf *pdf = new RooAddPdf(pdfname,pdfname,RooArgList(*pdf1, *pdf2), RooArgList(*frac1));
  
  return pdf;
}

RooArgList MLPdfFactory::makeHhDe3GaussianParams(TString basename)
{
  RooRealVar *offset1 = new RooRealVar(basename+"_offset1",basename+"_offset1",0);
  RooRealVar *sigma1 = new RooRealVar(basename+"_sigma1",basename+"_sigma1",1);
  RooRealVar *offset2 = new RooRealVar(basename+"_offset2",basename+"_offset2",0);
  RooRealVar *sigma2 = new RooRealVar(basename+"_sigma2",basename+"_sigma2",1);
  RooRealVar *offset3 = new RooRealVar(basename+"_offset3",basename+"_offset3",0);
  RooRealVar *sigma3 = new RooRealVar(basename+"_sigma3",basename+"_sigma3",1);
  RooRealVar *frac1  = new RooRealVar(basename+"_frac1",basename+"_frac1",0.1);
  RooRealVar *frac2  = new RooRealVar(basename+"_frac2",basename+"_frac2",0.1);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass1 = new RooRealVar(basename+"_mass1",basename+"_mass1",0);
  RooRealVar *mass2 = new RooRealVar(basename+"_mass2",basename+"_mass2",0);

  RooArgList pars(*offset1,*sigma1, *offset2, *sigma2, *offset3, *sigma3);
  pars.add(RooArgList(*frac1, *frac2, *recoMass, *mass1, *mass2));
  return pars;
}

RooAbsPdf* MLPdfFactory::makeHhDe3GaussianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 8, parameters, 11, args, 0)) return 0;
   
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px1      = (RooRealVar*)obs.at(1);
  RooRealVar *py1      = (RooRealVar*)obs.at(2);
  RooRealVar *pz1      = (RooRealVar*)obs.at(3);
  RooRealVar *px2      = (RooRealVar*)obs.at(4);
  RooRealVar *py2      = (RooRealVar*)obs.at(5);
  RooRealVar *pz2      = (RooRealVar*)obs.at(6);
  RooRealVar *gamma    = (RooRealVar*)obs.at(7);
  RooRealVar *offset1   = (RooRealVar*)parameters.at(0);
  RooRealVar *sigma1    = (RooRealVar*)parameters.at(1);
  RooRealVar *offset2   = (RooRealVar*)parameters.at(2);
  RooRealVar *sigma2    = (RooRealVar*)parameters.at(3);
  RooRealVar *offset3   = (RooRealVar*)parameters.at(4);
  RooRealVar *sigma3    = (RooRealVar*)parameters.at(5);
  RooRealVar *frac1     = (RooRealVar*)parameters.at(6);
  RooRealVar *frac2     = (RooRealVar*)parameters.at(7);  
  RooRealVar *massReco = (RooRealVar*)parameters.at(8);
  RooRealVar *mass1    = (RooRealVar*)parameters.at(9);
  RooRealVar *mass2    = (RooRealVar*)parameters.at(10);

  RooFormulaVar *p1    = new RooFormulaVar(pdfname+"_p1", pdfname+"_p1", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px1,*py1,*pz1));
  RooFormulaVar *p2    = new RooFormulaVar(pdfname+"_p2", pdfname+"_p2", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px2,*py2,*pz2));
  RooArgList vars1(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset1);
  RooFormulaVar *mean1 = new RooFormulaVar(pdfname+"_mean1", pdfname+"_mean1", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars1);
  RooArgList vars2(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset2);
  RooFormulaVar *mean2 = new RooFormulaVar(pdfname+"_mean2", pdfname+"_mean2", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars2);
  RooArgList vars3(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset3);
  RooFormulaVar *mean3 = new RooFormulaVar(pdfname+"_mean3", pdfname+"_mean3", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars3);
  RooGaussian *pdf1 = new RooGaussian(pdfname+"_g1", pdfname+"_g1", *theobs, *mean1, *sigma1);
  RooGaussian *pdf2 = new RooGaussian(pdfname+"_g2", pdfname+"_g2", *theobs, *mean2, *sigma2);
  RooGaussian *pdf3 = new RooGaussian(pdfname+"_g3", pdfname+"_g3", *theobs, *mean3, *sigma3);
  RooAbsPdf *pdf = new RooAddPdf(pdfname,pdfname,RooArgList(*pdf1, *pdf2, *pdf3), RooArgList(*frac1, *frac2));
  
  return pdf;
}



RooArgList MLPdfFactory::makeHhDeCryBallParams(TString basename)
{
  RooRealVar *offset = new RooRealVar(basename+"_offset",basename+"_offset",0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma",basename+"_sigma",1);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass1 = new RooRealVar(basename+"_mass1",basename+"_mass1",0);
  RooRealVar *mass2 = new RooRealVar(basename+"_mass2",basename+"_mass2",0);
  RooRealVar *alpha = new RooRealVar(basename+"_alpha", basename+"alpha", 1);
  RooRealVar *N     = new RooRealVar(basename+"_N", basename+"N", 1);
  
  return RooArgList(*offset,*sigma, *alpha, *N, *recoMass, *mass1, *mass2);
}

RooAbsPdf *MLPdfFactory::makeHhDeCryBallPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args) 
{
  if (!checkArgLength(pdfname, obs, 8, parameters, 7, args, 0)) return 0;
   
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px1      = (RooRealVar*)obs.at(1);
  RooRealVar *py1      = (RooRealVar*)obs.at(2);
  RooRealVar *pz1      = (RooRealVar*)obs.at(3);
  RooRealVar *px2      = (RooRealVar*)obs.at(4);
  RooRealVar *py2      = (RooRealVar*)obs.at(5);
  RooRealVar *pz2      = (RooRealVar*)obs.at(6);
  RooRealVar *gamma    = (RooRealVar*)obs.at(7);
  RooRealVar *offset   = (RooRealVar*)parameters.at(0);
  RooRealVar *sigma    = (RooRealVar*)parameters.at(1);
  RooRealVar *alpha    = (RooRealVar*)parameters.at(2);
  RooRealVar *N        = (RooRealVar*)parameters.at(3);
  RooRealVar *massReco = (RooRealVar*)parameters.at(4);
  RooRealVar *mass1    = (RooRealVar*)parameters.at(5);
  RooRealVar *mass2    = (RooRealVar*)parameters.at(6);

  RooFormulaVar *p1    = new RooFormulaVar(pdfname+"_p1", pdfname+"_p1", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px1,*py1,*pz1));
  RooFormulaVar *p2    = new RooFormulaVar(pdfname+"_p2", pdfname+"_p2", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px2,*py2,*pz2));
  RooArgList vars(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset);
  RooFormulaVar *mean = new RooFormulaVar(pdfname+"_mean", pdfname+"_mean", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars);
  RooAbsPdf *pdf = new RooCBShape(pdfname, pdfname, *theobs, *mean, *sigma, *alpha, *N);

  return pdf;
  
}

RooArgList MLPdfFactory::makeHhDeCryBallGParams(TString basename)
{
  RooRealVar *cboffset = new RooRealVar(basename+"_cboffset",basename+"_cboffset",0);
  RooRealVar *cbsigma = new RooRealVar(basename+"_cbsigma",basename+"_cbsigma",1);
  RooRealVar *alpha = new RooRealVar(basename+"_alpha", basename+"alpha", 1);
  RooRealVar *N     = new RooRealVar(basename+"_N", basename+"N", 1.05);
  RooRealVar *goffset = new RooRealVar(basename+"_goffset",basename+"_goffset",0);
  RooRealVar *gsigma = new RooRealVar(basename+"_gsigma",basename+"_gsigma",1);
  RooRealVar *fraccb = new RooRealVar(basename+"_fraccb",basename+"_fraccb",0.5);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass1 = new RooRealVar(basename+"_mass1",basename+"_mass1",0);
  RooRealVar *mass2 = new RooRealVar(basename+"_mass2",basename+"_mass2",0);
 
  RooArgList pars(*cboffset, *cbsigma, *alpha);
  pars.add(RooArgList(*N, *goffset, *gsigma, *fraccb, *recoMass, *mass1, *mass2));
  return pars;
}

RooAbsPdf* MLPdfFactory::makeHhDeCryBallGPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 8, parameters, 10, args, 0)) return 0;
 
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px1      = (RooRealVar*)obs.at(1);
  RooRealVar *py1      = (RooRealVar*)obs.at(2);
  RooRealVar *pz1      = (RooRealVar*)obs.at(3);
  RooRealVar *px2      = (RooRealVar*)obs.at(4);
  RooRealVar *py2      = (RooRealVar*)obs.at(5);
  RooRealVar *pz2      = (RooRealVar*)obs.at(6);
  RooRealVar *gamma    = (RooRealVar*)obs.at(7);
  RooRealVar *cboffset = (RooRealVar*)parameters.at(0);
  RooRealVar *cbsigma  = (RooRealVar*)parameters.at(1);
  RooRealVar *alpha    = (RooRealVar*)parameters.at(2);
  RooRealVar *N        = (RooRealVar*)parameters.at(3);
  RooRealVar *goffset  = (RooRealVar*)parameters.at(4);
  RooRealVar *gsigma   = (RooRealVar*)parameters.at(5);
  RooRealVar *fraccb   = (RooRealVar*)parameters.at(6);
  RooRealVar *massReco = (RooRealVar*)parameters.at(7);
  RooRealVar *mass1    = (RooRealVar*)parameters.at(8);
  RooRealVar *mass2    = (RooRealVar*)parameters.at(9);

  RooFormulaVar *p1    = new RooFormulaVar(pdfname+"_p1", pdfname+"_p1", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px1,*py1,*pz1));
  RooFormulaVar *p2    = new RooFormulaVar(pdfname+"_p2", pdfname+"_p2", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px2,*py2,*pz2));
  RooArgList gvars(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *goffset);
  RooFormulaVar *gmean = new RooFormulaVar(pdfname+"_gmean", pdfname+"_gmean", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",gvars);
  RooArgList cbvars(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *cboffset);
  RooFormulaVar *cbmean = new RooFormulaVar(pdfname+"_cbmean", pdfname+"_cbmean", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",cbvars);
  
  RooGaussian *gaus = new RooGaussian(pdfname+"_gaus",pdfname+"_gaus",*theobs, *gmean, *gsigma);
  RooCBShape  *cb   = new RooCBShape(pdfname+"_cb",pdfname+"_cb",*theobs, *cbmean, *cbsigma, *alpha, *N);
  RooAbsPdf *pdf = new RooAddPdf(pdfname,pdfname,RooArgList(*cb,*gaus),*fraccb);

  return pdf;
}


RooArgList MLPdfFactory::makeExclJetParams(TString basename)
{
  RooRealVar *njets = new RooRealVar(basename+"_njets", basename+"_njets", 0);

  return RooArgList(*njets);
}

RooAbsPdf *MLPdfFactory::makeExclJetPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 1, args, 0)) return 0;

  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *njets    = (RooRealVar*)parameters.at(0);

  return new MLExclJetPdf(pdfname,pdfname, *theobs, *njets);
}

RooArgList MLPdfFactory::makeBtagParams(TString basename)
{
  RooRealVar *njet   = new RooRealVar(basename+"_njets", basename+"_njets", 0);
  RooRealVar *nb     = new RooRealVar(basename+"_nb", basename+"_nb", 0);
  RooRealVar *eb     = new RooRealVar(basename+"_eb", basename+"_eb", 0);
  RooRealVar *enob   = new RooRealVar(basename+"_enob", basename+"_enob", 0);

  return RooArgList(*njet, *nb, *eb, *enob);
}

RooAbsPdf *MLPdfFactory::makeBtagPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;

  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *njets    = (RooRealVar*)parameters.at(0);
  RooRealVar *nb       = (RooRealVar*)parameters.at(1);
  RooRealVar *eb       = (RooRealVar*)parameters.at(2);
  RooRealVar *enob     = (RooRealVar*)parameters.at(3);

  return new VecbosBtagPdf(pdfname, pdfname, *theobs, *njets, *nb, *eb, *enob);
}

RooArgList MLPdfFactory::makeHhDeCruijffParams(TString basename)
{
  RooRealVar *offset = new RooRealVar(basename+"_offset",basename+"_offset",0);
  RooRealVar *sigmaL  = new RooRealVar(basename+"_sigmaL",basename+"_sigmaL",1);
  RooRealVar *sigmaR  = new RooRealVar(basename+"_sigmaR",basename+"_sigmaR",1);
  RooRealVar *alphaL  = new RooRealVar(basename+"_alphaL",basename+"_alphaL",1);
  RooRealVar *alphaR  = new RooRealVar(basename+"_alphaR",basename+"_alphaR",1);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass1   = new RooRealVar(basename+"_mass1",basename+"_mass1",0);
  RooRealVar *mass2   = new RooRealVar(basename+"_mass2",basename+"_mass2",0);

  return RooArgList(*offset,*sigmaL,*sigmaR,*alphaL,*alphaR,*recoMass,*mass1,*mass2);
  
}

RooAbsPdf *MLPdfFactory::makeHhDeCruijffPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 8, parameters, 8, args, 0)) return 0;

  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px1      = (RooRealVar*)obs.at(1);
  RooRealVar *py1      = (RooRealVar*)obs.at(2);
  RooRealVar *pz1      = (RooRealVar*)obs.at(3);
  RooRealVar *px2      = (RooRealVar*)obs.at(4);
  RooRealVar *py2      = (RooRealVar*)obs.at(5);
  RooRealVar *pz2      = (RooRealVar*)obs.at(6);
  RooRealVar *gamma    = (RooRealVar*)obs.at(7);
  RooRealVar *offset   = (RooRealVar*)parameters.at(0);
  RooRealVar *sigmaL    = (RooRealVar*)parameters.at(1);
  RooRealVar *sigmaR    = (RooRealVar*)parameters.at(1);
  RooRealVar *alphaL    = (RooRealVar*)parameters.at(2);
  RooRealVar *alphaR    = (RooRealVar*)parameters.at(2);
  RooRealVar *massReco = (RooRealVar*)parameters.at(4);
  RooRealVar *mass1    = (RooRealVar*)parameters.at(5);
  RooRealVar *mass2    = (RooRealVar*)parameters.at(6);

  RooFormulaVar *p1    = new RooFormulaVar(pdfname+"_p1", pdfname+"_p1", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px1,*py1,*pz1));
  RooFormulaVar *p2    = new RooFormulaVar(pdfname+"_p2", pdfname+"_p2", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px2,*py2,*pz2));
  RooArgList vars(*gamma, *massReco, *p1, *mass1, *p2, *mass2, *offset);
  RooFormulaVar *mean = new RooFormulaVar(pdfname+"_mean", pdfname+"_mean", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2)+sqrt(@1*@1+@4*@4)-sqrt(@5*@5+@4*@4))+@6",vars);
  RooAbsPdf *pdf = new RhhCruijffPdf(pdfname, pdfname, *theobs, *mean, *sigmaL, *sigmaR, *alphaL, *alphaR);

  return pdf;


}


RooArgList MLPdfFactory::makeDeDoubleGaussianParams(TString basename)
{
  RooRealVar *mean   = new RooRealVar(basename+"_mean",basename+"_mean",0);
  RooRealVar *offset = new RooRealVar(basename+"_offset",basename+"_offset",0);
  RooRealVar *sigma1 = new RooRealVar(basename+"_sigma1",basename+"_sigma1",1);
  RooRealVar *sigma2 = new RooRealVar(basename+"_sigma2",basename+"_sigma2",1);
  RooRealVar *frac   = new RooRealVar(basename+"_frac",basename+"_frac",1);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass  = new RooRealVar(basename+"_mass",basename+"_mass",0);
  
  return RooArgList(*mean, *offset,*sigma1, *sigma2, *frac,*recoMass, *mass);

}

RooAbsPdf* MLPdfFactory::makeDeDoubleGaussianPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 5, parameters, 7, args, 0)) return 0;
   
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px       = (RooRealVar*)obs.at(1);
  RooRealVar *py       = (RooRealVar*)obs.at(2);
  RooRealVar *pz       = (RooRealVar*)obs.at(3);
  RooRealVar *gamma    = (RooRealVar*)obs.at(4);
  RooRealVar *mean     = (RooRealVar*)parameters.at(0);
  RooRealVar *offset   = (RooRealVar*)parameters.at(1);
  RooRealVar *sigma1   = (RooRealVar*)parameters.at(2);
  RooRealVar *sigma2   = (RooRealVar*)parameters.at(3);
  RooRealVar *frac     = (RooRealVar*)parameters.at(4);
  RooRealVar *massReco = (RooRealVar*)parameters.at(5);
  RooRealVar *mass     = (RooRealVar*)parameters.at(6);

  RooFormulaVar *p     = new RooFormulaVar(pdfname+"_p", pdfname+"_p", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px,*py,*pz));
  RooArgList vars(*gamma, *massReco, *p, *mass, *mean);
  RooFormulaVar *mean1 = new RooFormulaVar(pdfname+"_mean1", pdfname+"_mean1", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2))+@4",vars);
  RooFormulaVar *mean2 = new RooFormulaVar(pdfname+"_mean2", pdfname+"_mean2", "@0+@1",RooArgList(*mean1,*offset));

  RooGaussian *g1 = new RooGaussian(pdfname+"_g1",pdfname+"_g1", *theobs, *mean1, *sigma1);
  RooGaussian *g2 = new RooGaussian(pdfname+"_g2",pdfname+"_g2", *theobs, *mean2, *sigma2);
  RooAddPdf   *pdf = new RooAddPdf(pdfname, pdfname, RooArgList(*g1, *g2), *frac);

  return pdf;
}

RooArgList MLPdfFactory::makeHh0DeCryBallParams(TString basename)
{
  RooRealVar *offset = new RooRealVar(basename+"_offset",basename+"_offset",0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma",basename+"_sigma",1);
  RooRealVar *recoMass = new RooRealVar(basename+"_recoMass",basename+"_recoMass",0);
  RooRealVar *mass1 = new RooRealVar(basename+"_mass1",basename+"_mass1",0);
  RooRealVar *alpha = new RooRealVar(basename+"_alpha", basename+"alpha", 1);
  RooRealVar *N     = new RooRealVar(basename+"_N", basename+"N", 1);
                                                                                                                             
  return RooArgList(*offset,*sigma, *alpha, *N, *recoMass, *mass1);
}

RooAbsPdf* MLPdfFactory::makeHh0DeCryBallPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 5, parameters, 6, args, 0)) return 0;
  
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  RooRealVar *px1      = (RooRealVar*)obs.at(1);
  RooRealVar *py1      = (RooRealVar*)obs.at(2);
  RooRealVar *pz1      = (RooRealVar*)obs.at(3);
  RooRealVar *gamma    = (RooRealVar*)obs.at(4);
  RooRealVar *offset   = (RooRealVar*)parameters.at(0);
  RooRealVar *sigma    = (RooRealVar*)parameters.at(1);
  RooRealVar *alpha    = (RooRealVar*)parameters.at(2);
  RooRealVar *N        = (RooRealVar*)parameters.at(3);
  RooRealVar *massReco = (RooRealVar*)parameters.at(4);
  RooRealVar *mass1    = (RooRealVar*)parameters.at(5);
  
  RooFormulaVar *p1    = new RooFormulaVar(pdfname+"_p1", pdfname+"_p1", "sqrt(@0*@0+@1*@1+@2*@2)",RooArgList(*px1,*py1,*pz1));
  RooArgList vars(*gamma, *massReco, *p1, *mass1, *offset);
  RooFormulaVar *mean = new RooFormulaVar(pdfname+"_mean", pdfname+"_mean", "@0*(sqrt(@1*@1+@2*@2)-sqrt(@3*@3+@2*@2))+@4",vars);
  RooAbsPdf *pdf = new RooCBShape(pdfname, pdfname, *theobs, *mean, *sigma, *alpha, *N);
  return pdf;
}



RooAbsPdf*  MLPdfFactory::makeInvArgusPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 2, args, 0)) return 0;
  RooRealVar *pstar  = (RooRealVar*)obs.at(0);
  RooRealVar *cutoff = (RooRealVar*)parameters.at(0);
  RooRealVar *shape  = (RooRealVar*)parameters.at(1);
  RooFormulaVar *mesFormula = new RooFormulaVar(pdfname+"_mesfunc",pdfname+"_mesfunc","sqrt(@0*@0-@1*@1)", RooArgSet(*cutoff, *pstar));
  RooAbsPdf *pdf = new RooArgusBG(pdfname,pdfname,*mesFormula, *cutoff, *shape);
  return pdf;
}

RooArgList MLPdfFactory::makeLandauParams(TString basename)
{
  RooRealVar *mean  = new RooRealVar(basename+"_mean" ,basename+"_mean",0);
  RooRealVar *sigma = new RooRealVar(basename+"_sigma",basename+"_sigma",1);
  return RooArgList(*mean,*sigma);
}

RooAbsPdf *MLPdfFactory::makeLandauPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 2, args, 0)) return 0;
  return new RooLandau(pdfname,pdfname,*(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1));
}

RooArgList MLPdfFactory::makeLineParams(TString basename)
{
  RooRealVar *p1 = new RooRealVar(basename+"_p1",basename+"_p1",1);
  return RooArgList(*p1);
}

RooAbsPdf* MLPdfFactory::makeLinePdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 1, args, 0)) return 0;
  return new RooPolynomial(pdfname, pdfname, *(RooAbsReal*)obs.at(0), RooArgList(*(RooAbsReal*)parameters.at(0)));
}

RooArgList MLPdfFactory::makePoly2Params(TString basename)
{
  RooRealVar *p1 = new RooRealVar(basename+"_p1",basename+"_p1",1);
  RooRealVar *p2 = new RooRealVar(basename+"_p2",basename+"_p2",1);
  return RooArgList(*p1,*p2);
}

RooAbsPdf* MLPdfFactory::makePoly2Pdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 2, args, 0)) return 0;
  return new RooPolynomial(pdfname, pdfname, *(RooAbsReal*)obs.at(0), RooArgList(*(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1)));
}

RooArgList MLPdfFactory::makeBadSigDTParams(TString basename) 
{
  RooRealVar *tauB0 = new RooRealVar(basename+"_tauB0",basename+"_tauB0",1);
  RooRealVar *dm = new RooRealVar(basename+"_dm",basename+"_dm",1);
  RooRealVar *D = new RooRealVar(basename+"_D",basename+"_D",1);
  RooRealVar *dD = new RooRealVar(basename+"_dD",basename+"_dD",0);
  RooRealVar *C  = new RooRealVar(basename+"_C",basename+"_C",0);
  RooRealVar *mu = new RooRealVar(basename+"_mu",basename+"_mu",0);
  RooRealVar *cp = new RooRealVar(basename+"_cp",basename+"_cp",1);
  
  return RooArgList(*tauB0,*dm,*D,*dD,*C,*mu,*cp);
}

RooAbsPdf* MLPdfFactory::makeBadSigDTPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 7, args, 0)) return 0;

  RooCategory *tag  = (RooCategory*) obs.at(0);

  RooRealVar  *tau         = (RooRealVar*)  parameters.at(0);
  RooRealVar  *dm          = (RooRealVar*)  parameters.at(1);
  RooRealVar  *D           = (RooRealVar*)  parameters.at(2);
  RooRealVar  *dD          = (RooRealVar*)  parameters.at(3);
  RooRealVar  *C           = (RooRealVar*)  parameters.at(4);
  RooRealVar  *mu          = (RooRealVar*)  parameters.at(5);
  RooRealVar  *cp          = (RooRealVar*)  parameters.at(6);

  RooGenericPdf *pdf = new RooGenericPdf(pdfname,pdfname, 
					    "((@3*@0*@1+(1+@0*@2/2))-@7*(@0*@1+@3*(1+@0*@2/2))/(1+(@5*@6)^2))*(1/(1-(@3*@7*(1/(1+(@5*@6)^2)))))/2" ,
					    RooArgSet(*tag, *D, *dD, *mu, *cp, *tau, *dm, *C) );

  return pdf;
}

RooArgList MLPdfFactory::makeBadBkgDTParams(TString basename) 
{
  RooRealVar *mu = new RooRealVar(basename+"_mu",basename+"_mu",0);
  return RooArgList(*mu);
}

RooAbsPdf* MLPdfFactory::makeBadBkgDTPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  RooCategory *tag  = (RooCategory*) obs.at(0);
  RooRealVar  *mu   = (RooRealVar*)  parameters.at(0);
  
  //RooArgSet *emptyset=new RooArgSet();
  // RooGenericPdf *pdf = new RooGenericPdf(pdfname, pdfname, "1"   ,*emptyset );
  RooGenericPdf *pdf = new RooGenericPdf(pdfname, pdfname, "1+@0*@1", RooArgSet(*tag,*mu)) ;
  return pdf;
}

RooAbsPdf* MLPdfFactory::makeNoPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 0, args, 0)) return 0;
  return new RooPolynomial(pdfname, pdfname, *(RooAbsReal*)obs.at(0), RooArgList());
}


RooAbsPdf* MLPdfFactory::makePrebuiltPdf(const TList &args)
{
  // This is a special function that really just returns the pdf that you gave it.  Useful when you're using MLFit
  // to build PDFs, but one of your PDFs is already built (like the DIRC PDF from the Dirc builder).
  
  if (args.GetSize() != 1) {
    cout << "MLPdfFactory Error: Trying to return a prebuilt PDF, but my arg list is: " << endl;
    args.Print();
    return 0;
  }
  return (RooAbsPdf*)args.At(0);
  
}



//=============================================
// Resolution functions...
//=============================================
RooArgList MLPdfFactory::makeScaled2GParams(TString basename)
{
  RooRealVar *corebias = new RooRealVar(basename+"_corebias",basename+"_corebias",0);
  RooRealVar *coresigma = new RooRealVar(basename+"_coresigma",basename+"_coresigma",1);
  RooRealVar *outbias = new RooRealVar(basename+"_outbias",basename+"_outbias",0);
  RooRealVar *outsigma = new RooRealVar(basename+"_outsigma",basename+"_outsigma",1);
  RooRealVar *fout = new RooRealVar(basename+"_fout",basename+"_fout",0);
  return RooArgList(*corebias,*coresigma,*outbias,*outsigma,*fout);
}

RooResolutionModel* MLPdfFactory::makeScaled2GRF(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args) 
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 5, args, 0)) return 0;
  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooRealVar *dte = (RooRealVar*)obs.at(1);
  
  RooRealVar *corebias  = (RooRealVar*)parameters.at(0);
  RooRealVar *coresigma = (RooRealVar*)parameters.at(1);
  RooRealVar *outbias   = (RooRealVar*)parameters.at(2);
  RooRealVar *outsigma  = (RooRealVar*)parameters.at(3);
  RooRealVar *fout      = (RooRealVar*)parameters.at(4);
  
  RooGaussModel *core = new RooGaussModel(pdfname+"_core",pdfname+"_core",*dt,*corebias,*coresigma,*dte);
  RooGaussModel *out = new RooGaussModel(pdfname+"_out",pdfname+"_out",*dt,*outbias,*outsigma);
  RooAddModel *rf = new RooAddModel(pdfname,pdfname,RooArgList(*out,*core),RooArgList(*fout));
  return rf;	   
}

RooArgList MLPdfFactory::makeScaled3GParams(TString basename)
{
  RooRealVar *corebias = new RooRealVar(basename+"_corebias",basename+"_corebias",0);
  RooRealVar *coresigma = new RooRealVar(basename+"_coresigma",basename+"_coresigma",1);
  RooRealVar *tailbias = new RooRealVar(basename+"_tailbias",basename+"_tailbias",0);
  RooRealVar *tailsigma = new RooRealVar(basename+"_tailsigma",basename+"_tailsigma",1);
  RooRealVar *outbias = new RooRealVar(basename+"_outbias",basename+"_outbias",0);
  RooRealVar *outsigma = new RooRealVar(basename+"_outsigma",basename+"_outsigma",1);
  RooRealVar *fout = new RooRealVar(basename+"_fout",basename+"_fout",0);
  RooRealVar *ftail = new RooRealVar(basename+"_ftail",basename+"_ftail",0);
  return RooArgList(*corebias,*coresigma,*tailbias,*tailsigma,*outbias,*outsigma,*fout,*ftail);
}

RooResolutionModel* MLPdfFactory::makeScaled3GRF(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args) 
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 8, args, 0)) return 0;
  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooRealVar *dte = (RooRealVar*)obs.at(1);
  
  RooRealVar *corebias  = (RooRealVar*)parameters.at(0);
  RooRealVar *coresigma = (RooRealVar*)parameters.at(1);
  RooRealVar *tailbias  = (RooRealVar*)parameters.at(2);
  RooRealVar *tailsigma = (RooRealVar*)parameters.at(3);
  RooRealVar *outbias   = (RooRealVar*)parameters.at(4);
  RooRealVar *outsigma  = (RooRealVar*)parameters.at(5);
  RooRealVar *fout      = (RooRealVar*)parameters.at(6);
  RooRealVar *ftail     = (RooRealVar*)parameters.at(7);
  
  RooGaussModel *core = new RooGaussModel(pdfname+"_core",pdfname+"_core",*dt,*corebias,*coresigma,*dte);
  RooGaussModel *tail = new RooGaussModel(pdfname+"_tail",pdfname+"_tail",*dt,*tailbias,*tailsigma,*dte);
  RooGaussModel *out = new RooGaussModel(pdfname+"_out",pdfname+"_out",*dt,*outbias,*outsigma);
  RooAddModel *rf = new RooAddModel(pdfname,pdfname,RooArgList(*out,*tail,*core),RooArgList(*fout,*ftail));
  return rf;	   
}


RooArgList MLPdfFactory::make3GRFParams(TString basename)
{
  RooRealVar *corebias = new RooRealVar(basename+"_corebias",basename+"_corebias",0);
  RooRealVar *coresigma = new RooRealVar(basename+"_coresigma",basename+"_coresigma",1);
  RooRealVar *tailbias = new RooRealVar(basename+"_tailbias",basename+"_tailbias",0);
  RooRealVar *tailsigma = new RooRealVar(basename+"_tailsigma",basename+"_tailsigma",1);
  RooRealVar *outbias = new RooRealVar(basename+"_outbias",basename+"_outbias",0);
  RooRealVar *outsigma = new RooRealVar(basename+"_outsigma",basename+"_outsigma",1);
  RooRealVar *fout = new RooRealVar(basename+"_fout",basename+"_fout",0);
  RooRealVar *ftail = new RooRealVar(basename+"_ftail",basename+"_ftail",0);
  return RooArgList(*corebias,*coresigma,*tailbias,*tailsigma,*outbias,*outsigma,*fout,*ftail);
}


RooResolutionModel* MLPdfFactory::make3GRF(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args) 
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 8, args, 0)) return 0;
  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooRealVar *dte = (RooRealVar*)obs.at(1);
  
  RooRealVar *corebias  = (RooRealVar*)parameters.at(0);
  RooRealVar *coresigma = (RooRealVar*)parameters.at(1);
  RooRealVar *tailbias  = (RooRealVar*)parameters.at(2);
  RooRealVar *tailsigma = (RooRealVar*)parameters.at(3);
  RooRealVar *outbias   = (RooRealVar*)parameters.at(4);
  RooRealVar *outsigma  = (RooRealVar*)parameters.at(5);
  RooRealVar *fout      = (RooRealVar*)parameters.at(6);
  RooRealVar *ftail     = (RooRealVar*)parameters.at(7);
  
  RooGaussModel *core = new RooGaussModel(pdfname+"_core",pdfname+"_core",*dt,*corebias,*coresigma,*dte);
  RooGaussModel *tail = new RooGaussModel(pdfname+"_tail",pdfname+"_tail",*dt,*tailbias,*tailsigma,*dte);
  RooGaussModel *out = new RooGaussModel(pdfname+"_out",pdfname+"_out",*dt,*outbias,*outsigma);
  RooAddModel *rf = new RooAddModel(pdfname+"_3GRF",pdfname+"_3GRF",RooArgList(*out,*tail,*core),RooArgList(*fout,*ftail));
  return rf;	   
}

RooArgList MLPdfFactory::make2GRFParams(TString basename)
{
  RooRealVar *corebias = new RooRealVar(basename+"_corebias",basename+"_corebias",0);
  RooRealVar *coresigma = new RooRealVar(basename+"_coresigma",basename+"_coresigma",1);
  RooRealVar *outbias = new RooRealVar(basename+"_outbias",basename+"_outbias",0);
  RooRealVar *outsigma = new RooRealVar(basename+"_outsigma",basename+"_outsigma",1);
  RooRealVar *fout = new RooRealVar(basename+"_fout",basename+"_fout",0);
  return RooArgList(*corebias,*coresigma,*outbias,*outsigma,*fout);
}


RooResolutionModel* MLPdfFactory::make2GRF(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args) 
{
  if (!checkArgLength(pdfname, obs, 2, parameters, 5, args, 0)) return 0;
  RooRealVar *dt  = (RooRealVar*)obs.at(0);
  RooRealVar *dte = (RooRealVar*)obs.at(1);
  
  RooRealVar *corebias  = (RooRealVar*)parameters.at(0);
  RooRealVar *coresigma = (RooRealVar*)parameters.at(1);
  RooRealVar *outbias   = (RooRealVar*)parameters.at(2);
  RooRealVar *outsigma  = (RooRealVar*)parameters.at(3);
  RooRealVar *fout      = (RooRealVar*)parameters.at(4);
  
  RooGaussModel *core = new RooGaussModel(pdfname+"_core",pdfname+"_core",*dt,*corebias,*coresigma,*dte);
  RooGaussModel *out = new RooGaussModel(pdfname+"_out",pdfname+"_out",*dt,*outbias,*outsigma);
RooAddModel *rf = new RooAddModel(pdfname+"_2GRF",pdfname+"_2GRF",RooArgList(*out,*core),RooArgList(*fout));
  return rf;	   
}


RooArgList MLPdfFactory::makeTottiParams(TString basename)
{
  RooRealVar *mean   = new RooRealVar(basename+"_mean",basename+"_mean",0.1);
  RooRealVar *sigma0 = new RooRealVar(basename+"_sigma0",basename+"_sigma0",0.1);
  RooRealVar *sigma1 = new RooRealVar(basename+"_sigma1",basename+"_sigma1",0.1);
  RooRealVar *alpha  = new RooRealVar(basename+"_alpha",basename+"_alpha",0.1);
  RooRealVar *N      = new RooRealVar(basename+"_N",basename+"_N",0.1);
  RooRealVar *f1     = new RooRealVar(basename+"_f1",basename+"_f1",0.1); 
  return RooArgList(*mean,*sigma0,*sigma1,*alpha,*N,*f1);
}

RooAbsPdf* MLPdfFactory::makeTottiPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters,
				      const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 6, args, 0)) return 0;
  RooRealVar *theobs = (RooRealVar*)obs.at(0);
  RooRealVar *mean   = (RooRealVar*)parameters.at(0);
  RooRealVar *sigma0 = (RooRealVar*)parameters.at(1);
  RooRealVar *sigma1 = (RooRealVar*)parameters.at(2);
  RooRealVar *alpha  = (RooRealVar*)parameters.at(3);
  RooRealVar *N      = (RooRealVar*)parameters.at(4);
  RooRealVar *f1     = (RooRealVar*)parameters.at(5);
  
  RooAbsPdf *gaussian = new RooGaussian(pdfname+"_gaussian",pdfname+"_gaussian",*theobs,*mean,*sigma0);
  RooAbsPdf *cryball = new RooCBShape(pdfname+"_cryball",pdfname+"_cryball",*theobs,*mean,*sigma1,*alpha,*N);
  RooAbsPdf *pdf = new RooAddPdf(pdfname, pdfname, RooArgList(*gaussian,*cryball), *f1);
  return pdf;
}


RooArgList MLPdfFactory::makeCruijffParams(TString basename)
{
  RooRealVar *mean   = new RooRealVar(basename+"_mean",basename+"_mean",0.1);
  RooRealVar *sigmaL = new RooRealVar(basename+"_sigmaL",basename+"_sigmaL",0.1);
  RooRealVar *sigmaR = new RooRealVar(basename+"_sigmaR",basename+"_sigmaR",0.1);
  RooRealVar *alphaL = new RooRealVar(basename+"_alphaL",basename+"_alphaL",0.1);
  RooRealVar *alphaR = new RooRealVar(basename+"_alphaR",basename+"_alphaR",0.1);
  return RooArgList(*mean,*sigmaL,*sigmaR,*alphaL,*alphaR);
}

RooAbsPdf* MLPdfFactory::makeCruijffPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters,
					const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 5, args, 0)) return 0;
  return new RhhCruijffPdf(pdfname, pdfname, *(RooAbsReal*)obs.at(0), 
			   *(RooAbsReal*)parameters.at(0), 
			   *(RooAbsReal*)parameters.at(1), 
			   *(RooAbsReal*)parameters.at(2),
			   *(RooAbsReal*)parameters.at(3),
			   *(RooAbsReal*)parameters.at(4));
}

RooArgList MLPdfFactory::makeDoubleCruijffParams(TString basename)
{
  RooRealVar *mean1   = new RooRealVar(basename+"_mean1", basename+"_mean1", 0);
  RooRealVar *mean2   = new RooRealVar(basename+"_mean2", basename+"_mean2", 0);
  RooRealVar *sigmaL1 = new RooRealVar(basename+"_sigmaL1", basename+"sigmaL1", 1);
  RooRealVar *sigmaR1 = new RooRealVar(basename+"_sigmaR1", basename+"sigmaR1", 1);
  RooRealVar *sigmaL2 = new RooRealVar(basename+"_sigmaL2", basename+"sigmaL2", 1);
  RooRealVar *sigmaR2 = new RooRealVar(basename+"_sigmaR2", basename+"sigmaR2", 1);
  RooRealVar *alphaL1 = new RooRealVar(basename+"_alphaL1", basename+"_alphaL1",0.1);
  RooRealVar *alphaR1 = new RooRealVar(basename+"_alphaR1", basename+"_alphaR1",0.1);
  RooRealVar *alphaL2 = new RooRealVar(basename+"_alphaL2", basename+"_alphaL2",0.1);
  RooRealVar *alphaR2 = new RooRealVar(basename+"_alphaR2", basename+"_alphaR2",0.1);
  RooRealVar *f1      = new RooRealVar(basename+"_f1", basename+"_f1",0.5);
  
  RooArgList set1(*mean1,*sigmaL1,*sigmaR1,*alphaL1,*alphaR1);
  RooArgList set2(*mean2,*sigmaL2,*sigmaR2,*alphaL2,*alphaR2,*f1); 

  set1.add(set2);
  return set1;
}
   
		    

RooAbsPdf* MLPdfFactory::makeDoubleCruijffPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 11, args, 0)) return 0;

   RhhCruijffPdf* pdf1 = new RhhCruijffPdf(pdfname+"_1", pdfname+"_1", *(RooAbsReal*)obs.at(0), 
					  *(RooAbsReal*)parameters.at(0), 
					  *(RooAbsReal*)parameters.at(1), 
					  *(RooAbsReal*)parameters.at(2),
					  *(RooAbsReal*)parameters.at(3),
					  *(RooAbsReal*)parameters.at(4) );

   RhhCruijffPdf* pdf2 = new RhhCruijffPdf(pdfname+"_2", pdfname+"_2", *(RooAbsReal*)obs.at(0), 
					  *(RooAbsReal*)parameters.at(5), 
					  *(RooAbsReal*)parameters.at(6), 
					  *(RooAbsReal*)parameters.at(7),
					  *(RooAbsReal*)parameters.at(8),
					  *(RooAbsReal*)parameters.at(9) );

   return new RooAddPdf( pdfname, pdfname, RooArgList(*pdf1, *pdf2), 
			*(RooAbsReal*)parameters.at(10) );
}

RooArgList MLPdfFactory::makeCrystalCruijffParams(TString basename)
{
  RooRealVar *mean    = new RooRealVar(basename+"_mean",basename+"_mean",0.1);
  RooRealVar *sigma   = new RooRealVar(basename+"_sigma",basename+"_sigma",0.1);
  RooRealVar *alpha   = new RooRealVar(basename+"_alpha",basename+"_alpha",0.1);
  RooRealVar *alphaCB = new RooRealVar(basename+"_alphaCB",basename+"_alphaCB",0.1);
  RooRealVar *nCB     = new RooRealVar(basename+"_nCB",basename+"_nCB",0.1);
  return RooArgList(*mean,*sigma,*alpha,*alphaCB,*nCB);
}

RooAbsPdf* MLPdfFactory::makeCrystalCruijffPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters,
                                               const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 5, args, 0)) return 0;
  return new RhhCrystalCruijffPdf(pdfname, pdfname, *(RooAbsReal*)obs.at(0),
                                  *(RooAbsReal*)parameters.at(0),
                                  *(RooAbsReal*)parameters.at(1),
                                  *(RooAbsReal*)parameters.at(2),
                                  *(RooAbsReal*)parameters.at(3),
                                  *(RooAbsReal*)parameters.at(4));
}

RooArgList MLPdfFactory::makeStepFunctionParams(TString basename, const TList &args)
{
  RooArgList parameters ;
  TObject* obj = args.At(0);
  if(!(obj->IsA()->InheritsFrom(TAxis::Class()))) {
    cout << "MLPdfFactory Error making parameters for StepFunction pdf" << endl ;
  } else {
    TAxis* axis = static_cast<TAxis*>(obj) ;
    int nbins = axis->GetNbins() ;
    int npars = nbins-1 ;
    double dx = axis->GetXmax()-axis->GetXmin() ;
    char parname[256] ;
    for(int i=1; i<=npars; ++i) {
      sprintf(parname,"%s_bin%d",basename.Data(),i);
      double ddx = axis->GetBinWidth(i) ;
      RooRealVar* var = new RooRealVar(parname,parname,1./dx,0,1./ddx) ;
      parameters.add(*var) ;
    }
  }
  return parameters ;
}

RooAbsPdf* MLPdfFactory::makeStepFunctionPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters,
				       const TList &args)
{
  RooAbsPdf* pdf(0) ;
  TObject* obj = args.At(0);
   if(!(obj->IsA()->InheritsFrom(TAxis::Class()))) {
    cout << "MLPdfFactory Error making StepFunction pdf" << endl ;
  } else {
    TAxis* axis = static_cast<TAxis*>(obj) ;
    pdf = new RooParametricStepFunction(pdfname,pdfname,*(RooAbsReal*)obs.at(0),parameters,
					*const_cast<TArrayD*>(axis->GetXbins()),axis->GetNbins()) ;
  }
  return pdf ;
}

RooArgList MLPdfFactory::makeBinnedPdfParams(TString basename, const TList &args)
{
  RooArgList parameters ;
  TObject* obj = args.At(0);
  if(!(obj->IsA()->InheritsFrom(TAxis::Class()))) {
    cout << "MLPdfFactory Error making parameters for binned pdf" << endl ;
  } else {
    TAxis* axis = static_cast<TAxis*>(obj) ;
    int nbins = axis->GetNbins() ;
    int npars = nbins-1 ;
    char parname[256] ;
    for(int i=1; i<=npars; ++i) {
      sprintf(parname,"%s_bin%d",basename.Data(),i);
      RooRealVar* var = new RooRealVar(parname,parname,1./(npars+1),0,1) ;
      parameters.add(*var) ;
    }
  }
  return parameters ;
}

RooAbsPdf* MLPdfFactory::makeBinnedPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters,
				       const TList &args)
{
  RooAbsPdf* pdf(0) ;
  TObject* obj = args.At(0);
   if(!(obj->IsA()->InheritsFrom(TAxis::Class()))) {
    cout << "MLPdfFactory Error making binned pdf" << endl ;
  } else {
    TAxis* axis = static_cast<TAxis*>(obj) ;
    pdf = new RhhBinnedPdf(pdfname,pdfname,
			   *(RooAbsReal*)obs.at(0),parameters,*(axis->GetXbins())) ;
  }
  return pdf ;
}

RooArgList MLPdfFactory::makeBreitWignerParams(TString basename)
{
  RooRealVar *mean = new RooRealVar(basename+"_mean",basename+"_mean",0);
  RooRealVar *width = new RooRealVar(basename+"_width",basename+"_width",1);
  return RooArgList(*mean,*width);
}

RooAbsPdf* MLPdfFactory::makeBreitWignerPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 2, args, 0)) return 0;
  return new RooBreitWigner(pdfname, pdfname, *(RooAbsReal*)obs.at(0), *(RooAbsReal*)parameters.at(0), *(RooAbsReal*)parameters.at(1));
}

RooArgList MLPdfFactory::makeBreitWignerAndPolyParams(TString basename)
{
  RooRealVar *mean     = new RooRealVar(basename+"_mean",basename+"_mean",0);
  RooRealVar *width    = new RooRealVar(basename+"_width",basename+"_width",1);
  RooRealVar *p1       = new RooRealVar(basename+"_p1",basename+"_p1",1);
  // Fraction of Breit-Wigner in the full pdf
  RooRealVar *fraction = new RooRealVar(basename+"_frac",basename+"_frac",1);
  return RooArgList(*mean, *width, *p1, *fraction);
}

RooAbsPdf* MLPdfFactory::makeBreitWignerAndPolyPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;
   RooAbsPdf* BW = new RooBreitWigner(pdfname+"_bw", pdfname+"bw", 
				      *(RooAbsReal*)obs.at(0), 
				      *(RooAbsReal*)parameters.at(0), 
				      *(RooAbsReal*)parameters.at(1) );
   RooAbsPdf* poly = new RooPolynomial(pdfname+"_poly", pdfname+"poly", 
				       *(RooAbsReal*)obs.at(0), 
				       RooArgList( *(RooAbsReal*)parameters.at(2) ) );

   return new RooAddPdf(pdfname, pdfname, RooArgSet(*BW,*poly), RooArgSet( *(RooAbsReal*)parameters.at(3) ) );
}

RooArgList MLPdfFactory::makeKeysParams(TString basename, const TList &args)
{
  return RooArgList();
}
RooArgList MLPdfFactory::makeBreitWignerAndExpoParams(TString basename)
{
  RooRealVar *mean     = new RooRealVar(basename+"_mean",basename+"_mean",0);
  RooRealVar *width    = new RooRealVar(basename+"_width",basename+"_width",1);
  RooRealVar *exp       = new RooRealVar(basename+"_exp",basename+"_exp",1);
  // Fraction of Breit-Wigner in the full pdf
  RooRealVar *fraction = new RooRealVar(basename+"_frac",basename+"_frac",1);
  return RooArgList(*mean, *width, *exp, *fraction);
}

RooAbsPdf* MLPdfFactory::makeBreitWignerAndExpoPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;
   RooAbsPdf* BW = new RooBreitWigner(pdfname+"_bw", pdfname+"bw", 
				      *(RooAbsReal*)obs.at(0), 
				      *(RooAbsReal*)parameters.at(0), 
				      *(RooAbsReal*)parameters.at(1) );
   RooAbsPdf* expo = new RooExponential(pdfname+"_expo", pdfname+"_expo",
					*dynamic_cast<RooRealVar*>(obs.at(0)),
					*dynamic_cast<RooRealVar*>(parameters.at(2))) ;
   return new RooAddPdf(pdfname, pdfname, RooArgSet(*BW,*expo), RooArgSet( *(RooAbsReal*)parameters.at(3) ) );
}

RooArgList MLPdfFactory::make2DKeysParams(TString basename, const TList &args)
{
  return RooArgList();
}

RooAbsPdf* MLPdfFactory::makeKeysPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters,const TList &args)
{ 
  if (!checkArgLength(pdfname, obs, 1, parameters, 0, args, 2)) return 0;
  
  TObjString *file = (TObjString*)args.At(0);
  TObjString *nameinfile = (TObjString*)args.At(1);
  RooRealVar *theobs   = (RooRealVar*)obs.at(0);
  TFile *f = new TFile(file->String());
  TObject *obj = f->Get(nameinfile->String());
  RooDataSet* dataset(0);
  if (obj->InheritsFrom("RooDataSet"))
     dataset = (RooDataSet*)obj;
  else
     dataset = new RooDataSet("datasetKeys", "input dataset for keys pdf",(TTree*)obj,*theobs);
  RooAbsPdf *pdf = new RooKeysPdf(pdfname, pdfname, *theobs, *dataset );
  return pdf;
}

RooAbsPdf* MLPdfFactory::make2DKeysPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters,const TList &args)
{ 
  if (!checkArgLength(pdfname, obs, 2, parameters, 0, args, 2)) return 0;
  
  TObjString *file = (TObjString*)args.At(0);
  TObjString *nameinfile = (TObjString*)args.At(1);
  RooRealVar *theobs1   = (RooRealVar*)obs.at(0);
  RooRealVar *theobs2   = (RooRealVar*)obs.at(1);
  TFile *f = new TFile(file->String());
  TObject *obj = f->Get(nameinfile->String());
  RooDataSet* dataset(0);
  if (obj->InheritsFrom("RooDataSet"))
     dataset = (RooDataSet*)obj;
  else
     dataset = new RooDataSet("dataset2DKeys", "input dataset for 2D keys pdf",
			      (TTree*)obj,RooArgSet(RooArgList(*theobs1, *theobs2)));
  RooAbsPdf *pdf = new Roo2DKeysPdf(pdfname, pdfname, *theobs1, *theobs2 , *dataset );
  return pdf;
}

RooArgList MLPdfFactory::makeHistParams(TString basename, const TList &args)
{
  return RooArgList();
}

RooAbsPdf*  
MLPdfFactory::makeHistPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  assert( args.GetSize() >=1 ) ;
  TFile* f(0) ;
  RooDataHist* ds(0) ;
  bool isowner(false) ;
  TObject* obj = args.At(0) ;
  if( obj->InheritsFrom("RooDataHist") ) {
    ds = static_cast<RooDataHist*>(obj) ;
  } else if(obj->InheritsFrom("TObjString") ) {
    assert(  args.GetSize() >=2 ) ;
    TObjString *file = dynamic_cast<TObjString*>(obj);
    TObjString *nameinfile = dynamic_cast<TObjString*>(args.At(1));
    if(file && nameinfile) {
      f = new TFile(file->String());
      if(f->IsOpen()) {
	obj = f->Get(nameinfile->String());
	if( obj ) {
	  if( obj->InheritsFrom("RooDataHist") ) 
	    ds = dynamic_cast<RooDataHist*>(obj) ;
	  else if( obj->InheritsFrom("RooAbsData") ) {
	    RooAbsData* origset = dynamic_cast<RooAbsData*>(obj) ;
	    // apply a selection if argument given
	    TString selection="" ;
	    if( args.GetSize()>=3 && (obj = args.At(2))->InheritsFrom("TObjString")) 
	      selection = dynamic_cast<TObjString*>(obj)->String() ;
	    // reduce to correct range
	    RooAbsData* reducedset = origset->reduce(obs,selection) ;
	    cout << " dataset from size " << origset->numEntries() << " to size " 
		 << reducedset->numEntries() << "." << endl ;
	    // now project into RooDataHist
	    ds = new RooDataHist("dstmp","dstmp",obs,*reducedset) ;
	    isowner = true ;
	    delete reducedset ;
	  }
	} else cout << "cannot find object with name " << nameinfile << " in file " << file << endl ;
      } else cout << "cannot open file " << file << endl ;
    }
  } 
  RooHistPdf* pdf(0) ;
  const int interpolationorder=0 ;
  if(ds) pdf = new RooHistPdf(pdfname,pdfname,obs,*ds,interpolationorder) ;
  //if(isowner) delete ds ;
  //if(f) delete f ;
  return pdf ;
}

RooArgList MLPdfFactory::makeExpoParams(TString basename, const TList &args)
{
  RooRealVar *exp = new RooRealVar(basename+"_exp",basename+"_exp",0);
  return RooArgList(*exp);
}

RooAbsPdf*  
MLPdfFactory::makeExpoPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 1, args, 0)) return 0;
  return new RooExponential(pdfname, pdfname,
			    *dynamic_cast<RooRealVar*>(obs.at(0)),
			    *dynamic_cast<RooRealVar*>(parameters.at(0))) ;
}

RooArgList MLPdfFactory::makeModExpParams(TString basename, const TList &args)
{
  RooRealVar *slope = new RooRealVar(basename+"_slope",basename+"_slope",-0.1);
  RooRealVar *shape = new RooRealVar(basename+"_shape",basename+"_shape",0.01);
  return RooArgList(*slope,*shape);
}

RooAbsPdf*
MLPdfFactory::makeModExpPdf(TString  pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 2, args, 0)) return 0;
  return new RhhModExpPdf(pdfname, pdfname, *(RooAbsReal*)obs.at(0),
                          *(RooAbsReal*)parameters.at(0),
                          *(RooAbsReal*)parameters.at(1));
}



// DIJET

RooArgList MLPdfFactory::makeNonResParams(TString basename)
{
  RooRealVar *sqrtS = new RooRealVar(basename+"_sqrtS", basename+"_sqrtS", 7000.);
  RooRealVar *p0    = new RooRealVar(basename+"_p0", basename+"_p0", 0.);
  RooRealVar *p1    = new RooRealVar(basename+"_p1", basename+"_p1", 1.);
  RooRealVar *p2    = new RooRealVar(basename+"_p2", basename+"_p2", 1.);
  RooRealVar *p3    = new RooRealVar(basename+"_p3", basename+"_p3", 1.);

  return RooArgList(*sqrtS, *p0, *p1, *p2, *p3);
}

RooAbsPdf* MLPdfFactory::makeNonResPdf(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 5, args, 0)) return 0;
  return new MLFitNonResPdf(pdfname, pdfname, *(RooAbsReal*)obs.at(0), 
			    *(RooAbsReal*)parameters.at(0), 
			    *(RooAbsReal*)parameters.at(1), 
			    *(RooAbsReal*)parameters.at(2), 
			    *(RooAbsReal*)parameters.at(3), 
			    *(RooAbsReal*)parameters.at(4)); 
}

RooArgList MLPdfFactory::makeNonResParamsSYSA(TString basename)
{
  RooRealVar *sqrtS = new RooRealVar(basename+"_sqrtS", basename+"_sqrtS", 7000.);
  RooRealVar *p0    = new RooRealVar(basename+"_p0", basename+"_p0", 0.);
  RooRealVar *p1    = new RooRealVar(basename+"_p1", basename+"_p1", 1.);
  RooRealVar *p2    = new RooRealVar(basename+"_p2", basename+"_p2", 1.);

  return RooArgList(*sqrtS, *p0, *p1, *p2);
}

RooAbsPdf* MLPdfFactory::makeNonResPdfSYSA(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;
  return new MLFitNonResSYSAPdf(pdfname, pdfname, *(RooAbsReal*)obs.at(0), 
				*(RooAbsReal*)parameters.at(0), 
				*(RooAbsReal*)parameters.at(1), 
				*(RooAbsReal*)parameters.at(2), 
				*(RooAbsReal*)parameters.at(3)); 
}

RooArgList MLPdfFactory::makeNonResParamsSYSB(TString basename)
{
  RooRealVar *sqrtS = new RooRealVar(basename+"_sqrtS", basename+"_sqrtS", 7000.);
  RooRealVar *p0    = new RooRealVar(basename+"_p0", basename+"_p0", 0.);
  RooRealVar *p1    = new RooRealVar(basename+"_p1", basename+"_p1", 1.);
  RooRealVar *p2    = new RooRealVar(basename+"_p2", basename+"_p2", 1.);

  return RooArgList(*sqrtS, *p0, *p1, *p2);
}

RooAbsPdf* MLPdfFactory::makeNonResPdfSYSB(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;
  return new MLFitNonResSYSBPdf(pdfname, pdfname, *(RooAbsReal*)obs.at(0), 
				   *(RooAbsReal*)parameters.at(0), 
				   *(RooAbsReal*)parameters.at(1), 
				   *(RooAbsReal*)parameters.at(2),
				   *(RooAbsReal*)parameters.at(3)); 
}

RooArgList MLPdfFactory::makeNonResParamsSYSC(TString basename)
{
  RooRealVar *sqrtS = new RooRealVar(basename+"_sqrtS", basename+"_sqrtS", 7000.);
  RooRealVar *p0    = new RooRealVar(basename+"_p0", basename+"_p0", 0.);
  RooRealVar *p1    = new RooRealVar(basename+"_p1", basename+"_p1", 1.);
  RooRealVar *p2    = new RooRealVar(basename+"_p2", basename+"_p2", 1.);

  return RooArgList(*sqrtS, *p0, *p1, *p2);
}

RooAbsPdf* MLPdfFactory::makeNonResPdfSYSC(TString pdfname, const RooArgList &obs, RooArgList parameters, const TList &args)
{
  if (!checkArgLength(pdfname, obs, 1, parameters, 4, args, 0)) return 0;
  return new MLFitNonResSYSCPdf(pdfname, pdfname, *(RooAbsReal*)obs.at(0),
                                   *(RooAbsReal*)parameters.at(0),
                                   *(RooAbsReal*)parameters.at(1),
                                   *(RooAbsReal*)parameters.at(2),
                                   *(RooAbsReal*)parameters.at(3));
}
