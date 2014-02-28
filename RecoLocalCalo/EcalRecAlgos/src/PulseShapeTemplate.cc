#include "RecoLocalCalo/EcalRecAlgos/interface/PulseShapeTemplate.cc"

//--------------------------------------------------------------------------------------------------
ECALShapeConvGaussian::ECALShapeConvGaussian(RooRealVar &m, TH1D* hist,
					     RooRealVar *mean0, RooRealVar *sigma0) {  

  if(mean0)  { mean  = mean0;  }
  else       { mean  = new RooRealVar("mean","mean",0,-10,10); }
  if(sigma0) { sigma = sigma0; }
  else       { sigma = new RooRealVar("sigma","sigma",1);    }
  gaus  = new RooGaussian("gaus","gaus",m,*mean,*sigma);

  inHist = (TH1D*)hist->Clone("inHist");
  
  dataHist = new RooDataHist("dataHist","dataHist",RooArgSet(m),inHist);
  histPdf  = new RooHistPdf("histPdf","histPdf",m,*dataHist,8);
  model    = new RooFFTConvPdf("model","model",m,*histPdf,*gaus);
}

ECALShapeConvGaussian::~ECALShapeConvGaussian()
{
  delete mean;
  delete sigma;
  delete gaus;
  delete inHist;
  delete dataHist;
  delete histPdf;
}

  
