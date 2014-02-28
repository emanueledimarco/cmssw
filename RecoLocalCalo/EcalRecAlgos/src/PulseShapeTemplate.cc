#include "RecoLocalCalo/EcalRecAlgos/interface/PulseShapeTemplate.hh"

//--------------------------------------------------------------------------------------------------
ECALShapeConvGaussian::ECALShapeConvGaussian(const char *name, RooRealVar &m, TH1D* hist,
					     RooRealVar *mean0, RooRealVar *sigma0) {  

  if(mean0)  { mean  = mean0;  }
  else       { mean  = new RooRealVar("mean","mean",0,-10,10); }
  if(sigma0) { sigma = sigma0; }
  else       { sigma = new RooRealVar("sigma","sigma",1);    }
  
  char buf[200];
  sprintf(buf,"gaus_%s",name);
  gaus  = new RooGaussian(buf,buf,m,*mean,*sigma);

  sprintf(buf,"inHist_%s",name);
  inHist = (TH1D*)hist->Clone(buf);
  
  sprintf(buf,"dataHist_%s",name);
  dataHist = new RooDataHist(buf,buf,RooArgSet(m),inHist);
  sprintf(buf,"histPdf_%s",name);
  histPdf  = new RooHistPdf(buf,buf,m,*dataHist,8);
  sprintf(buf,"model_%s",name);
  model    = new RooFFTConvPdf(buf,buf,m,*histPdf,*gaus);
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

  
