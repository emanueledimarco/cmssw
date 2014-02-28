#include "TROOT.h"
#include "TH1D.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

class PulseShapeModel {
public:
  PulseShapeModel() : model(0) { }
  virtual ~PulseShapeModel() { delete model; }
  RooAbsPdf *model;
};

class ECALShapeConvGaussian : public PulseShapeModel {
public:
  ECALShapeConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, RooRealVar *mean0 = 0, RooRealVar *sigma0=0);
  ~ECALShapeConvGaussian();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
};

//--------------------------------------------------------------------------------------------------
CMCTemplateConvGaussian::CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, 
                                                 RooRealVar *mean0, RooRealVar *sigma0)
{  
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];  
  
  if (pass) {
    if(mean0)  { mean  = mean0;  }
    else       { sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10); }
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5);    }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  } else {
    if(mean0)  { mean  = mean0;  }
    else       { sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10); }
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0,2.5);    }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  }

  sprintf(vname,"inHist_%s",hist->GetName());
  inHist = (TH1D*)hist->Clone(vname);
  
  sprintf(vname,"dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
  sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,8);
//     sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,10); //use this if fit doesn't converge well...
//      sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,4); //use this if fit doesn't converge well...
//    sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,1); //use this if fit doesn't converge well...
//      sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,0); //use this if fit doesn't converge well...
  sprintf(vname,"signal%s",name);   model    = new RooFFTConvPdf(vname,vname,m,*histPdf,*gaus);
}

CMCTemplateConvGaussian::~CMCTemplateConvGaussian()
{
  delete mean;
  //delete sigma;
  delete gaus;
  delete inHist;
  delete dataHist;
  delete histPdf;
}

  
