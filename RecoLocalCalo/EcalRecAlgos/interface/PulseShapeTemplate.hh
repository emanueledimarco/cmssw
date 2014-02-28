#include "TROOT.h"
#include "TH1D.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"

class PulseShapeModel {
public:
  PulseShapeModel() : model(0) { }
  virtual ~PulseShapeModel() { delete model; }
  RooAbsPdf *model;
};

class ECALShapeConvGaussian : public PulseShapeModel {
public:
  ECALShapeConvGaussian(const char *name, RooRealVar &m, TH1D* hist, RooRealVar *mean0 = 0, RooRealVar *sigma0=0);
  ~ECALShapeConvGaussian();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
};

