#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"

EcalPulseCovariance::EcalPulseCovariance() {
  int N = EcalPulseShape::TEMPLATESAMPLES*(EcalPulseShape::TEMPLATESAMPLES+1)/2;
  for(int k=0; k<N; ++k) covval[k] = 0.;
}
