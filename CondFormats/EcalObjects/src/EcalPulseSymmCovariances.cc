#include "CondFormats/EcalObjects/interface/EcalPulseSymmCovariances.h"

EcalPulseSymmCovariance::EcalPulseSymmCovariance() {
  int N = TEMPLATESAMPLES * (TEMPLATESAMPLES + 1) / 2;
  for (int k = 0; k < N; ++k)
    covval[k] = 0.;
}

EcalPh2PulseSymmCovariance::EcalPh2PulseSymmCovariance() {
  int N = TEMPLATESAMPLES * (TEMPLATESAMPLES + 1) / 2;
  for (int k = 0; k < N; ++k)
    covval[k] = 0.;
}
