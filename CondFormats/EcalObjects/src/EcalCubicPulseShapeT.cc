#include "CondFormats/EcalObjects/interface/EcalCubicPulseShapeT.h"

template <class P>
EcalCubicPulseShapeT<P>::EcalCubicPulseShapeT() {
  for (int s = 0; s < TEMPLATESAMPLES; ++s)
    pdfval[s] = 0.;
}

template struct EcalCubicPulseShapeT<ecalPh2>;
template struct EcalCubicPulseShapeT<ecalPh1>;
