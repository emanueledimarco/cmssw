#ifndef CondFormats_EcalObjects_EcalCubicPulseShapes_h
#define CondFormats_EcalObjects_EcalCubicPulseShapes_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

template <class P>
struct EcalCubicPulseShapeT {
public:
  static constexpr int TEMPLATESAMPLES = static_cast<int>(P::kPulseShapeTemplateSampleSize);
  static constexpr int PARSPERSAMPLE = static_cast<int>(P::kParsPerTemplateSample);

  EcalCubicPulseShapeT();

  float parameters[TEMPLATESAMPLES*PARSPERSAMPLE];

  float pdfval(int ipar) const { 
    int baseIndex = (ipar / PARSPERSAMPLE) * PARSPERSAMPLE;
    return parameters[baseIndex];
  }

  const float* splinepars(int ipar) const {
    int baseIndex = (ipar / PARSPERSAMPLE) * PARSPERSAMPLE;
    return &parameters[baseIndex + 1];
  }

  COND_SERIALIZABLE;
};



using EcalPh1CubicPulseShape = EcalCubicPulseShapeT<ecalPh1>;

typedef EcalCondObjectContainer<EcalPh1CubicPulseShape> EcalPh1CubicPulseShapesMap;
typedef EcalPh1CubicPulseShapesMap::const_iterator EcalPh1CubicPulseShapesMapIterator;
typedef EcalPh1CubicPulseShapesMap EcalPh1CubicPulseShapes;



using EcalPh2CubicPulseShape = EcalCubicPulseShapeT<ecalPh2>;

typedef EcalCondObjectContainer<EcalPh2CubicPulseShape> EcalPh2CubicPulseShapesMap;
typedef EcalPh2CubicPulseShapesMap::const_iterator EcalPh2CubicPulseShapesMapIterator;
typedef EcalPh2CubicPulseShapesMap EcalPh2CubicPulseShapes;

#endif
