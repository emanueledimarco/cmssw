#ifndef CondFormats_EcalObjects_EcalCubicPulseShapes_h
#define CondFormats_EcalObjects_EcalCubicPulseShapes_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

template <class P>
struct EcalCubicPulseShapeT {
public:
  static constexpr int TEMPLATESAMPLES = static_cast<int>(P::kPulseShapeTemplateSampleSize);

  EcalCubicPulseShapeT();

  float pdfval[TEMPLATESAMPLES];

  float val(int isample) const { return pdfval[isample]; }

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
