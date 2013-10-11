// $Id: HiggsToGammaGammaLinkDef.h,v 1.1 2011/10/22 10:49:11 meridian Exp $

#ifndef HIGGSANALYSIS_LINKDEF_H
#define HIGGSANALYSIS_LINKDEF_H


#include "../interface/GBRTree.h"
#include "../interface/GBRForest.h"



#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class GBRTree+; 
#pragma link C++ class GBRForest+; 

#endif
