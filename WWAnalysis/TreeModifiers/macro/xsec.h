#ifndef XSEC_H
#define XSEC_H

#include <map>
#include "WWAnalysis/TreeModifiers/interface/XSecProvider.h"
XSecProvider xsecProvider;

float evt_7TeV(int id){
  return xsecProvider.getEvt7TeV(id);
}

float evt_8TeV(int id){
  return xsecProvider.getEvt8TeV(id);
}

float get7TeVggHXsec(float mass) {
  return xsecProvider.get7TeVggHXsec(mass);
}    

float get8TeVggHXsec(float mass) {
  return xsecProvider.get8TeVggHXsec(mass);
}    

float get7TeVVBFXsec(float mass) {
  return xsecProvider.get7TeVVBFXsec(mass);
}    

float get8TeVVBFXsec(float mass) {
  return xsecProvider.get8TeVVBFXsec(mass);
}    

float get7TeVWHiXsec(float mass) {
  return xsecProvider.get7TeVWHiXsec(mass);
}

float get8TeVWHiXsec(float mass) {
  return xsecProvider.get8TeVWHiXsec(mass);
}

float get7TeVZHiXsec(float mass) {
  return xsecProvider.get7TeVZHiXsec(mass);
}

float get8TeVZHiXsec(float mass) {
  return xsecProvider.get8TeVZHiXsec(mass);
}

float get7TeVttHXsec(float mass) {
  return xsecProvider.get7TeVttHXsec(mass);
}

float get8TeVttHXsec(float mass) {
  return xsecProvider.get8TeVttHXsec(mass);
}

float get7TeVggHXsecByChannel(float mass, int ch) {
  return xsecProvider.get7TeVggHXsecByChannel(mass,ch);
}

float get8TeVggHXsecByChannel(float mass, int ch) {
  return xsecProvider.get8TeVggHXsecByChannel(mass,ch);
}

float get7TeVVBFXsecByChannel(float mass, int ch) {
  return xsecProvider.get7TeVVBFXsecByChannel(mass,ch);
}

float get8TeVVBFXsecByChannel(float mass, int ch) {
  return xsecProvider.get8TeVVBFXsecByChannel(mass,ch);
}

float get7TeVWHiXsecByChannel(float mass, int ch) {
  return xsecProvider.get7TeVWHiXsecByChannel(mass,ch);
}

float get8TeVWHiXsecByChannel(float mass, int ch) {
  return xsecProvider.get8TeVWHiXsecByChannel(mass,ch);
}

float get7TeVZHiXsecByChannel(float mass, int ch) {
  return xsecProvider.get7TeVZHiXsecByChannel(mass,ch);
}

float get8TeVZHiXsecByChannel(float mass, int ch) {
  return xsecProvider.get8TeVZHiXsecByChannel(mass,ch);
}

float get7TeVttHXsecByChannel(float mass, int ch) {
  return xsecProvider.get7TeVttHXsecByChannel(mass,ch);
}

float get8TeVttHXsecByChannel(float mass, int ch) {
  return xsecProvider.get8TeVttHXsecByChannel(mass,ch);
}

float get7TeVggHQCDUp(float mass) {
  return xsecProvider.get7TeVggHQCDUp(mass);
}    

float get7TeVggHQCDDown(float mass) {
  return xsecProvider.get7TeVggHQCDDown(mass);
}    

float get7TeVVBFQCDUp(float mass) {
  return xsecProvider.get7TeVVBFQCDUp(mass);
}   

float get7TeVVBFQCDDown(float mass) {
  return xsecProvider.get7TeVVBFQCDDown(mass);
}   

float get7TeVWHiQCDUp(float mass) {
  return xsecProvider.get7TeVWHiQCDUp(mass);
}

float get7TeVWHiQCDDown(float mass) {
  return xsecProvider.get7TeVWHiQCDDown(mass);
}

float get7TeVZHiQCDUp(float mass) {
  return xsecProvider.get7TeVZHiQCDUp(mass);
}

float get7TeVZHiQCDDown(float mass) {
  return xsecProvider.get7TeVZHiQCDDown(mass);
}

float get7TeVttHQCDUp(float mass) {
  return xsecProvider.get7TeVttHQCDUp(mass);
}

float get7TeVttHQCDDown(float mass) {
  return xsecProvider.get7TeVttHQCDDown(mass);
}

float get8TeVggHQCDUp(float mass) {
  return xsecProvider.get8TeVggHQCDUp(mass);
}    

float get8TeVggHQCDDown(float mass) {
  return xsecProvider.get8TeVggHQCDDown(mass);
}    

float get8TeVVBFQCDUp(float mass) {
  return xsecProvider.get8TeVVBFQCDUp(mass);
}   

float get8TeVVBFQCDDown(float mass) {
  return xsecProvider.get8TeVVBFQCDDown(mass);
}   

float get8TeVWHiQCDUp(float mass) {
  return xsecProvider.get8TeVWHiQCDUp(mass);
}

float get8TeVWHiQCDDown(float mass) {
  return xsecProvider.get8TeVWHiQCDDown(mass);
}

float get8TeVZHiQCDUp(float mass) {
  return xsecProvider.get8TeVZHiQCDUp(mass);
}

float get8TeVZHiQCDDown(float mass) {
  return xsecProvider.get8TeVZHiQCDDown(mass);
}

float get8TeVttHQCDUp(float mass) {
  return xsecProvider.get8TeVttHQCDUp(mass);
}

float get8TeVttHQCDDown(float mass) {
  return xsecProvider.get8TeVttHQCDDown(mass);
}

float get7TeVggHPDFUp(float mass) {
  return xsecProvider.get7TeVggHPDFUp(mass);
}   

float get7TeVggHPDFDown(float mass) {
  return xsecProvider.get7TeVggHPDFDown(mass);
}   

float get7TeVVBFPDFUp(float mass) {
  return xsecProvider.get7TeVVBFPDFUp(mass);
}   

float get7TeVVBFPDFDown(float mass) {
  return xsecProvider.get7TeVVBFPDFDown(mass);
}   

float get7TeVWHiPDFUp(float mass) {
  return xsecProvider.get7TeVWHiPDFUp(mass);
}

float get7TeVWHiPDFDown(float mass) {
  return xsecProvider.get7TeVWHiPDFDown(mass);
}

float get7TeVZHiPDFUp(float mass) {
  return xsecProvider.get7TeVZHiPDFUp(mass);
}

float get7TeVZHiPDFDown(float mass) {
  return xsecProvider.get7TeVZHiPDFDown(mass);
}

float get7TeVttHPDFUp(float mass) {
  return xsecProvider.get7TeVttHPDFUp(mass);
}

float get7TeVttHPDFDown(float mass) {
  return xsecProvider.get7TeVttHPDFDown(mass);
}

float get8TeVggHPDFUp(float mass) {
  return xsecProvider.get8TeVggHPDFUp(mass);
}

float get8TeVggHPDFDown(float mass) {
  return xsecProvider.get8TeVggHPDFDown(mass);
}

float get8TeVVBFPDFUp(float mass) {
  return xsecProvider.get8TeVVBFPDFUp(mass);
}    

float get8TeVVBFPDFDown(float mass) {
    return xsecProvider.get8TeVVBFPDFDown(mass);
}

float get8TeVWHiPDFUp(float mass) {
  return xsecProvider.get8TeVWHiPDFUp(mass);
}

float get8TeVWHiPDFDown(float mass) {
    return xsecProvider.get8TeVWHiPDFDown(mass);
}

float get8TeVZHiPDFUp(float mass) {
  return xsecProvider.get8TeVZHiPDFUp(mass);
}

float get8TeVZHiPDFDown(float mass) {
    return xsecProvider.get8TeVZHiPDFDown(mass);
}

float get8TeVttHPDFUp(float mass) {
  return xsecProvider.get8TeVttHPDFUp(mass);
}

float get8TeVttHPDFDown(float mass) {
    return xsecProvider.get8TeVttHPDFDown(mass);
}

float get7TeVBkgXsec(int id) {
  return xsecProvider.get7TeVBkgXsec(id);
}

float get8TeVBkgXsec(int id) {
  return xsecProvider.get8TeVBkgXsec(id);
}



#endif
