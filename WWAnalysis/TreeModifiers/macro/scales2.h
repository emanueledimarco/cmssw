#ifndef SCALES_H
#define SCALES_H

#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include "pu.h"

#include "WWAnalysis/TreeModifiers/interface/LeptSfProvider.h"
#include "WWAnalysis/TreeModifiers/interface/XSecProvider.h"
#include "WWAnalysis/TreeModifiers/interface/MassErrCorrProvider.h"

/////////////////////// Lepton efficiency scale factors /////////////////////////////////////////////

bool is2011;
bool ismc;
bool initialized;

LeptSfProvider sfProvider;
XSecProvider xsecProvider;
MassErrCorrProvider ebeProvider;

float getMuonIDSF(float pt, float eta) {
  return sfProvider.getMuonIDSF(pt,eta);
}

float getMuonIPSF(float pt, float eta) {
  return sfProvider.getMuonIDSF(pt,eta);
}

float getMuonIsoSF(float pt, float eta) {
  return  sfProvider.getMuonIsoSF(pt,eta);
}

float getMuonSF(float pt, float eta) {
  return sfProvider.getMuonSF(pt, eta);
}

float getElectronSF(float pt, float eta) {
  return sfProvider.getElectronSF(pt,eta);
}

float getPR(float pt, float eta, float id) {
  return  sfProvider.getPR(pt,eta,id);
}

float getSF(float pt, float eta, float id) {
  return  sfProvider.getSF(pt,eta,id);
}

///////////// Cross-sections and related uncertainties /////////////////////////////

float evt_8TeV(int id) {
    return xsecProvider.getEvt8TeV(id);
}

float evt_7TeV(int id) {
    return xsecProvider.getEvt7TeV(id);
}


float getBkgXsec(int id) {
    return is2011 ? xsecProvider.get7TeVBkgXsec(id)   : xsecProvider.get8TeVBkgXsec(id);
}

float getXsecggH(float mass) {
    return is2011 ? xsecProvider.get7TeVggHXsec(mass) : xsecProvider.get8TeVggHXsec(mass);
}

float getXsecVBF(float mass) {
    return is2011 ? xsecProvider.get7TeVVBFXsec(mass) : xsecProvider.get8TeVVBFXsec(mass);
}

float getXsecWHi(float mass) {
    return is2011 ? xsecProvider.get7TeVWHiXsec(mass) : xsecProvider.get8TeVWHiXsec(mass);
}

float getXsecZHi(float mass) {
    return is2011 ? xsecProvider.get7TeVZHiXsec(mass) : xsecProvider.get8TeVZHiXsec(mass);
}

float getXsecttH(float mass) {
    return is2011 ? xsecProvider.get7TeVttHXsec(mass) : xsecProvider.get8TeVttHXsec(mass);
}

float getXsecggHByChannel(float mass, int ch, bool separate2e2mu=false) {
    return is2011 ? xsecProvider.get7TeVggHXsecByChannel(mass, ch, separate2e2mu) : xsecProvider.get8TeVggHXsecByChannel(mass, ch, separate2e2mu);
}

float getXsecVBFByChannel(float mass, int ch, bool separate2e2mu=false) {
    return is2011 ? xsecProvider.get7TeVVBFXsecByChannel(mass, ch, separate2e2mu) : xsecProvider.get8TeVVBFXsecByChannel(mass, ch, separate2e2mu);
}

float getXsecWHiByChannel(float mass, int ch, bool separate2e2mu=false) {
    return is2011 ? xsecProvider.get7TeVWHiXsecByChannel(mass, ch, separate2e2mu) : xsecProvider.get8TeVWHiXsecByChannel(mass, ch, separate2e2mu);
}

float getXsecZHiByChannel(float mass, int ch, bool separate2e2mu=false) {
    return is2011 ? xsecProvider.get7TeVZHiXsecByChannel(mass, ch, separate2e2mu) : xsecProvider.get8TeVZHiXsecByChannel(mass, ch, separate2e2mu);
}

float getXsecttHByChannel(float mass, int ch, bool separate2e2mu=false) {
    return is2011 ? xsecProvider.get7TeVttHXsecByChannel(mass, ch, separate2e2mu) : xsecProvider.get8TeVttHXsecByChannel(mass, ch, separate2e2mu);
}

float getggHQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVggHQCDUp(mass)   : xsecProvider.get8TeVggHQCDUp(mass);
    else    return is2011 ? xsecProvider.get7TeVggHQCDDown(mass) : xsecProvider.get8TeVggHQCDDown(mass);

}

float getVBFQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVVBFQCDUp(mass)   : xsecProvider.get8TeVVBFQCDUp(mass);
    else    return is2011 ? xsecProvider.get7TeVVBFQCDDown(mass) : xsecProvider.get8TeVVBFQCDDown(mass);

}

float getWHiQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVWHiQCDUp(mass)   : xsecProvider.get8TeVWHiQCDUp(mass);
    else    return is2011 ? xsecProvider.get7TeVWHiQCDDown(mass) : xsecProvider.get8TeVWHiQCDDown(mass);

}

float getZHiQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVZHiQCDUp(mass)   : xsecProvider.get8TeVZHiQCDUp(mass);
    else    return is2011 ? xsecProvider.get7TeVZHiQCDDown(mass) : xsecProvider.get8TeVZHiQCDDown(mass);

}

float getttHQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVttHQCDUp(mass)   : xsecProvider.get8TeVttHQCDUp(mass);
    else    return is2011 ? xsecProvider.get7TeVttHQCDDown(mass) : xsecProvider.get8TeVttHQCDDown(mass);

}

float getggHPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVggHPDFUp(mass)   : xsecProvider.get8TeVggHPDFUp(mass);
    else    return is2011 ? xsecProvider.get7TeVggHPDFDown(mass) : xsecProvider.get8TeVggHPDFDown(mass);

}

float getVBFPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVVBFPDFUp(mass)   : xsecProvider.get8TeVVBFPDFUp(mass);
    else    return is2011 ? xsecProvider.get7TeVVBFPDFDown(mass) : xsecProvider.get8TeVVBFPDFDown(mass);

}

float getWHiPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVWHiPDFUp(mass)   : xsecProvider.get8TeVWHiPDFUp(mass);
    else    return is2011 ? xsecProvider.get7TeVWHiPDFDown(mass) : xsecProvider.get8TeVWHiPDFDown(mass);

}

float getZHiPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVZHiPDFUp(mass)   : xsecProvider.get8TeVZHiPDFUp(mass);
    else    return is2011 ? xsecProvider.get7TeVZHiPDFDown(mass) : xsecProvider.get8TeVZHiPDFDown(mass);

}

float getttHPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? xsecProvider.get7TeVttHPDFUp(mass)   : xsecProvider.get8TeVttHPDFUp(mass);
    else    return is2011 ? xsecProvider.get7TeVttHPDFDown(mass) : xsecProvider.get8TeVttHPDFDown(mass);

}

float getQQZZPDFUncertainty7TeV(float m4l) {
    return 1. + 0.0035*sqrt(m4l-30.);
}

float getGGZZPDFUncertainty7TeV(float m4l) {
    return 1. + 0.0066*sqrt(m4l-10.);
}

float getQQZZQCDScaleUncertainty7TeV(float m4l) {
    return 1. + 0.01*sqrt((m4l-20.)/13.);
}

float getGGZZQCDScaleUncertainty7TeV(float m4l) {
    return 1.04 + 0.1*sqrt((m4l+40.)/40.);
}



///////////// Pileup reweighting /////////////////////////////

float getPUWeight(float numsim, int mode=1) {
    return is2011 ? getPUWeight2011(numsim) : getPUWeight2012(numsim, mode);
}

/////////////////////// Event by event errors /////////////////////////////////////////////

float getMassErrCorr(float pt[4], float eta[4], float lid[4], float lmassErr[4], float phoMassErr[2]) {
  return ebeProvider.getMassErrCorr(pt,eta,lid,lmassErr,phoMassErr);
}

///////////// Initialize all global variables /////////////////////////////

void setmc(bool isMC) { ismc = isMC; }

void init(bool is7TeV) {

    if (initialized) {
        if (is7TeV == is2011) return;            
        else {
            is2011 = is7TeV;
	    sfProvider.init(is2011);
            ebeProvider.init(is2011,ismc);
        }
    }
    else {
        initialized = true;
        is2011 = is7TeV;
	sfProvider.init(is2011);
	xsecProvider.initXsec();
	xsecProvider.initQCDScale();
	xsecProvider.initPDF();
        initpuweights();
        ebeProvider.init(is2011,ismc);
    }

}

#endif
