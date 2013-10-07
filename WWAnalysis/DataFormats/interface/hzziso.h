#ifndef HZZISO_H
#define HZZISO_H

struct HZZIso {

    static float getEcalIsoCustom(float isoval, float rho, float pt, float abseta, bool isMu, bool isMC) {
        float EAmuB = !isMC?0.090:0.089; 
        float EAmuE = !isMC?0.045:0.054; 
        float EAelB = !isMC?0.093:0.117; 
        float EAelE = !isMC?0.022:0.030;

        float elcorr = !isMC?0.124:0.073;
        float mucorr = !isMC?0.127:0.076;

        float iso = isoval;

        if (isMu) iso -= rho*(abseta < 1.479 ? EAmuB : EAmuE);

        else iso -= rho*(abseta < 1.479 ? EAelB : EAelE);

        if (isMu  && abseta < 1.479) iso -= mucorr;
        if (!isMu && abseta < 1.479) iso -= elcorr;

        return (iso>0.0 ? iso : 0.0);
    }

    static float getHcalIsoCustom(float isoval, float rho, float pt, float abseta, bool isMu, bool isMC) {
        float EAmuB = !isMC?0.035:0.036; 
        float EAmuE = !isMC?0.058:0.064; 
        float EAelB = !isMC?0.061:0.067; 
        float EAelE = !isMC?0.160:0.174;

        float iso = isoval;

        if (isMu) iso -= rho*(abseta < 1.479 ? EAmuB : EAmuE);
        else iso -= rho*(abseta < 1.479 ? EAelB : EAelE);

        return (iso>0.0 ? iso : 0.0);
    }


    static float getEcalIso(float isoval, float rho, float pt, float abseta, bool isMu) {
        float EAmuB = 0.074;
        float EAmuE = 0.045;
        float EAelB = 0.101;
        float EAelE = 0.046;

        float iso = isoval;

        if (isMu) iso -= rho*(abseta < 1.479 ? EAmuB : EAmuE);
        else iso -= rho*(abseta < 1.479 ? EAelB : EAelE);

        return iso;
    }

    static float getHcalIso(float isoval, float rho, float pt, float abseta, bool isMu) {
        float EAmuB = 0.022;
        float EAmuE = 0.030;
        float EAelB = 0.021;
        float EAelE = 0.040;

        float iso = isoval;

        if (isMu) iso -= rho*(abseta < 1.479 ? EAmuB : EAmuE);
        else iso -= rho*(abseta < 1.479 ? EAelB : EAelE);

        return iso;
    }

};

#endif
