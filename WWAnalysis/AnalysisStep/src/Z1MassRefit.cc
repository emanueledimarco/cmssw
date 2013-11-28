#include <cmath>
#include <cstdio>
#include <TMatrixTSym.h>
#include <TVectorT.h>
#include <TLorentzVector.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include "WWAnalysis/AnalysisStep/interface/Z1MassRefit.h"
#include "WWAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h"

/// Tools and private stuff
namespace {
    struct {
        TLorentzVector p4[3], p1s, p2s;
        double scaleErr[2];
    } inputs;

    double z1Gauss(const double *x) {
        inputs.p1s = inputs.p4[0]; inputs.p2s = inputs.p4[1];
        inputs.p1s *= (1.0 + x[0]);
        inputs.p2s *= (1.0 + x[1]);
        double term1 = x[0]/inputs.scaleErr[0];
        double term2 = x[1]/inputs.scaleErr[1];
        double term3 = ((inputs.p1s + inputs.p2s + inputs.p4[2]).M() - 91.2)/(2.5/2.33);
        return term1*term1 + term2*term2 + term3*term3;
    }
}

Z1MassRefit::Z1MassRefit(const CompositeCandMassResolution & pscaleTool, double maxChi2In) :
    maxChi2In_(maxChi2In),
    pscaleTool_(pscaleTool)
{
}

Z1MassRefit::~Z1MassRefit()
{
}

Z1MassRefit::Result
Z1MassRefit::refit(const reco::Candidate &z) 
{
    Result ret;

    const reco::Candidate &l1 = *z.daughter(0);
    const reco::Candidate &l2 = *z.daughter(1);
    inputs.p4[0] = TLorentzVector(l1.px(), l1.py(), l1.pz(), l1.energy());
    inputs.p4[1] = TLorentzVector(l2.px(), l2.py(), l2.pz(), l2.energy());
    inputs.scaleErr[0] = pscaleTool_.getPScaleError(l1);
    inputs.scaleErr[1] = pscaleTool_.getPScaleError(l2);
    if (z.numberOfDaughters() > 2) {
        const reco::Candidate &p1 = *z.daughter(2);
        inputs.p4[2] = TLorentzVector(p1.px(), p1.py(), p1.pz(), p1.energy());
    } else {
        inputs.p4[2] = TLorentzVector(0.,0.,0.,0.);
    }

    double zeros[2] =  { 0., 0. };
    ret.chi2In  =  z1Gauss(zeros);

    if (ret.chi2In > maxChi2In_) {
        ret.applied = false;
    } else {
        ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
        min.SetMaxFunctionCalls(1000000);
        min.SetMaxIterations(100000);
        min.SetTolerance(0.001);

        ROOT::Math::Functor f(&z1Gauss,2); 
        min.SetFunction(f);


        min.SetLimitedVariable(0,"x",0.,0.01,-8.,8.);
        min.SetLimitedVariable(1,"y",0.,0.01,-8.,8.);
        min.Minimize();

        const double *xs = min.X(); 
        ret.applied = true;
        ret.l1 = xs[0];
        ret.l2 = xs[1];
        ret.cov[0] = min.CovMatrix(0,0);
        ret.cov[1] = min.CovMatrix(0,1);
        ret.cov[2] = min.CovMatrix(1,1);
        ret.chi2Out =  z1Gauss(xs);
    }
    return ret;
}
