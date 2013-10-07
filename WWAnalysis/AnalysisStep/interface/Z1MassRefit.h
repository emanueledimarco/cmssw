#ifndef WWAnalysis_AnalysisStep_Z1MassRefit_h
#define WWAnalysis_AnalysisStep_Z1MassRefit_h

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
struct CompositeCandMassResolution;

class Z1MassRefit {
    public:
        Z1MassRefit(const CompositeCandMassResolution & pscaleTool, double maxChi2In=10) ;
        ~Z1MassRefit();

        struct Result {
            bool   applied;
            double l1, l2;
            double cov[3];
            double chi2In, chi2Out;
        };
        
        Result refit(const reco::Candidate &z) ;

    private:
        double maxChi2In_;
        const CompositeCandMassResolution & pscaleTool_;
};
#endif
