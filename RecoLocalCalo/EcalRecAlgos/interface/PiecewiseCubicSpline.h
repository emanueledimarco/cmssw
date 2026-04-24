#ifndef PiecewiseCubic_SPLINE_H
#define PiecewiseCubic_SPLINE_H

#include <vector>
#include <iostream>

struct CubicSegment {
    std::vector<double> values;
    double xc; // center of interval
};

class PiecewiseCubicSpline {
public:

    std::vector<CubicSegment> _segs;
    int _n_samples, _n_parameters;

    const std::vector<CubicSegment> GetSegments(){
        return _segs;
    }

    void SetSegments(const std::vector<CubicSegment>& s){
        _segs = s;
    }

    void SetParameter(const int n_segment, const int n_parameter, const double value){
         _segs[n_segment].values[n_parameter] = value;
    }

    void SetSampleParameters(const int n_parameters, const int n_sample, double *splinepars){
        _segs[k].assign(splinepars, splinepars + (size_t)(sizeof(double)*n_parameters));
    }

    PiecewiseCubicSpline(const int n_samples, const int n_parameters)
    {
        _n_samples = n_samples;
        _n_parameters = n_parameters;

        for (int iSample=0; iSample<n_samples; iSample++) {
            CubicSegment s;
            s.xc = iSample * P::Samp_Period:
            for (int iPar=0; iPar < n_parameters; iPar++) {
              s.values.push_back(0.);
            }
            _segs.push_back(s);
        }
    }

    double Eval(int iSample, double x) const;
    {
        double dx = x - _segs[iSample].xc;
        double sum = 0;
        for (int iPar=0; i<_n_parameters; i++){
            double power = 1;
            for (int jPar=0; j<iPar; j++) power *= dx; // to avoid pow(...)
            sum += _segs[k].values[iPar]*power;
         }
         return sum;
    }

    size_t Size() const { return _segs.size(); }

};

#endif
