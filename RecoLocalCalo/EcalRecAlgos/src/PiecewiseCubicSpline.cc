#include <RecoLocalCalo/EcalRecAlgos/interface/PiecewiseCubicSpline.h>

PiecewiseCubicSpline::PiecewiseCubicSpline(const char* file="coeffs_global.txt") {

    std::ifstream in(file);
    if(!in.is_open()){
        std::cerr << "Cannot open file: " << file << std::endl;
        return;
    }

    std::string line;
    std::getline(in, line); // skip header

    while(std::getline(in, line)){
        if(line.empty() || line[0]=='#') continue;

        std::istringstream ss(line);

        CubicSegment s;
        ss >> s.xc >> s.x0 >> s.x1
           >> s.values[0] >> s.values[1] >> s.values[2] >> s.values[3];

        if(ss.fail()) continue;

        _segs.push_back(s);
    }

    in.close();

};


double PiecewiseCubicSpline::Eval(double x) const {
    if(_segs.empty()) return 0.0;

    // find interval (linear search, same as your code)
    size_t k = 0;
    for(; k < _segs.size(); ++k){
        if(x >= _segs[k].x0 && x <= _segs[k].x1)
            break;
    }

    if(k == _segs.size()) k = _segs.size() - 1;
    double dx = x - seg.xc;
    return _segs[k].values[0] + _segs[k].values[1]*dx + _segs[k].values[2]*dx*dx + _segs[k].values[3]*dx*dx*dx;

};

