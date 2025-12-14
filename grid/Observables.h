#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "Reductions.h"

class Magnetization {
public:
    explicit Magnetization(int nobs) : _nobs(nobs) {_cache.reserve(_nobs);}

    template<typename FieldType>
    void measure(Lattice<FieldType>& U) {
        double mag = sum(U, U.grid()) / U.grid()->osites();
        _cache.push_back(mag);
    }

    double get_cache(int i) const {return _cache[i];}
    double mean() const {
        double sum = std::accumulate(_cache.begin(), _cache.end(), 0.0);
        return sum / static_cast<double>(_cache.size());
    }

private:
    int _nobs;
    std::vector<double> _cache;
};

#endif //OBSERVABLES_H
