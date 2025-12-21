#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <vector>
#include <numeric>
#include <cmath>

#include "Reductions.h"


class Magnetization {
public:
    explicit Magnetization(int nobs) : _nobs(nobs) {_cache.reserve(_nobs);}

    template<typename FieldType, int d>
    void measure(const Lattice<FieldType, d>& U) {
        // divide by spatial lattice sites and euclidean times sites
        double mag = sum(U, U.grid()) / U.grid()->osites();
        _cache.emplace_back(mag);
    }

    double get_cache(int i) const {return _cache[i];}
    double mean(double burnin) const {
        if (_cache.empty()) return 0.0;

        int start = std::ceil(burnin*_nobs);
        double sum = 0;

        for(int i = start; i < _nobs; ++i) sum += _cache[i];
        return sum / static_cast<double>(_cache.size());
    }

    double abs_mean(double burnin) const {
        if(_cache.empty()) return 0.0;

        int start = std::ceil(burnin*_nobs);
        double sum = 0;

        for(int i = start; i < _nobs; ++i) sum += std::abs(_cache[i]);
        return sum / static_cast<double>(_cache.size());
    }

    std::vector<double> _cache;
private:
    int _nobs;
};

#endif //OBSERVABLES_H
