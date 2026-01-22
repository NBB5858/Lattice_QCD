#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <vector>

#include "ObservableBase.h"
#include "../openmp/Reductions.h"

class Magnetization : public ObservableBase {
public:
    explicit Magnetization(int nobs) : ObservableBase(nobs) {}

    template<typename FieldType, int d>
    void measure(const Lattice<FieldType, d>& U) {
        auto nsites = static_cast<double>(U.grid()->nsites());
        double mag = sum(U, U.grid()) / nsites;
        _cache.emplace_back(mag);
    }
};

#endif //OBSERVABLES_H
