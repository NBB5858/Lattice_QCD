#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <vector>

#include "ObservableBase.h"
#include "../Reductions.h"

class Magnetization : public ObservableBase {
public:
    explicit Magnetization(int nobs) : ObservableBase(nobs) {}

    template<typename FieldType, int d>
    void measure(const Lattice<FieldType, d>& U) {
        auto osites = static_cast<double>(U.grid()->osites());
        double mag = sum(U, U.grid()) / osites;
        _cache.emplace_back(mag);
    }
};

#endif //OBSERVABLES_H
