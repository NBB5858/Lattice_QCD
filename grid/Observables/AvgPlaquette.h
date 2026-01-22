#ifndef AVGPLAQUETTE_H
#define AVGPLAQUETTE_H

#include <cstddef>
#include <vector>
#include <cmath>

#include "ObservableBase.h"
#include "../openmp/Reductions.h"

class AvgPlaquette : public ObservableBase {
public:
    explicit AvgPlaquette(int nobs) : ObservableBase(nobs) {}

    template<typename FieldType, int d>
    void measure(const Lattice<FieldType, d>& U) {
        const float Nd = static_cast<float>(U.grid()->dims().size());
        const float nsites = static_cast<float>(U.grid()->nsites());
        const float Nplaq = nsites * Nd * (Nd - 1.0f) / 2.0f;
        float N = static_cast<float>(FieldTraits<FieldType>::groupdim);

        float avg_plaq = (1.0f/N) * sum(plaqReTraceSum(U), U.grid()) / Nplaq;

        _cache.emplace_back(avg_plaq);
    }
};



#endif //AVGPLAQUETTE_H
