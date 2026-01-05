#ifndef AVGPLAQUETTE_H
#define AVGPLAQUETTE_H

#include "ObservableBase.h"
#include "../Reductions.h"

class AvgPlaquette : public ObservableBase {
public:
    explicit AvgPlaquette(int nobs) : ObservableBase(nobs) {}

    template<typename FieldType, int d>
    void measure(const Lattice<FieldType, d>& U) {
        const double Nd = static_cast<double>(U.grid()->dims().size());
        const double os = static_cast<double>(U.grid()->osites());
        const double Nplaq = os * Nd * (Nd - 1.0) / 2.0;
        double N = static_cast<double>(FieldTraits<FieldType>::groupdim);

        double avg_plaq = (1.0/N) * sum(plaqReTraceSum(U), U.grid()) / Nplaq;

        _cache.emplace_back(avg_plaq);
    }
};



#endif //AVGPLAQUETTE_H
