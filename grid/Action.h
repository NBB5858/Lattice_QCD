#ifndef ACTION_H
#define ACTION_H

#include "Fields/ScalarField.h"
#include "Fields/U1Field.h"
#include "Reductions.h"
#include "Lattice.h"

class ScalarAction {
public:

    ScalarAction(double mass2, double lambda) : _mass2(mass2), _lambda(lambda) {}

    template<int d>
    double S(const Lattice<ScalarField, d>& U) const {
        return sum(0.5*dot(grad(U), grad(U))
                      + 0.5 * _mass2 * U * U
                      + (_lambda / 24.0) * U * U * U * U,
                    U.grid());
    }

    template<int d>
    auto Force(const Lattice<ScalarField, d>& U) const  {
        return -1.0*box(U) + _mass2 * U + (_lambda / 6.0) * U * U * U;
    }

private:
    double _mass2;
    double _lambda;
};


template<typename GaugeField>
class WilsonAction {
public:
    explicit WilsonAction(double beta) : _beta(beta) {};

    template<int d>
    double S(const Lattice<GaugeField, d>& U) const {
        const double Nd = static_cast<double>(U.grid()->dims().size());
        const double os = static_cast<double>(U.grid()->osites());
        const double Nplaq = os * Nd * (Nd - 1.0) / 2.0;
        return _beta * ( Nplaq - (1.0/N) * sum(plaqReTraceSum(U), U.grid()) );
    }

    template<int d>
    auto Force(const Lattice<GaugeField, d>& U) const {
        return (_beta / N) * TA( U * adj(stapleSum(U)) );
    }

private:
    double _beta;
    double N = static_cast<double>(FieldTraits<GaugeField>::groupdim);
};








#endif //ACTION_H