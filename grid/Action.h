#ifndef ACTION_H
#define ACTION_H

#include "Fields/ScalarField.h"
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


template<typename GroupElement>
class GaugeAction {
public:
    explicit GaugeAction(double beta) : _beta(beta) {};

    template<int d>
    double S(const Lattice<VectorField<GroupElement, d>, d>& U) const {
        auto n = static_cast<double>(GroupElement::groupdim);
        return _beta * ( 1 - sum(plaqTraceSum(U), U.grid()).val() / n );
    }

    template<int d>
    auto Force(const Lattice<VectorField<GroupElement, d>, d>& U) const {





    }

private:
    double _beta;
};








#endif //ACTION_H