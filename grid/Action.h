#ifndef ACTION_H
#define ACTION_H

#include <cmath>

#include "ScalarField.h"
#include "Reductions.h"
#include "Lattice.h"

class ScalarAction {
public:

    ScalarAction(double mass2, double lambda) : _mass2(mass2), _lambda(lambda) {}

    double S( const Lattice<ScalarField>& U) const {
        return sum(0.5*dot(grad(U), grad(U))
                      + 0.5 * _mass2 * U * U
                      + (_lambda / 24.0) * U * U * U * U,
                    U.grid());
    }

    auto Force( const Lattice<ScalarField>& U ) const  {
        return -1.0*box(U) + _mass2 * U + (_lambda / 6.0) * U * U * U;
    }

private:
    double _mass2;
    double _lambda;
};





#endif //ACTION_H
