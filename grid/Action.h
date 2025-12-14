#ifndef ACTION_H
#define ACTION_H

#include <cmath>

#include "ScalarField.h"
#include "Reductions.h"
#include "Lattice.h"

class ScalarAction {
public:

    ScalarAction(double mass, double lambda) : _mass(mass), _lambda(lambda) {};

    double S( const Lattice<ScalarField>& U) const {
        return sum(0.5*dot(grad(U), grad(U))
                      + 0.5 * pow(_mass, 2) * U * U
                      + (_lambda / 24.0) * U * U * U * U,
                    U.grid());
    }

    // auto Force( const Lattice<ScalarField>& U ) const  {
    //     return -1.0*box(U) + pow(_mass, 2) * U + (_lambda / 6.0) * U * U * U;
    // }

    auto Force( const Lattice<ScalarField>& U ) const  {
        return -1.0*box(U) + pow(_mass, 2) * U + (_lambda / 6.0) * U * U * U;
    }

    // if this just evaluates into a double, what is the poiint of doing lazy evaluaion
    // upon assignment into a lattice? when would that ever be used?
    // I suppose when you are updatng the lattice itself after an acceptance? but not really
    // grids force is also very differe, has p
    // ScalarField Force( const Lattice<ScalarField>& U ) const  {
    //     return -1.0*box(U) + pow(_mass, 2) * U + (_lambda / 6.0) * U * U * U;
    // }


    // Lattice<ScalarField>& generate_momentum(RNG rng) const {
    //
    // }

private:
    double _mass;
    double _lambda;
};





#endif //ACTION_H
