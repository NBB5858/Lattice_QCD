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


// how do I even define a plaquette
// given a lattice, want to define all the plaquettes
// this is tricky because

// how would I define a single plaquette?


// we'll define the field U, which is defined at each site. it's a vector field, U_{\mu}
// U_{1} for example contains the value of U in the direction of the bond
// probably convenient to store U_{-1} as well, which is tricky because need to coordinate it's value
// with the adjacent


// U[ss](1)U[ss+mu](-0)U(-1)U(0)

//

// U(1)

// gauge fields are different because they
// have inverse
// is it an inverse for each direction? I think it is...

// at each point, for each plane, we define a plaquette

//

//
// template<typename GaugeField>
// GaugeField plaquette(int ss, int mu, int nu, const Lattice<VectorField<GaugeField>>& U, GridBase* grid ) {
//
//     int ss;
//     int right = grid->_stencil.neighbors[ss].p[mu];
//     int up = grid->_stencil.neighbors[ss].p[nu];
//     int right_up = grid->_stencil.neighbors[right].p[nu];
//
//     // each site stores just a field, not a field and it's inverse
//     U(ss)[mu] * U(right)[nu] * inverse(U(up)[mu]) * inverse(U(ss)[nu]);
//
// }

// how do I ensure this is parallelized?
// (1) I want to compute point by point, orientation by orientation plaquettes in parallel
// (3) this will give a lattice of plaquettes, which I need to trace over, also in parallel
// (4) once the traces are done, I need to sum over the traces

// will have something like
// GetPlaquettes(Lattice) --> Calculates all plaquettes site by site. Plaquette is the matrix product,
// no trace.  Returns... a lattice where at each site we have plaquettes for every plane
// so will return a tensor field U_{\mu \nu}, which I think is already fine, have this already.


// Vectorfield of vector field of group element

// should this be something which is evaluated lazily?
// I think so; but grid doesn't do this
// can I do this more simply?
// it's probably fine to do it that way
//

// template<typename GroupElement>
// VectorField<VectorField<GroupElement>> PlaquetteLattice(Lattice<VectorField<GroupElement>>& U) {
//
// }


// plaquette(U) would return the plaquetteified lattice;
//




//
//
// template<typename GaugeField>
// class GaugeAction {
// public:
//     GaugeAction(double beta) : _beta(beta) {};
//
//
//     double S(const Lattice<GaugeField>& U) const {
//
//     }
//
//
//
// private:
//     double _beta;
// };








#endif //ACTION_H
