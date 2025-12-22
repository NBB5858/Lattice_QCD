#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include "pragma_shorts.h"
#include "Grid.h"
#include "Expression.h"

// functions for integrating densities
template<int d>
inline double sum( const Lattice<ScalarField, d>& expr, const GridBase<d>* grid ) {
    double ret = 0.0;
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr).val();
    });
    return ret;
}

template<typename Op, typename T1, int d>
double sum( const LatticeUnaryExpression<Op, T1>& expr, const GridBase<d>* grid ) {
    double ret = 0.0;
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr).val();
    });
    return ret;
}

template<typename Op, typename T1, int d>
double sum( const LatticeExtendedUnaryExpression<Op, T1>& expr, const GridBase<d>* grid ) {
    double ret = 0.0;
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return ret;
}


template<typename Op, typename T1, typename T2, int d>
double sum( const LatticeBinaryExpression<Op, T1, T2>& expr, const GridBase<d>* grid ) {
    double ret = 0.0;
    thread_sum(ss, grid->osites(), ret, {
        // relies on knowing eval will return a scalar field?
        // is currently OKish because you only call sum on things that are scalars anyway
        // may become a problem later
        ret += eval(ss, expr).val();
    });
    return ret;
}


// what does summing over plaquettes actually mean?
// each site has a tensor like U_mn
// sum over each site and add component by component

// Lattice<VectorField<ScalarField>> V;
//
// sum(V);

// ret should be the same field type that we put in
// for a





#endif //REDUCTIONS_H
