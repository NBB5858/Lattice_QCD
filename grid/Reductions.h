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

// template<typename Op, typename T1>
// double sum( const LatticeExtendedUnaryExpression<Op, T1>& expr, const GridBase* grid ) {
//     double ret = 0.0;
//     thread_sum(ss, grid->osites(), ret, {
//         ret += eval(ss, expr);
//     });
//     return ret;
// }


template<typename Op, typename T1, typename T2, int d>
double sum( const LatticeBinaryExpression<Op, T1, T2>& expr, const GridBase<d>* grid ) {
    double ret = 0.0;
    thread_sum(ss, grid->osites(), ret, {
        // relies on knowing eval will return a scalar field
        ret += eval(ss, expr).val();
    });
    return ret;
}

#endif //REDUCTIONS_H
