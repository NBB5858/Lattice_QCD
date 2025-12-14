#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include "pragma_shorts.h"
#include "Grid.h"
#include "Expression.h"

// functions for integrating densities
inline double sum( const Lattice<ScalarField>& expr, const GridBase* grid ) {
    double ret = 0.0;
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr).val();
    });
    return ret;
}

template<typename Op, typename T1>
double sum( const LatticeUnaryExpression<Op, T1>& expr, const GridBase* grid ) {
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


template<typename Op, typename T1, typename T2>
double sum( const LatticeBinaryExpression<Op, T1, T2>& expr, const GridBase* grid ) {
    double ret = 0.0;
    thread_sum(ss, grid->osites(), ret, {
        // relies on knowing eval will return a scalar field
        ret += eval(ss, expr).val();
    });
    return ret;
}

#endif //REDUCTIONS_H
