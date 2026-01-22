#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include "pragma_shorts.h"
#include "../Grid.h"
#include "../Lattice/Expression.h"

template<typename FieldType, int d>
inline auto sum(const Lattice<FieldType, d>& expr, const GridBase<d>* grid)
    -> decltype(simd::sum_lanes(eval(0, expr)))
{
    using R = decltype(eval(0, expr));
    R ret{};
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });

    return simd::sum_lanes(ret);
}

template<typename Op, typename T1, int d>
inline auto sum(const LatticeUnaryExpression<Op, T1>& expr, const GridBase<d>* grid)
    -> decltype(simd::sum_lanes(eval(0, expr)))
{
    using R = decltype(eval(0, expr));
    R ret{};
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return simd::sum_lanes(ret);
}

template<typename Op, typename T1, int d>
inline auto sum(const LatticeExtendedUnaryExpression<Op, T1>& expr, const GridBase<d>* grid)
    -> decltype(simd::sum_lanes(eval(0, expr)))
{
    using R = decltype(eval(0, expr));
    R ret{};
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return simd::sum_lanes(ret);
}

template<typename Op, typename T1, typename T2, int d>
inline auto sum(const LatticeBinaryExpression<Op, T1, T2>& expr, const GridBase<d>* grid)
    -> decltype(simd::sum_lanes(eval(0, expr)))
{
    using R = decltype(eval(0, expr));
    R ret{};
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return simd::sum_lanes(ret);
}


#endif //REDUCTIONS_H