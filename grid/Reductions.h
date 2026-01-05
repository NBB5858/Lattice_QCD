#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include "pragma_shorts.h"
#include "Grid.h"
#include "Expression.h"


template<typename FieldType, int d>
inline auto sum(const Lattice<FieldType, d>& expr, const GridBase<d>* grid)
    -> decltype(eval(0, expr))
{
    using R = decltype(eval(0, expr));
    R ret{};  // 0 for double, (0,0) for complex<double>
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return ret;
}

template<typename Op, typename T1, int d>
inline auto sum(const LatticeUnaryExpression<Op, T1>& expr, const GridBase<d>* grid)
    -> decltype(eval(0, expr))
{
    using R = decltype(eval(0, expr));
    R ret{};
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return ret;
}

template<typename Op, typename T1, int d>
inline auto sum(const LatticeExtendedUnaryExpression<Op, T1>& expr, const GridBase<d>* grid)
    -> decltype(eval(0, expr))
{
    using R = decltype(eval(0, expr));
    R ret{};
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return ret;
}

template<typename Op, typename T1, typename T2, int d>
inline auto sum(const LatticeBinaryExpression<Op, T1, T2>& expr, const GridBase<d>* grid)
    -> decltype(eval(0, expr))
{
    using R = decltype(eval(0, expr));
    R ret{};
    thread_sum(ss, grid->osites(), ret, {
        ret += eval(ss, expr);
    });
    return ret;
}

//
// // For things like int dx V_{\mu} (x)...
// template<typename Op, typename T1, int d>
// auto sum( const LatticeUnaryExpression<Op, T1>& expr, const GridBase<d>* grid ) {
//     using IndexStructure = decltype(eval(0, expr));
//     IndexStructure ret{};
//
//     thread_sum(ss, grid->osites(), ret, {
//         ret += eval(ss, expr);
//     });
//     return ret;
// }
//
// template<typename Op, typename T1, int d>
// auto sum( const LatticeExtendedUnaryExpression<Op, T1>& expr, const GridBase<d>* grid ) {
//     using IndexStructure = decltype(eval(0, expr));
//     IndexStructure ret{};
//
//     thread_sum(ss, grid->osites(), ret, {
//         ret += eval(ss, expr);
//     });
//     return ret;
// }
//
// template<typename Op, typename T1, typename T2, int d>
// auto sum( const LatticeBinaryExpression<Op, T1, T2>& expr, const GridBase<d>* grid ) {
//     using IndexStructure = decltype(eval(0, expr));
//     IndexStructure ret{};
//
//     thread_sum(ss, grid->osites(), ret, {
//         ret += eval(ss, expr);
//     });
//     return ret;
// }


#endif //REDUCTIONS_H