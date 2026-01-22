#ifndef EXPRESSION_H
#define EXPRESSION_H

#include "../simd/simd_types.h"
#include "../simd/CMatOverloads.h"
#include "../simd/LinkOverloads.h"

class LatticeBase;

template<typename FieldType, int d>
class Lattice;

class LatticeExpression{};

template<typename Op, typename T1>
class LatticeUnaryExpression;

template<typename Op, typename T1>
class LatticeExtendedUnaryExpression;

template<typename Op, typename T1, typename T2>
class LatticeBinaryExpression;

template <typename T> using is_lattice = std::is_base_of<LatticeBase, T>;
template <typename T> using is_lattice_expr = std::is_base_of<LatticeExpression, T>;

/********************
 * Unary operations
 ********************/
#define GridUnOpClass(name, ret)					\
struct name {								\
template<class _arg> static auto func(const _arg& a) -> decltype(ret) { return ret; } \
};

GridUnOpClass(UnarySquare, square(a));
GridUnOpClass(UnaryTrace, trace(a));
GridUnOpClass(UnaryAdj, adj(a));
GridUnOpClass(UnaryTA, TA(a));
GridUnOpClass(UnaryU1Exp, U1exp(a));
GridUnOpClass(UnarySU2Exp, SU2exp(a));
GridUnOpClass(UnaryGenericExp, genericExp(a));


#define GRID_UNOP(name)  name

#define GRID_DEF_UNOP(op, name)						\
template <typename T1, typename std::enable_if<is_lattice<T1>::value||is_lattice_expr<T1>::value,T1>::type * = nullptr> \
inline auto op(const T1 &arg) ->decltype(LatticeUnaryExpression<GRID_UNOP(name),T1>(GRID_UNOP(name)(), arg)) \
{ return LatticeUnaryExpression<GRID_UNOP(name),T1>(GRID_UNOP(name)(), arg); }

GRID_DEF_UNOP(square, UnarySquare);
GRID_DEF_UNOP(trace, UnaryTrace);
GRID_DEF_UNOP(adj, UnaryAdj);
GRID_DEF_UNOP(TA, UnaryTA);
GRID_DEF_UNOP(U1exp, UnaryU1Exp);
GRID_DEF_UNOP(SU2exp, UnarySU2Exp);
GRID_DEF_UNOP(genericExp, UnaryGenericExp);

/********************
 * Extended Unary operations
 ********************/
#define GridExtUnOpClass(name, ret)					\
struct name {								\
template<class _arg> static auto func(const uint64_t ss, const _arg& a) -> decltype(ret) { return ret; } \
};

GridExtUnOpClass(UnaryGradSquared, grad_squared(ss, a));
GridExtUnOpClass(UnaryBox, box(ss, a));
GridExtUnOpClass(UnaryPlaqReTraceSum, PlaqReTraceSum(ss, a));
GridExtUnOpClass(UnaryStapleSum, StapleSum(ss, a));

#define GRID_EXT_UNOP(name)  name

#define GRID_DEF_EXTENDED_UNOP(op, name)						\
template <typename T1, typename std::enable_if<is_lattice<T1>::value||is_lattice_expr<T1>::value,T1>::type * = nullptr> \
inline auto op(const T1 &arg) ->decltype(LatticeExtendedUnaryExpression<GRID_EXT_UNOP(name),T1>(GRID_EXT_UNOP(name)(), arg)) \
{ return LatticeExtendedUnaryExpression<GRID_EXT_UNOP(name),T1>(GRID_EXT_UNOP(name)(), arg); }

GRID_DEF_EXTENDED_UNOP(grad_squared, UnaryGradSquared);
GRID_DEF_EXTENDED_UNOP(box, UnaryBox);
GRID_DEF_EXTENDED_UNOP(plaqReTraceSum, UnaryPlaqReTraceSum);
GRID_DEF_EXTENDED_UNOP(stapleSum, UnaryStapleSum);

/********************
 * Binary operations
 ********************/
#define GridBinOpClass(name, combination) \
struct name { \
template <class _left, class _right> \
static auto func(const _left &lhs, const _right &rhs) -> decltype(combination) \
{ return combination; } \
};

GridBinOpClass(BinaryAdd, lhs +rhs);
GridBinOpClass(BinarySub, lhs -rhs);
GridBinOpClass(BinaryMul, lhs *rhs);
GridBinOpClass(BinaryDiv, lhs /rhs);
GridBinOpClass(BinaryAnd, lhs &rhs);
GridBinOpClass(BinaryOr, lhs | rhs);
GridBinOpClass(BinaryAndAnd, lhs &&rhs);
GridBinOpClass(BinaryOrOr, lhs || rhs);
GridBinOpClass(BinaryDot, Dot(lhs, rhs));

#define GRID_BINOP(name)  name

#define GRID_BINOP_LEFT(op, name) \
template <typename T1, typename T2, \
typename std::enable_if<is_lattice<T1>::value || \
                        is_lattice_expr<T1>::value, T1>::type * = nullptr > \
inline auto op( const T1 &lhs, const T2 &rhs ) \
-> decltype(LatticeBinaryExpression<GRID_BINOP(name), T1, T2>(GRID_BINOP(name)(), lhs, rhs)) \
{ return     LatticeBinaryExpression<GRID_BINOP(name),T1,T2>(GRID_BINOP(name)(),lhs,rhs); }

GRID_BINOP_LEFT(operator+, BinaryAdd);
GRID_BINOP_LEFT(operator-, BinarySub);
GRID_BINOP_LEFT(operator*, BinaryMul);
GRID_BINOP_LEFT(operator/, BinaryDiv);
GRID_BINOP_LEFT(operator&, BinaryAnd);
GRID_BINOP_LEFT(operator|, BinaryOr);
GRID_BINOP_LEFT(operator&&, BinaryAndAnd);
GRID_BINOP_LEFT(operator||, BinaryOrOr);
GRID_BINOP_LEFT(dot, BinaryDot);


#define GRID_BINOP_RIGHT(op, name)					\
template <typename T1, typename T2,					\
typename std::enable_if<!is_lattice<T1>::value&&!is_lattice_expr<T1>::value,T1>::type * = nullptr, \
typename std::enable_if< is_lattice<T2>::value|| is_lattice_expr<T2>::value,T2>::type * = nullptr> \
inline auto op(const T1 &lhs, const T2 &rhs)				\
->decltype(LatticeBinaryExpression<GRID_BINOP(name),T1,T2>(GRID_BINOP(name)(),lhs, rhs)) \
{ return     LatticeBinaryExpression<GRID_BINOP(name),T1,T2>(GRID_BINOP(name)(),lhs, rhs); }

GRID_BINOP_RIGHT(operator+, BinaryAdd);
GRID_BINOP_RIGHT(operator-, BinarySub);
GRID_BINOP_RIGHT(operator*, BinaryMul);
GRID_BINOP_RIGHT(operator/, BinaryDiv);
GRID_BINOP_RIGHT(operator&, BinaryAnd);
GRID_BINOP_RIGHT(operator|, BinaryOr);
GRID_BINOP_RIGHT(operator&&, BinaryAndAnd);
GRID_BINOP_RIGHT(operator||, BinaryOrOr);
GRID_BINOP_RIGHT(dot, BinaryDot);

/********************
 * Evaluation
 ********************/
template <class FieldType, int d>
auto eval(const uint64_t ss, const Lattice<FieldType, d> &arg){ return arg(ss); }

template <typename T>
auto eval(const uint64_t /*ss*/, T arg) { return arg; }

template <typename Op, typename T1>
auto eval(const uint64_t ss, const LatticeUnaryExpression<Op, T1> &expr) {
    return expr.op.func(eval(ss, expr.arg));
}

template <typename Op, typename T1>
auto eval(const uint64_t ss, const LatticeExtendedUnaryExpression<Op, T1> &expr){
    return expr.op.func(ss, expr.arg);
}

template <typename Op, typename T1, typename T2>
auto eval(const uint64_t ss, const LatticeBinaryExpression<Op, T1, T2> &expr){
    return expr.op.func(eval(ss, expr.arg1), eval(ss, expr.arg2));
}

/********************
 * Expressions
 ********************/
// arithmetic types are not persistent and should be stored as values
// should stored references to lattice and lattice expressions
template<typename T>
using  wrap = std::conditional_t< std::is_arithmetic_v<T> || is_lattice_expr<T>::value, T, const T& >;


template <typename Op, typename T1>
class LatticeUnaryExpression : public LatticeExpression {
public:
    Op op;
    wrap<T1> arg;
    LatticeUnaryExpression(Op op, const T1 &arg) : op(op), arg(arg) {}
};

template <typename Op, typename T1>
class LatticeExtendedUnaryExpression : public LatticeExpression {
public:
    Op op;
    wrap<T1> arg;
    LatticeExtendedUnaryExpression(Op op, const T1 &arg) : op(op), arg(arg) {}
};

template <typename Op, typename T1, typename T2>
class LatticeBinaryExpression : public LatticeExpression {
public:
    Op op;
    wrap<T1> arg1;
    wrap<T2> arg2;
    LatticeBinaryExpression(Op op, const T1 &arg1, const T2 &arg2) : op(op), arg1(arg1), arg2(arg2) {}
};

#endif //EXPRESSION_H
