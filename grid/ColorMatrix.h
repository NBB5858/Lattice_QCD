#ifndef COLORMATRIX_H
#define COLORMATRIX_H


#include <Eigen/Dense>
#include <complex>

template<int Nc>
using CMat = Eigen::Matrix<std::complex<double>, Nc, Nc>;


template <class Derived>
inline auto adj(const Eigen::MatrixBase<Derived>& A)
    -> decltype(A.adjoint())
{
    return A.adjoint();
}

template<class T, std::enable_if_t<std::is_arithmetic_v<T> || std::is_same_v<T, std::complex<double>>,int> = 0>
inline T adj(T x) { return std::conj(x); }

template <class Derived>
inline auto trace(const Eigen::MatrixBase<Derived>& A)
    -> typename Derived::Scalar
{
    return A.trace();
}
template<class T, std::enable_if_t<std::is_arithmetic_v<T> || std::is_same_v<T, std::complex<double>>,int> = 0>
inline T trace(T x) { return x; }

template <class Derived>
inline typename Derived::PlainObject TA(const Eigen::MatrixBase<Derived>& M) {
    using Mat = typename Derived::PlainObject;
    const int Nc = static_cast<int>(M.rows());

    Mat AH = 0.5 * (M.derived() - M.adjoint());   // anti-Hermitian part
    if (Nc == 1) return AH;                       // U(1): no traceless projection

    const std::complex<double> tr = AH.trace();
    return AH - (tr / static_cast<double>(Nc)) * Mat::Identity();
}

inline CMat<1> U1exp(const CMat<1>& A) {
    CMat<1> out;
    out(0,0) = std::exp(A(0,0));
    return out;
}

template<int Nc>
inline CMat<Nc> Identity() {
    return CMat<Nc>::Identity();
}

template<int Nc>
inline CMat<Nc> Zero() {
    return CMat<Nc>::Zero();
}


#endif //COLORMATRIX_H
